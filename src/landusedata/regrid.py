import tempfile
import os
import numpy as np
import xesmf as xe
import xarray as xr
# Keep for debugging if needed:
# import sys
# sys.exit(4)


# This function determines if we are using a SCRIP
# or a mesh file to create weights
def _detect_se_grid_type(filepath):
    ds = xr.open_dataset(filepath)
    if 'nodeCoords' in ds and 'elementConn' in ds:
        return 'mesh'
    elif 'grid_center_lat' in ds and 'grid_corner_lat' in ds:
        return 'scrip'
    else:
        raise ValueError(
            f"Unrecognized grid file format: {filepath}. "
            f"Expected ESMF mesh (nodeCoords, elementConn) or "
            f"SCRIP (grid_center_lat, grid_corner_lat)."
        )

def _scrip_to_mesh_tempfile(scrip_file):
    ds = xr.open_dataset(scrip_file)
    ncells   = ds.dims['grid_size']
    ncorners = ds.dims['grid_corners']

    # Treat each corner as a unique node — avoids deduplication complexity
    # Original size is (grid_size, grid_corners=4)
    # Flattening makes each cell "i" have coordinates:
    # [i*4+1, i*4+2, i*4+3, i*4+4]
    node_lons = ds['grid_corner_lon'].values.reshape(-1)
    node_lats = ds['grid_corner_lat'].values.reshape(-1)
    nNodes = len(node_lons)

    # 1-indexed connectivity: cell i owns nodes [i*ncorners+1 .. (i+1)*ncorners]
    conn = (np.arange(nNodes, dtype=np.int32) + 1).reshape(ncells, ncorners)

    mesh_ds = xr.Dataset({
        'nodeCoords': xr.DataArray(
            np.stack([node_lons, node_lats], axis=1),
            dims=['nodeCount', 'coordDim'],
        ),
        'elementConn': xr.DataArray(
            conn,
            dims=['elementCount', 'maxNodePElement'],
        ),
        'numElementConn': xr.DataArray(
            np.full(ncells, ncorners, dtype=np.int8),
            dims=['elementCount'],
        ),
        'centerCoords': xr.DataArray(
            np.stack([ds['grid_center_lon'].values,
                      ds['grid_center_lat'].values], axis=1),
            dims=['elementCount', 'coordDim'],
        ),
    })

    with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as f:
        tmp_path = f.name
    mesh_ds.to_netcdf(tmp_path)
    return tmp_path

    

def RegridConservative(ds_to_regrid, ds_regrid_target, regridder_weights_file, regrid_reuse, \
                       regrid_method="conservative", se_grid_file=None):

    # define the regridder transformation
    regridder = GenerateRegridder(ds_to_regrid, ds_regrid_target, regridder_weights_file, \
                                  regrid_reuse, regrid_method=regrid_method, \
                                  se_grid_file=se_grid_file)

    # Loop through the variables to regrid
    ds_regrid = RegridLoop(ds_to_regrid, regridder)

    return (ds_regrid, regridder)

def GenerateRegridder(ds_to_regrid, ds_regrid_target, regridder_weights_file,
                      regrid_reuse, regrid_method="conservative", se_grid_file=None):

    # This is a little confusing, but some things to help clear it up
    #
    # "regridder" is an esmf method that needs to be built, this is the thing
    # that we use to transfer the actual data arrays from one grid to another
    #
    # The mapping that regridder will use, which has weights and stuff to go
    # from one grid to another, can be saved in a file so that its ready to use
    # at later points in the call sequence, rebuilding this mapping takes time. we
    # save that in a netcdf file (yes, we don't save in a data array because
    # we want to use the esmf libraries, which want an nc file)
    #
    # For equal angle grids, creating this mapping and regridder is straight forward
    # and just needs grid cell centers. However, if its not (like spectral element),
    # then we need to use a file that defines the grid cell corners too. This
    # is the mesh or scrip file.
    
    print("\nDefining regridder, method: ", regrid_method)

    # This checks to make sure that if we have a non equal angle grid
    # that we must have a mesh/script file
    if se_grid_file is None and 'lat' not in ds_regrid_target.dims:
        raise ValueError(
            "Target grid has no 'lat' dimension — it appears to be a spectral element"
            "or unstructured grid. Provide --se_grid_file with a SCRIP or ESMF mesh "
            "file describing the target grid."
        )
    
    if se_grid_file is not None:
        grid_type = _detect_se_grid_type(se_grid_file)
        tmp_mesh = None
        try:
            mesh_path = se_grid_file if grid_type == 'mesh' \
                        else _scrip_to_mesh_tempfile(se_grid_file)
            if grid_type != 'mesh':
                tmp_mesh = mesh_path

            if regrid_reuse:
                regridder = xe.Regridder(ds_to_regrid, mesh_path, regrid_method,
                                         weights=regridder_weights_file)
            else:
                regridder = xe.Regridder(ds_to_regrid, mesh_path, regrid_method)
                regridder.to_netcdf(regridder_weights_file)
                print("regridder saved to file: ", regridder_weights_file)
        finally:
            if tmp_mesh is not None:
                os.unlink(tmp_mesh)

    elif regrid_reuse:
        regridder = xe.Regridder(ds_to_regrid, ds_regrid_target, regrid_method,
                                 weights=regridder_weights_file)
    else:
        regridder = xe.Regridder(ds_to_regrid, ds_regrid_target, regrid_method)
        regridder.to_netcdf(regridder_weights_file)
        print("regridder saved to file: ", regridder_weights_file)

    return regridder


def RegridLoop(ds_to_regrid, regridder):

    # To Do: implement this with dask
    print("\nRegridding")

    # Loop through the variables one at a time to conserve memory
    ds_varnames = list(ds_to_regrid.variables.keys())
    varlen = len(ds_to_regrid.variables)
    first_var = True
    for i, var in enumerate(ds_to_regrid):
        msg_text = " variable {}/{}: {}\n".format(i+1, varlen, var)

        # Skip time variable and only regrid variables that match the lat/lon shape.
        do_regrid = not "time" in var
        do_regrid = do_regrid and \
            ds_to_regrid[var][0].shape == (ds_to_regrid.lat.shape[0], ds_to_regrid.lon.shape[0])
        if not do_regrid :
            print("skipping" + msg_text)
            continue

        print("regridding" + msg_text)

        # For the first non-coordinate variable, copy and regrid the dataset as a whole.
        # This makes sure to correctly include the lat/lon in the regridding.
        if first_var:
            ds_regrid = ds_to_regrid[var].to_dataset() # convert data array to dataset
            ds_regrid = regridder(ds_regrid)
            first_var = False

        # Once the first variable has been included, then we can regrid by variable
        else:
            # Somehow regridder() can affect its input in-place! We don't want this,
            # so we pass it in as a copy.
            ds_regrid[var] = regridder(ds_to_regrid[var].copy())

    return(ds_regrid)
