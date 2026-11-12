import numpy as np
import xarray as xr
import esmpy as esmf

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

def _build_latlon_grid(ds,src_mask=None):

    # build an esmf "grid" fiel type from lat/lon data
    # since these grids are equal angle, we can do simple math to retrieve the corners
    lat = ds['lat'].values
    lon = ds['lon'].values
    nlat, nlon = len(lat), len(lon)

    # Compute cell borders from centers (equal-angle assumption)
    # we don't need 2d fields for this. 
    dlat = lat[1] - lat[0]
    dlon = lon[1] - lon[0]
    lat_b = np.concatenate([[lat[0] - dlat/2], (lat[:-1] + lat[1:])/2, [lat[-1] + dlat/2]])
    lon_b = np.concatenate([[lon[0] - dlon/2], (lon[:-1] + lon[1:])/2, [lon[-1] + dlon/2]])

    grid = esmf.Grid(
        np.array([nlon, nlat], dtype=np.int32),
        num_peri_dims=1,
        coord_sys=esmf.CoordSys.SPH_DEG )

    grid.add_coords(staggerloc=esmf.StaggerLoc.CENTER)
    grid.add_coords(staggerloc=esmf.StaggerLoc.CORNER)
    
    # ESMF is lon-major so arrays are (nlon, nlat) — transpose from numpy convention
    lon_2d, lat_2d = np.meshgrid(lon, lat)
    grid.get_coords(0,staggerloc=esmf.StaggerLoc.CENTER)[...] = lon_2d.T
    grid.get_coords(1,staggerloc=esmf.StaggerLoc.CENTER)[...] = lat_2d.T

    lon_b_2d, lat_b_2d = np.meshgrid(lon_b[:-1], lat_b)
    grid.get_coords(0, staggerloc=esmf.StaggerLoc.CORNER)[...] = lon_b_2d.T
    grid.get_coords(1, staggerloc=esmf.StaggerLoc.CORNER)[...] = lat_b_2d.T


    if src_mask is not None:
        grid.add_item(esmf.GridItem.MASK, staggerloc=esmf.StaggerLoc.CENTER)
        grid.get_item(esmf.GridItem.MASK, staggerloc=esmf.StaggerLoc.CENTER)[...] = \
            src_mask.astype(np.int32).T
    
    return grid


def _build_target_mesh(se_grid_file):
    grid_type = _detect_se_grid_type(se_grid_file)
    if grid_type == 'scrip':
        return esmf.Mesh(filename=se_grid_file,
                         filetype=esmf.FileFormat.SCRIP,
                         convert_to_dual=False)
    else:
        return esmf.Mesh(filename=se_grid_file,
                         filetype=esmf.FileFormat.ESMFMESH,
                         convert_to_dual=False)

class BuildRegridder:

    def __init__(self, src_grid, dst, weights_file=None, reuse_weights=False,
                 target_lat=None, target_lon=None, mesh_lat=None, mesh_lon=None,
                 use_src_mask=False, src_mask=None):
        
        self._dst_is_mesh = isinstance(dst, esmf.Mesh)

        self._src_field = esmf.Field(src_grid, staggerloc=esmf.StaggerLoc.CENTER)
        if self._dst_is_mesh:
            self._dst_field = esmf.Field(dst, meshloc=esmf.MeshLoc.ELEMENT)
        else:
            self._dst_field = esmf.Field(dst, staggerloc=esmf.StaggerLoc.CENTER)

        if reuse_weights and weights_file is not None:
            self._regrid = esmf.RegridFromFile(self._src_field, self._dst_field, weights_file)
        else:
            self._regrid = esmf.Regrid(
                self._src_field, self._dst_field,
                regrid_method=esmf.RegridMethod.CONSERVE,
                norm_type=esmf.NormType.DSTAREA,
                unmapped_action=esmf.UnmappedAction.IGNORE,
                filename=weights_file ,
                src_mask_values=np.array([0], dtype=np.int32) if use_src_mask else None )

        self._src_mask = src_mask
        self._target_lat = target_lat
        self._target_lon = target_lon
        self._mesh_lat = mesh_lat
        self._mesh_lon = mesh_lon
        
    def _apply(self, data_2d):
        self._src_field.data[...] = data_2d.T  # numpy (lat,lon) → esmf (lon,lat)
        self._regrid(self._src_field, self._dst_field)
        return self._dst_field.data.copy()

    def _regrid_variable(self, da):
        other_dims = [d for d in da.dims if d not in ('lat', 'lon')]
        da = da.transpose(*other_dims, 'lat', 'lon')
        data = da.values
        extra_dims = data.shape[:-2]
        extra_dim_names = list(da.dims[:-2])
        spatial_dims = ['ncells'] if self._dst_is_mesh else ['lat', 'lon']

        if self._dst_is_mesh and self._mesh_lat is not None:
            coords = {'lat': ('ncells', self._mesh_lat), 'lon': ('ncells', self._mesh_lon)}
        elif not self._dst_is_mesh and self._target_lat is not None:
            coords = {'lat': self._target_lat, 'lon': self._target_lon}
        else:
            coords = {}
            
        if extra_dims:
            out = np.empty(extra_dims + self._dst_field.data.shape)
            for idx in np.ndindex(extra_dims):
                out[idx] = self._apply(data[idx])
            return xr.DataArray(out, dims=extra_dim_names + spatial_dims, coords=coords)
        else:
            return xr.DataArray(self._apply(data), dims=spatial_dims, coords=coords)


    def __call__(self, ds):
        if isinstance(ds, xr.DataArray):
            return self._regrid_variable(ds)
        result = {var: self._regrid_variable(ds[var]) for var in ds.data_vars}
        return xr.Dataset(result)

    
    
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
                      regrid_reuse, regrid_method="conservative", se_grid_file=None,
                      src_mask=None):

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
            "Target grid has no 'lat' dimension — it appears to be a spectral element "
            "or unstructured grid. Provide --se_grid_file with a SCRIP or ESMF mesh "
            "file describing the target grid."
        )

    # construct an ESMF "grid file" from source coordinates
    src_grid = _build_latlon_grid(ds_to_regrid,src_mask=src_mask)

    # construct an ESMF "mesh" or "grid" from the target coordinates
    if se_grid_file is not None:
        target_meshorgrid = _build_target_mesh(se_grid_file)
        mesh_lon = target_meshorgrid.get_coords(0, meshloc=esmf.MeshLoc.ELEMENT)
        mesh_lat = target_meshorgrid.get_coords(1, meshloc=esmf.MeshLoc.ELEMENT)
        target_lat = None
        target_lon = None
        
    else:
        target_meshorgrid = _build_latlon_grid(ds_regrid_target)
        mesh_lat = None
        mesh_lon = None
        target_lat = ds_regrid_target['lat'].values
        target_lon = ds_regrid_target['lon'].values
        
    regridder = BuildRegridder(src_grid, target_meshorgrid,
                               weights_file=regridder_weights_file,
                               reuse_weights=regrid_reuse,
                               target_lat=target_lat,
                               target_lon=target_lon,
                               mesh_lat=mesh_lat,
                               mesh_lon=mesh_lon,
                               use_src_mask=(src_mask is not None),
                               src_mask=src_mask)

    if not regrid_reuse:
        print("regridder saved to file: ", regridder_weights_file)

    return regridder



def RegridLoop(ds_to_regrid, regridder):

    # To Do: implement this with dask
    print("\nRegridding")

    # Loop through the variables one at a time to conserve memory
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
