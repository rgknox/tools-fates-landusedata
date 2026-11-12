import xarray as xr
from collections import Counter

COMMON_CAP_THRESHOLD = 1.0e-15


def ImportRegridTarget(filename):
    dataset = xr.open_dataset(filename)

    # Check the file type
    dim_list = list(dataset.dims)
    if all(item not in dim_list for item in ['lsmlat','gridcell']):
        raise TypeError("incorrect file, must be surface dataset")

    # Prepare the the regrid dataset
    dataset = _RegridTargetPrep(dataset)

    return dataset


def _RegridTargetPrep(regrid_target):
    """
    Rename the CLM surface data file regrid target lat/lon dimensions

    The CLM surface datasets use a dimensions name that is not recognizable
    by esmf.  As such, this function renames the dimensions.  It also adds
    lat/lon coordinates based on the LONGXY and LATIXY variables.
    """

    # Skip if lat and lon are already the dimension names
    if "lat" in regrid_target.dims and "lon" in regrid_target.dims:
        return regrid_target

    # Drop dimension variables if they already exist
    regrid_target = regrid_target.drop_vars("lat", errors="ignore")
    regrid_target = regrid_target.drop_vars("lon", errors="ignore")

    # Rename dimensions and add coordinates
    if 'lsmlat' in regrid_target.dims:
        regrid_target = regrid_target.rename({'lsmlat':'lat','lsmlon':'lon'})
        regrid_target['lon'] = regrid_target.LONGXY.isel(lat=0)
        regrid_target['lat'] = regrid_target.LATIXY.isel(lon=0)
    elif 'gridcell' in regrid_target.dims:
        regrid_target['lon'] = regrid_target.LONGXY
        regrid_target['lat'] = regrid_target.LATIXY

    return regrid_target


# Open the LUH2 static data file
def ImportLUH2StaticFile(filename):
    dataset = xr.open_dataset(filename)
    if "lat_bnds" in dataset.var():
        dataset = dataset.rename({"lat_bnds": "lat_bounds", "lon_bnds": "lon_bounds"})
    # Check to see if the imported dataset has correct variables

    listcheck = ['ptbio', 'fstnf', 'carea', 'icwtr', 'ccode', 'lat_bounds', 'lon_bounds']
    if Counter(list(dataset.var())) != Counter(listcheck):
        raise TypeError("incorrect file, must be LUH2 static file")

    # dataset = dataset.drop(labels=['lat_bnds','lon_bnds'])

    # Convert all data from single to double precision
    dataset = dataset.astype("float64")
    return dataset


# Define the land/ocean mask based on the ice/water data
# from the LUH2 static data set
def DefineStaticMask(dataset):
    try:
        mask = xr.where(dataset.icwtr != 1, 1, 0)
        # mask = (1.-dataset.icwtr) / (1.-dataset.icwtr)
    except AttributeError:
        raise AttributeError("incorrect dataset, must be static luh2 dataset")
    else:
        return mask


# Surface dataset specific masking sub-function
def SetMaskRegridTarget(dataset):
    try:
        dataset["mask"] = xr.where(dataset.PCT_NATVEG > 0, 1, 0)
    except AttributeError:
        raise AttributeError("incorrect dataset, must be CLM surface dataset")
    else:
        return dataset


# TODO: add the follow common functions
# - write to files
# - read dataset into list of datasets
