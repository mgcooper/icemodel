# icemodel.netcdf

Purpose: Write icemodel output to NetCDF files and read it back, with CF
metadata for each channel.

Contents:

- File writing:
  - `makencfile` writes one NetCDF file for each simulation year from the
    saved `ice1` or `ice2` output files of a run. It does not support `met`:
    `getvarinfo` and `getchunksize` have no `met` branch, so a `met` call
    stops with an error.
  - `create` creates the file with its global attributes.
  - `defdimid` and `defdimvars` define the grid and time dimensions.
  - `defdatavars` defines the data variables and their attributes.
  - `writedims` writes the dimension values.
  - `writeice1` and `writeice2` write the `ice1` and `ice2` data.
  - `setfilename` builds the output file name.
  - `config` sets the NetCDF API preferences.
- Sizing:
  - `getdimdata` returns the depth and time dimensions of a run.
  - `getdimsize` returns the size of each named dimension.
  - `maxcells` returns the largest grid-cell count that fits a memory limit.
  - `getchunksize` returns the chunk sizes for the data variables.
  - `getvarinfo` loads one data file to get variable shapes and names.
- Reading and editing:
  - `ncread` reads a file into memory.
  - `getvardata` reads the data into arrays.
  - `nctype2mat` maps NetCDF data types to MATLAB types.
  - `redefatt` rewrites the `units` attribute of one variable in each file
    of a list.
- Metadata:
  - `getdefaults` returns these metadata lists for one file type:
    - variable names
    - standard names
    - long names
    - units
    - axes
  - `+defaults` holds the metadata lists and maps:
    - `varnames`
    - `standardnames`
    - `longnames`
    - `units`
    - `axes`
    - `variables` returns the metadata map for every icemodel channel.
      `icemodel.forcing`, `icemodel.plot`, and `icemodel.verification` read
      channel units and names from it.
    - `variable` returns one entry of that map.
    - `cfStandardNames` loads the CF Standard Name Table.
