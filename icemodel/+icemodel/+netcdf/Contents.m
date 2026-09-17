% +NETCDF
%
%   Contents file for +NETCDF and its subfolders.
%
%   +NETCDF
%   icemodel.netcdf.config                   - Configure icemodel.netcdf API preferences
%   icemodel.netcdf.create                   - Create a new NetCDF file with the given properties and global
%   icemodel.netcdf.defdatavars              - Define the icemodel data variables and attributes
%   icemodel.netcdf.defdimid                 - Define icemodel netcdf file dimensions
%   icemodel.netcdf.defdimvars               - Define icemodel netcdf grid and time dims and attributes
%   icemodel.netcdf.getchunksize
%   icemodel.netcdf.getdefaults
%   icemodel.netcdf.getdimdata               - Get dimensions of icemodel simulation data
%   icemodel.netcdf.getdimsize               - Return the size of each named dimension
%   icemodel.netcdf.getvardata               - Read icemodel data into memory and fill the arrays
%   icemodel.netcdf.getvarinfo               - Load one data file to get the shape and variable names
%   icemodel.netcdf.makencfile
%   icemodel.netcdf.maxcells                 - Calculate the maximum number of gridcells for icemodel nc file
%   icemodel.netcdf.ncread                   - Read icemodel nc file into memory
%   icemodel.netcdf.nctype2mat               - Map NetCDF data types to MATLAB data types
%   README.md
%   icemodel.netcdf.redefatt
%   icemodel.netcdf.setfilename
%   icemodel.netcdf.writedims                - Write dimensions to icemodel nc file
%   icemodel.netcdf.writeice1                - Write ice1 data to an icemodel nc file
%   icemodel.netcdf.writeice2                - Write ice2 data to icemodel nc file
%
%   +NETCDF/+DEFAULTS
%   icemodel.netcdf.defaults.axes            - Define the grid and time dimension units. Use empty char '' for variables
%   icemodel.netcdf.defaults.cfStandardNames - Load the official CF Standard Name Table as a set
%   icemodel.netcdf.defaults.longnames       - Define the grid and time dimension names
%   icemodel.netcdf.defaults.standardnames   - Not all values have standard names, so I constructed some
%   icemodel.netcdf.defaults.units           - Define the grid and time dimension units
%   icemodel.netcdf.defaults.variable        - Canonical {standard_name, long_name, unit, is_cf} for a channel
%   icemodel.netcdf.defaults.variables       - Canonical variable-metadata map for every icemodel channel
%   icemodel.netcdf.defaults.varnames        - Define the grid and time dimension names
%
%   updatecontents.m generated this file on 14 Sep 2026 at 20:13:59.
