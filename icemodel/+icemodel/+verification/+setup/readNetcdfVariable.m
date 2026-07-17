function values = readNetcdfVariable(pathname, varname, kwargs)
   %READNETCDFVARIABLE Read one ESM-SnowMIP-style NetCDF variable with NaN-fill.
   %
   %  values = icemodel.verification.setup.readNetcdfVariable(pathname, varname)
   %  values = icemodel.verification.setup.readNetcdfVariable(pathname, varname, ...
   %     optional=true, ntime=N)
   %
   %  Reads a NetCDF variable and converts upstream sentinel fill values
   %  (<= -900) to NaN. ESM-SnowMIP uses -9999 / -999.99 as missing-value
   %  sentinels in its met and observation files.
   %
   %  When optional=true and the variable is absent from the file, the
   %  function returns an ntime-element NaN column instead of erroring.
   %  This is used by builders that handle site-by-site channel
   %  variability (e.g. some ESM-SnowMIP sites lack 'ts' or 'albs').
   %
   %  Inputs
   %    pathname : string
   %        Absolute path to a NetCDF file.
   %    varname : string
   %        Variable name to read.
   %
   %  Name-value
   %    optional : logical (default false)
   %        Return NaN column if the variable is missing.
   %    ntime : double (default 0)
   %        Length of the NaN fill column when optional=true and
   %        the variable is missing.
   %
   %  Returns
   %    values : double column with sentinels replaced by NaN.
   %
   % See also: icemodel.verification.setup.readNetcdfTime,
   %  icemodel.verification.setup.readBestSnowDepth

   arguments
      pathname (1, 1) string
      varname  (1, 1) string
      kwargs.optional (1, 1) logical = false
      kwargs.ntime    (1, 1) double = 0
   end

   if kwargs.optional
      info = ncinfo(pathname);
      if ~any(string({info.Variables.Name}) == varname)
         values = nan(kwargs.ntime, 1);
         return
      end
   end

   values = double(ncread(pathname, char(varname)));
   values(values <= -900) = NaN;
   values = values(:);
end
