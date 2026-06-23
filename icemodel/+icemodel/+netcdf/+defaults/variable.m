function info = variable(name, kwargs)
   %VARIABLE Canonical {standard_name, long_name, unit, is_cf} for a channel.
   %
   %  info = icemodel.netcdf.defaults.variable(name)
   %  info = ... variable(name, validatecf=true)
   %
   % Looks NAME up in the single canonical variable-metadata map
   % (icemodel.netcdf.defaults.variables) and returns its metadata struct
   % with fields standard_name, long_name, unit, is_cf. Handles the indexed
   % ice-temperature string channels (ticeN -> land_ice_temperature [K]) and
   % thermistor-depth string channels (dticeN -> depth [m]) by pattern, so
   % an arbitrary-length string resolves without enumerating every index.
   %
   % NAME may be a scalar string/char (returns one struct) or a string array
   % (returns a struct array of the same size).
   %
   % Name-value
   %  validatecf : (false) assert every non-empty standard_name returned is
   %               in the official CF Standard Name Table
   %               (icemodel.netcdf.defaults.cfStandardNames); errors listing
   %               any standard_name we use that is NOT recognized by CF.
   %
   % An unknown channel errors: the canonical map must label every channel.
   %
   % See also: icemodel.netcdf.defaults.variables,
   %  icemodel.netcdf.defaults.cfStandardNames

   arguments
      name (1, :) string
      kwargs.validatecf (1, 1) logical = false
   end

   map = icemodel.netcdf.defaults.variables();

   info = repmat(emptyInfo(), 1, numel(name));
   for k = 1:numel(name)
      info(k) = resolveOne(name(k), map);
   end

   if kwargs.validatecf
      cf = icemodel.netcdf.defaults.cfStandardNames();
      % Preallocate to the worst case (every name invalid), then trim, so the
      % collected list never grows inside the loop.
      bad = strings(1, numel(info));
      n_bad = 0;
      for k = 1:numel(info)
         sn = string(info(k).standard_name);
         if strlength(sn) > 0 && info(k).is_cf && ~any(cf == sn)
            n_bad = n_bad + 1;
            bad(n_bad) = sn;
         end
      end
      bad = bad(1:n_bad);
      if ~isempty(bad)
         error('icemodel:netcdf:variable:nonCfStandardName', ...
            'standard_name(s) not in the official CF table: %s', ...
            strjoin(unique(bad), ', '))
      end
   end

   if isscalar(info)
      info = info(1);
   end
end

%% Local functions
function s = resolveOne(name, map)
   %RESOLVEONE Metadata for a single channel name, with pattern fallbacks.

   if isKey(map, name)
      s = map(name);
      return
   end

   % Indexed ice-temperature string: tice1..ticeN -> land_ice_temperature [K]
   if ~isempty(regexp(name, '^tice\d+$', 'once'))
      s = emptyInfo();
      s.standard_name = 'land_ice_temperature';
      s.long_name = sprintf('ice / firn temperature at thermistor %s', ...
         regexprep(name, '^tice', ''));
      s.unit = 'K';
      s.is_cf = true;
      return
   end

   % Indexed thermistor depths: dtice1..dticeN -> depth [m]
   if ~isempty(regexp(name, '^dtice\d+$', 'once'))
      s = emptyInfo();
      s.standard_name = 'depth';
      s.long_name = sprintf('depth of thermistor %s below the surface', ...
         regexprep(name, '^dtice', ''));
      s.unit = 'm';
      s.is_cf = true;   % "depth" is a CF standard name (and coordinate)
      return
   end

   error('icemodel:netcdf:variable:unknownChannel', ...
      ['no canonical metadata for channel "%s"; add it to ' ...
      'icemodel.netcdf.defaults.variables'], name)
end

function s = emptyInfo()
   %EMPTYINFO Prototype metadata struct.
   s = struct('standard_name', '', 'long_name', '', 'unit', '', ...
      'is_cf', false);
end
