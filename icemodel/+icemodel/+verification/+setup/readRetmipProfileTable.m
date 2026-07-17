function [data, metadata] = readRetmipProfileTable(filename)
   %READRETMIPPROFILETABLE Read a RetMIP initial profile table.
   %
   %  [data, metadata] = ...
   %     icemodel.verification.setup.readRetmipProfileTable(filename)
   %
   % Role
   %  Parser for RetMIP initial density/temperature/LWC profiles. It validates
   %  that the profile has a depth axis and at least one firn-state variable.

   arguments
      filename (1, 1) string
   end

   % RetMIP profile files use the same compact quoted semicolon form as the
   % surface protocol files, but they do not have a time axis.
   data = readSimpleDelimitedTable(filename);
   names = string(data.Properties.VariableNames);
   lower_names = lower(names);

   % RetMIP profile files must carry depth plus at least one initial state axis.
   has_depth = any(ismember(lower_names, ["depth", "z", "depth_m"]));
   has_state = any(contains(lower_names, ["density", "rho", ...
      "temperature", "temp", "lwc"]));
   if ~has_depth || ~has_state
      error('icemodel:verification:readRetmipProfileTable:badProfile', ...
         'RetMIP profile table needs depth and density/temperature/LWC columns');
   end

   metadata = struct( ...
      'filename', filename, ...
      'variables', names);
end

function data = readSimpleDelimitedTable(filename)
   %READSIMPLEDELIMITEDTABLE Parse RetMIP's small quoted .tab profile tables.
   lines = strip(readlines(filename));
   lines = lines(strlength(lines) > 0);
   if isempty(lines)
      error('icemodel:verification:readRetmipProfileTable:emptyFile', ...
         'RetMIP profile table is empty');
   end

   % Choose the delimiter from the header, then remove RetMIP's row-level quotes
   % before splitting numeric profile rows.
   delimiter = sprintf('\t');
   if contains(lines(1), ";")
      delimiter = ";";
   end
   lines = erase(lines, '"');
   names = split(lines(1), delimiter).';
   raw = strings(numel(lines) - 1, numel(names));
   for n = 2:numel(lines)
      raw(n - 1, :) = split(lines(n), delimiter).';
   end

   % Profile files are numeric after the header; keep nonnumeric columns as
   % strings if a future protocol extension adds labels.
   data = table();
   for k = 1:numel(names)
      values = raw(:, k);
      nums = str2double(values);
      if all(~isnan(nums))
         data.(matlab.lang.makeValidName(names(k))) = nums;
      else
         data.(matlab.lang.makeValidName(names(k))) = values;
      end
   end
   data.Properties.VariableNames = cellstr(names);
   data = canonicalProfileNames(data);
end

function data = canonicalProfileNames(data)
   %CANONICALPROFILENAMES Use verification variable names for stored columns.
   data = renameProfileAlias(data, ["depth", "z", "depth_m"], "depth");
   data = renameProfileAlias(data, ["density", "density_kgm3", "rho"], ...
      "density");
   data = renameProfileAlias(data, ["temperature", "temp", "t"], ...
      "subsurface_temperature");
   data = renameProfileAlias(data, ["lwc", "liquid_water_content"], "lwc");
end

function data = renameProfileAlias(data, aliases, canonical)
   %RENAMEPROFILEALIAS Rename the first accepted alias when needed.
   names = string(data.Properties.VariableNames);
   if ismember(canonical, names)
      return
   end
   lower_names = lower(names);
   hit = find(ismember(lower_names, lower(aliases)), 1);
   if ~isempty(hit)
      data = renamevars(data, names(hit), canonical);
   end
end
