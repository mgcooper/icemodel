function artifact = stampArtifactMetadata(artifact)
   %STAMPARTIFACTMETADATA Add variable metadata to staged artifact tables.
   %
   %  artifact = icemodel.verification.setup.stampArtifactMetadata(artifact)
   %
   % Role
   %  Verification targets can contain nested structs of tables/timetables
   %  rather than one met-style timetable. This helper recursively stamps every
   %  table-like payload with canonical units/descriptions where known, while
   %  leaving source-key/string columns with blank units instead of rejecting the
   %  whole observational artifact.
   %
   % See also: icemodel.forcing.helpers.stampMetadata

   if istable(artifact) || istimetable(artifact)
      prior = tableMetadata(artifact);
      artifact = icemodel.forcing.helpers.stampMetadata(artifact, ...
         strict=false);
      artifact = restoreExistingTableMetadata(artifact, prior);
      artifact = stampIntervalSmbMetadata(artifact);
      return
   end

   if ~isstruct(artifact)
      return
   end

   % Recurse through struct arrays so target bundles and per-profile maps are
   % handled by the same code path.
   for n = 1:numel(artifact)
      names = fieldnames(artifact(n));
      for k = 1:numel(names)
         field = names{k};
         artifact(n).(field) = ...
            icemodel.verification.setup.stampArtifactMetadata( ...
            artifact(n).(field));
         if strcmp(field, 'metadata')
            % Keep saved target metadata column-oriented for interactive MAT
            % inspection. Table payloads are stamped separately above.
            artifact(n).(field) = icemodel.forcing.helpers.columnizeMetadata( ...
               artifact(n).(field));
         end
      end
   end
end

function T = stampIntervalSmbMetadata(T)
   %STAMPINTERVALSMBMETADATA Mark interval SMB as an accumulated observation.
   names = string(T.Properties.VariableNames);
   if ~all(ismember(["start_date", "end_date", "smb"], names))
      return
   end

   idx = find(names == "smb", 1);
   T.Properties.VariableUnits{idx} = 'm w.e.';
   T.Properties.VariableDescriptions{idx} = ...
      'accumulated surface mass balance over the observation interval';
   if isprop(T.Properties.CustomProperties, "StandardNames")
      T.Properties.CustomProperties.StandardNames(idx) = "";
   end
end

function prior = tableMetadata(T)
   %TABLEMETADATA Capture existing per-variable metadata before canonical fill.
   nvars = width(T);
   prior = struct( ...
      'VariableNames', {string(T.Properties.VariableNames)}, ...
      'VariableUnits', metadataValues(T.Properties.VariableUnits, nvars), ...
      'VariableDescriptions', ...
      metadataValues(T.Properties.VariableDescriptions, nvars), ...
      'StandardNames', strings(1, nvars), ...
      'HasStandardNames', ...
      isprop(T.Properties.CustomProperties, "StandardNames"));
   if isprop(T.Properties.CustomProperties, "StandardNames")
      prior.StandardNames = metadataValues( ...
         T.Properties.CustomProperties.StandardNames, nvars);
   end
end

function values = metadataValues(values, nvars)
   %METADATAVALUES Return one metadata string per variable.
   values = string(values);
   if numel(values) < nvars
      values(end + 1:nvars) = "";
   else
      values = values(1:nvars);
   end
end

function T = restoreExistingTableMetadata(T, prior)
   %RESTOREEXISTINGTABLEMETADATA Keep source-specific metadata already present.
   names = string(T.Properties.VariableNames);
   for k = 1:numel(names)
      idx = find(prior.VariableNames == names(k), 1);
      if isempty(idx)
         continue
      end
      if prior.VariableUnits(idx) ~= ""
         T.Properties.VariableUnits{k} = char(prior.VariableUnits(idx));
      end
      if prior.VariableDescriptions(idx) ~= ""
         T.Properties.VariableDescriptions{k} = ...
            char(prior.VariableDescriptions(idx));
      end
      if prior.HasStandardNames ...
            && isprop(T.Properties.CustomProperties, "StandardNames")
         T.Properties.CustomProperties.StandardNames(k) = ...
            prior.StandardNames(idx);
      end
   end
end
