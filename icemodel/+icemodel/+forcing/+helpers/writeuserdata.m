function filenames = writeuserdata(Data, site, source, kwargs)
   %WRITEUSERDATA Save a Data timetable as met-swap userdata files.
   %
   %  filenames = icemodel.forcing.helpers.writeuserdata(Data, site, source)
   %  filenames = ... writeuserdata(_, outdir=..., naming="window", ...
   %     dt_out="1h", overwrite=true)
   %
   % Saves DATA under the icemodel met-swap ("userdata") naming convention.
   % naming="yearly" (default): one file per calendar year,
   %
   %    <site>_<source>_<YYYY>[_<cadence>].mat
   %
   % naming="window": one file spanning the full time axis,
   %
   %    <site>_<source>_<YYYYMMDD>_<YYYYMMDD>[_<cadence>].mat
   %
   % icemodel.loadmet resolves either form (a window file bracketing the run
   % year is preferred; otherwise the per-year file).
   %
   % consumed by icemodel.loadmet when opts.userdata / opts.uservars
   % request that variables of the met file be swapped with the
   % corresponding columns of the userdata file. Each .mat file holds
   % one variable named Data, a timetable carrying location metadata as
   % table CustomProperties (X, Y, Lat, Lon, Elev, Slope, ScalarUnits).
   % Existing targets are additive no-ops unless overwrite=true. Userdata
   % defaults to hourly output at this shared writer boundary; pass dt_out=""
   % only when an explicitly documented source must retain native cadence.
   % Explicit replacement and wider-window cleanup emit warnings.
   %
   % OUTDIR defaults to icemodel.getpath('userdata') (demo/data/input/
   % userdata when the demo or test config is active) and is created
   % when it does not exist.
   %
   % Inputs
   %  Data   - timetable of evaluation/forcing variables with location
   %           CustomProperties attached
   %  site   - site name encoded in the filename (e.g. "kanm")
   %  source - data source encoded in the filename (e.g. "merra")
   %
   % Outputs
   %  filenames - string column of selected full paths (including reused files)
   %
   % See also: icemodel.forcing.helpers.writemet, icemodel.loadmet,
   %  icemodel.setopts

   arguments
      Data timetable
      site (1, 1) string
      source (1, 1) string
      kwargs.outdir (1, 1) string = ""
      kwargs.requiremetadata (1, 1) logical = true
      kwargs.naming (1, 1) string ...
         {mustBeMember(kwargs.naming, ["window", "yearly"])} = "yearly"
      kwargs.dt_out (1, 1) string ...
         {mustBeMember(kwargs.dt_out, ["", "1h"])} = "1h"
      kwargs.overwrite (1, 1) logical = false
   end

   % Keep builders source-native; select the public userdata cadence only at
   % the final shared artifact boundary.
   Data = resampleUserdataTimestep(Data, kwargs.dt_out);
   requested_cadence_s = outputCadenceSeconds(Data);
   cadence_suffix = userdataCadenceSuffix(requested_cadence_s);
   cadence_matches = @(filename) ...
      icemodel.forcing.helpers.artifactCadenceMatches( ...
      filename, "Data", requested_cadence_s);

   % Userdata files can contain source-specific comparison columns, so stamp
   % known channels while leaving truly source-key columns blank.
   Data = icemodel.forcing.helpers.stampMetadata(Data, strict=false);
   identity_matches = @(filename) ...
      icemodel.forcing.helpers.artifactIdentityMatches( ...
      filename, Data, "Data");
   candidate_matches = @(filename) cadence_matches(filename) ...
      && identity_matches(filename);

   if kwargs.requiremetadata
      requireLocationMetadata(Data)
   end

   outdir = kwargs.outdir;
   if outdir == ""
      outdir = string(icemodel.getpath('userdata'));
   end
   % Stage into the per-source subfolder userdata/<source>/ so the flat userdata/
   % folder does not sprawl; the runtime resolves this subfolder first
   % (icemodel.loadmet.resolveUserdataFile).
   outdir = fullfile(outdir, char(source));
   if ~isfolder(outdir)
      mkdir(outdir)
   end

   switch kwargs.naming
      case "window"
         % One file spanning the full time axis. The encoded YYYYMMDD-YYYYMMDD
         % period lets icemodel.loadmet pick the file bracketing a run year.
         t1 = char(min(Data.Time), 'yyyyMMdd');
         t2 = char(max(Data.Time), 'yyyyMMdd');
         name = sprintf('%s_%s_%s_%s%s.mat', ...
            site, source, t1, t2, cadence_suffix);
         filenames = fullfile(outdir, name);

         % Validate the exact target before broad reuse. Legacy runtime lookup
         % still prefers the widest covering Data file, but a stale exact file
         % must not remain silently selectable when its identity or cadence
         % conflicts with this request.
         if isfile(filenames) && ~kwargs.overwrite
            assertReusableUserdata(filenames, cadence_matches, ...
               identity_matches, "window");
         end

         % A compatible broader source window still satisfies a narrower
         % ordinary request without resampling or duplicating native Data.
         if ~kwargs.overwrite
            enclosing = ...
                icemodel.forcing.helpers.findEnclosingWindowFile( ...
                outdir, site + "_" + source, cadence_suffix + ".mat", ...
                min(Data.Time), max(Data.Time), ...
                accept_candidate=candidate_matches);
            % Return the same path loadmet would select for met swapping.
            if enclosing ~= ""
               filenames = fullfile(outdir, enclosing);
               return
            end
         end
         wrote = savedata(filenames, Data, kwargs.overwrite);
         if wrote
            % A successful wider refresh supersedes only contained shorter
            % windows for this exact site/source naming class.
             icemodel.forcing.helpers.pruneSupersededWindowFiles( ...
                filenames, site + "_" + source, cadence_suffix + ".mat", ...
                accept_candidate=candidate_matches);
         end
      case "yearly"
         years_present = unique(year(Data.Time));
         filenames = strings(numel(years_present), 1);
         for n = 1:numel(years_present)
             yyyy = years_present(n);
             filenames(n) = fullfile(outdir, ...
                sprintf('%s_%s_%d%s.mat', ...
                site, source, yyyy, cadence_suffix));
             if isfile(filenames(n)) && ~kwargs.overwrite
                assertReusableUserdata(filenames(n), cadence_matches, ...
                   identity_matches, "year");
             end
             savedata(filenames(n), Data(year(Data.Time) == yyyy, :), ...
                kwargs.overwrite)
         end
   end
end

function assertReusableUserdata(filename, cadence_matches, ...
      identity_matches, span)
   %ASSERTREUSABLEUSERDATA Reject an exact-name collision with other provenance.
   if ~cadence_matches(filename)
      error('icemodel:forcing:writeuserdata:cadenceConflict', ...
         ['Existing userdata artifact %s has a different or unknown cadence. ' ...
         'Pass overwrite=true to replace that exact %s.'], filename, span)
   end
   if ~identity_matches(filename)
      error('icemodel:forcing:writeuserdata:identityConflict', ...
         ['Existing userdata artifact %s has conflicting source, product, ' ...
         'schema, sampling-method, or point metadata. Pass overwrite=true ' ...
         'to replace that exact %s.'], filename, span)
   end
end

function suffix = userdataCadenceSuffix(cadence_s)
   %USERDATACADENCESUFFIX Keep hourly legacy names and identify native variants.
   if ~isfinite(cadence_s) || cadence_s <= 0
      error('icemodel:forcing:writeuserdata:unknownCadence', ...
         'userdata output cadence must be finite before naming the artifact')
   end
   if abs(cadence_s - 3600) < 1e-6
      suffix = "";
      return
   end
   if abs(cadence_s - 900) < 1e-6
      suffix = "_15m";
   elseif abs(cadence_s - 1800) < 1e-6
      suffix = "_30m";
   else
      suffix = "_" + string(round(cadence_s)) + "s";
   end
end

function cadence_s = outputCadenceSeconds(Data)
   %OUTPUTCADENCESECONDS Read the writer-stamped requested artifact cadence.
   cadence_s = NaN;
   metadata = Data.Properties.UserData;
   if isstruct(metadata) ...
         && isfield(metadata, 'userdata_resample_output_cadence_seconds')
      cadence_s = double(metadata.userdata_resample_output_cadence_seconds);
   end
end

%% Local functions
function Data = resampleUserdataTimestep(Data, dt_out)
   %RESAMPLEUSERDATATIMESTEP Select hourly or explicit native userdata cadence.
   source = Data;
   cadence_s = NaN;
   if height(source) > 1
      cadence_s = median(seconds(diff(source.Time)), 'omitnan');
   end

   if dt_out == ""
      Data = recordResampleProvenance( ...
         Data, source, "explicit_native", cadence_s);
      return
   end

   % Uniform top-of-hour data already satisfy the public artifact contract.
   aligned = height(source) < 2 ...
      || (all(diff(source.Time) == hours(1)) ...
      && all(source.Time == dateshift(source.Time, 'start', 'hour')));
   if aligned
      Data = recordResampleProvenance( ...
         Data, source, "native_1h_unchanged", cadence_s);
      return
   end

   % Finer native averages are aggregated into clock-hour bins. Coarser native
   % products are interpolated to hourly support; no current writer caller uses
   % that path, but keeping it defined makes the public default unconditional.
   if cadence_s <= 3600
      Data = aggregateHourly(source);
      policy = "hourly_mean";
   else
      Data = interpolateHourly(source);
      policy = "hourly_linear";
   end
   Data = recordResampleProvenance(Data, source, policy, cadence_s);
end

function Data = aggregateHourly(source)
   %AGGREGATEHOURLY Average finer source samples in clock-hour bins.
   names = string(source.Properties.VariableNames);
   has_direction = ismember("wdir", names);
   if has_direction
      direction = timetable(source.Time, sind(source.wdir), cosd(source.wdir), ...
         'VariableNames', {'sin_wdir', 'cos_wdir'});
      direction = retime(direction, 'hourly', 'mean');
      source = removevars(source, 'wdir');
   end
   Data = retime(source, 'hourly', 'mean');
   if has_direction
      Data.wdir = mod(atan2d(direction.sin_wdir, direction.cos_wdir), 360);
      Data = movevars(Data, 'wdir', 'Before', find(names == "wdir", 1));
   end
end

function Data = interpolateHourly(source)
   %INTERPOLATEHOURLY Interpolate coarser source support onto clock-hour rows.
   new_time = (dateshift(source.Time(1), 'start', 'hour'):hours(1): ...
      dateshift(source.Time(end), 'start', 'hour'))';
   names = string(source.Properties.VariableNames);
   has_direction = ismember("wdir", names);
   if has_direction
      direction = timetable(source.Time, sind(source.wdir), cosd(source.wdir), ...
         'VariableNames', {'sin_wdir', 'cos_wdir'});
      direction = retime(direction, new_time, 'linear');
      source = removevars(source, 'wdir');
   end
   Data = retime(source, new_time, 'linear');
   if has_direction
      Data.wdir = mod(atan2d(direction.sin_wdir, direction.cos_wdir), 360);
      Data = movevars(Data, 'wdir', 'Before', find(names == "wdir", 1));
   end
end

function Data = recordResampleProvenance(Data, source, policy, cadence_s)
   %RECORDRESAMPLEPROVENANCE Stamp source/output cadence facts for artifact QA.
   metadata = source.Properties.UserData;
   if isempty(metadata)
      metadata = struct();
   elseif ~isstruct(metadata)
      error('icemodel:forcing:writeuserdata:badMetadata', ...
         'Data.Properties.UserData must be empty or a metadata struct')
   end
   metadata.userdata_resample_policy = policy;
   metadata.userdata_resample_source_row_count = height(source);
   metadata.userdata_resample_source_cadence_seconds = cadence_s;
   metadata.userdata_resample_output_cadence_seconds = 3600;
   if policy == "explicit_native"
      metadata.userdata_resample_output_cadence_seconds = cadence_s;
   end
   Data.Properties.UserData = metadata;
end

function requireLocationMetadata(Data)
   %REQUIRELOCATIONMETADATA Assert the Data-file CustomProperties contract.
   needed = ["X", "Y", "Lat", "Lon", "Elev", "Slope", "ScalarUnits"];
   have = string(fieldnames(Data.Properties.CustomProperties));
   missing = setdiff(needed, have);
   if ~isempty(missing)
      error('icemodel:forcing:writeuserdata:missingMetadata', ...
         'Data is missing CustomProperties metadata: %s', ...
         strjoin(missing, ', '));
   end
end

function wrote = savedata(filename, Data, overwrite)
   %SAVEDATA Save DATA to FILENAME as a variable named Data.
   exists = isfile(filename);
   wrote = false;
   % Ordinary repeated writes reuse the exact current artifact bytes.
   if exists && ~overwrite
      return
   end
   % Explicit replacement is intentionally visible to setup callers.
   if exists
      warning('icemodel:forcing:writeuserdata:overwrite', ...
         'Replacing existing userdata artifact %s.', filename);
   end
   % Persist one exact provenance record in both supported read locations. The
   % adapter may derive cadence/location facts absent from incoming UserData.
   S.artifact_metadata = icemodel.forcing.helpers.artifactMetadata(Data);
   Data.Properties.UserData = S.artifact_metadata;
   S.Data = Data;
   save(filename, '-struct', 'S')
   wrote = true;
end
