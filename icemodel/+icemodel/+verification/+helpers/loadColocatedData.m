function bundle = loadColocatedData(manifest, source, kwargs)
   %LOADCOLOCATEDDATA Assemble a timeseries bundle from staged per-source files.
   %
   %  bundle = icemodel.verification.helpers.loadColocatedData(manifest, "promice")
   %  bundle = icemodel.verification.helpers.loadColocatedData(manifest, "racmo")
   %
   % Inputs
   %  manifest         Resolved firn case manifest. Its `colocation` field
   %                   records which source legs were staged and the per-year
   %                   Data filenames under data/input/userdata.
   %  source           Source id to load. "racmo" loads the co-located RCM
   %                   reference Data (the standard use). "promice" loads the
   %                   eval target Data, used only for legacy PROMICE fixtures
   %                   staged before the bundled observations.mat eval contract;
   %                   freshly staged cases bundle observations.mat instead.
   %  input_data_root  Optional base input-data root. When blank, the standard
   %                   chain (inputDataRoot) resolves it.
   %
   % Outputs
   %  bundle   Struct with the verification target/reference contract:
   %             .format   "timeseries"
   %             .data     vertically concatenated per-year Data timetable
   %             .metadata provenance struct
   %
   % Role
   %  Operational helper for the firn lane. The forcing/reference side is never
   %  bundled: colocation is recorded as metadata (available sources + per-leg
   %  windows) pointing at individual per-source userdata files. This helper
   %  reconstitutes a comparison bundle on demand from those files, so
   %  comparecase/plotcase keep their timeseries contract for the RACMO
   %  reference (and the PROMICE-obs target on legacy fixtures) without a
   %  committed data bundle.
   %
   % See also: icemodel.verification.comparecase,
   %  icemodel.verification.helpers.resolveCandidateBundle,
   %  icemodel.verification.setup.importPromiceSites

   arguments
      manifest (1, 1) struct
      source (1, 1) string
      kwargs.input_data_root (1, 1) string = ""
   end

   ud_dir = fullfile(icemodel.verification.helpers.inputDataRoot( ...
      "input_data_root", kwargs.input_data_root), 'userdata');

   if ~isfield(manifest, 'colocation') || ~isfield(manifest.colocation, source)
      bundle = emptyBundle(sprintf('source "%s" not in manifest colocation', ...
         source));
      return
   end
   leg = manifest.colocation.(source);
   if ~isfield(leg, 'data_files') || isempty(leg.data_files)
      reason = 'leg not staged';
      if isfield(leg, 'reason') && strlength(string(leg.reason)) > 0
         reason = char(string(leg.reason));
      end
      bundle = emptyBundle(reason);
      return
   end

   files = string(leg.data_files);
   files = files(strlength(files) > 0);
   data = timetable.empty;
   for f = reshape(files, 1, [])
      pathname = fullfile(ud_dir, f);
      if ~isfile(pathname)
         continue
      end
      loaded = load(pathname, 'Data');
      if ~isfield(loaded, 'Data') || ~istimetable(loaded.Data)
         continue
      end
      if isempty(data)
         data = loaded.Data;
      else
         data = outerjoinrows(data, loaded.Data);
      end
   end

   if isempty(data)
      bundle = emptyBundle('no per-year Data files resolved on disk');
      return
   end

   bundle = struct('format', 'timeseries', 'data', data, ...
      'metadata', struct('source', char(source), ...
      'reconstituted_from', 'per-year userdata files'));
end

%% Local helpers
function bundle = emptyBundle(note)
   %EMPTYBUNDLE Explicit empty timeseries bundle with a reason note.
   bundle = struct('format', 'timeseries', 'data', timetable.empty, ...
      'metadata', struct('note', char(note)));
end

function combined = outerjoinrows(a, b)
   %OUTERJOINROWS Concatenate two Data timetables on a shared variable set.
   %
   % Per-year files share the same variable schema, so a row concatenation on
   % the union of times reconstructs the multi-year record. Reorder b's columns
   % to a's order and drop columns absent from a so vertcat succeeds.
   common = intersect(a.Properties.VariableNames, ...
      b.Properties.VariableNames, 'stable');
   combined = [a(:, common); b(:, common)];
   combined = sortrows(combined);
end
