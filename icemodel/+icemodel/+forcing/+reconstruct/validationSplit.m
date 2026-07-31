function split = validationSplit(years, kwargs)
   %VALIDATIONSPLIT Partition station years into selection and evaluation sets.
   %
   %  split = icemodel.forcing.reconstruct.validationSplit(2009:2023, ...
   %     station="kanm", seed=42)
   %  split = icemodel.forcing.reconstruct.validationSplit(2009:2023, ...
   %     station="kanm", seed=42, manifest_file="splits/kanm.json")
   %
   % Role
   %  Deterministic whole-year SELECTION / EVALUATION partition for the
   %  held-out validation protocol (gap-fill policy):
   %  fitting and method admission use selection years only; reported final
   %  numbers come from evaluation years only; the same observations never
   %  both select a method and grade it. Whole calendar years are the split
   %  unit so seasonal autocorrelation cannot leak across the boundary.
   %
   %  With manifest_file set, a persisted manifest WINS over recomputation:
   %  an existing, schema-valid file is replayed (POLICY B7 — the split is
   %  persisted, schema-verified, and replayed deterministically) only
   %  while it remains a disjoint,
   %  complete partition of the current record years. A missing file is
   %  created from this call's split. Callers change a persisted split only by
   %  deleting the manifest deliberately.
   %
   % Name-value
   %  station : string. Station identity recorded in the manifest; a loaded
   %     manifest must match it (guards against crossed files).
   %  seed : nonnegative integer. Required; drives the deterministic
   %     permutation so the split is reproducible without any ambient state.
   %  selection_fraction : scalar in (0, 1). Fraction of years assigned to
   %     the selection set (default 0.7); at least one year always lands in
   %     each set when two or more years exist.
   %  manifest_file : string. Optional JSON persistence path.
   %
   % Returns
   %  split : struct with fields station, seed, selection_fraction,
   %     years_selection, years_evaluation (ascending double rows).
   %
   % See also: icemodel.forcing.reconstruct.syntheticMissingness

   arguments
      years (1, :) double {mustBeInteger, mustBeNonempty}
      kwargs.station (1, 1) string
      kwargs.seed (1, 1) double {mustBeInteger, mustBeNonnegative}
      kwargs.selection_fraction (1, 1) double = ...
         icemodel.forcing.reconstruct.setopts().selection_fraction
      kwargs.manifest_file (1, 1) string = ""
   end
   mustBeInRange(kwargs.selection_fraction, 0, 1, 'exclusive')

   % A persisted manifest is the source of truth for replays.
   if kwargs.manifest_file ~= "" && isfile(kwargs.manifest_file)
      try
         split = jsondecode(fileread(kwargs.manifest_file));
      catch
         error('icemodel:reconstruct:validationSplit:invalidManifest', ...
            'split manifest is unreadable: %s', kwargs.manifest_file);
      end
      required = ["station", "seed", "selection_fraction", ...
         "years_selection", "years_evaluation"];
      if ~isstruct(split) || ~isscalar(split) ...
            || ~all(isfield(split, required))
         error('icemodel:reconstruct:validationSplit:invalidManifest', ...
            'split manifest has an invalid schema: %s', kwargs.manifest_file);
      end
      split.station = string(split.station);
      if ~isscalar(split.station)
         error('icemodel:reconstruct:validationSplit:invalidManifest', ...
            'split manifest station must be scalar: %s', ...
            kwargs.manifest_file);
      end
      if split.station ~= kwargs.station
         error('icemodel:reconstruct:validationSplit:stationMismatch', ...
            'manifest %s belongs to station %s, not %s', ...
            kwargs.manifest_file, split.station, kwargs.station);
      end
      split = validatePersistedSplit(split, years, kwargs.manifest_file);
      return
   end

   % Deterministic seeded permutation of the unique year list; at least one
   % year lands in each set so neither protocol side is ever empty.
   pool = unique(years);
   stream = RandStream('mt19937ar', 'Seed', kwargs.seed);
   order = pool(randperm(stream, numel(pool)));
   n_selection = round(kwargs.selection_fraction * numel(pool));
   n_selection = min(max(n_selection, 1), max(numel(pool) - 1, 1));
   split = struct( ...
      'station', kwargs.station, ...
      'seed', kwargs.seed, ...
      'selection_fraction', kwargs.selection_fraction, ...
      'years_selection', sort(order(1:n_selection)), ...
      'years_evaluation', sort(order(n_selection + 1:end)));

   if kwargs.manifest_file ~= ""
      % Persist through a plain JSON write; the read path above is the
      % replay contract.
      icemodel.helpers.ensureDirExists(fileparts(kwargs.manifest_file));
      fid = fopen(kwargs.manifest_file, 'w');
      cleaner = onCleanup(@() fclose(fid));
      fprintf(fid, '%s', jsonencode(split));
      clear cleaner
   end
end

function split = validatePersistedSplit(split, years, filename)
   %VALIDATEPERSISTEDSPLIT Require one complete current-year partition.
   selection = split.years_selection;
   evaluation = split.years_evaluation;
   valid_years = @(value) isnumeric(value) ...
      && (isempty(value) || isvector(value)) ...
      && all(isfinite(value)) && all(fix(value) == value) ...
      && numel(unique(value)) == numel(value);
   valid_scalars = isnumeric(split.seed) && isscalar(split.seed) ...
      && isfinite(split.seed) && split.seed >= 0 ...
      && fix(split.seed) == split.seed ...
      && isnumeric(split.selection_fraction) ...
      && isscalar(split.selection_fraction) ...
      && isfinite(split.selection_fraction) ...
      && split.selection_fraction > 0 && split.selection_fraction < 1;
   if ~isscalar(split.station) || ~valid_scalars ...
         || ~valid_years(selection) || ~valid_years(evaluation)
      error('icemodel:reconstruct:validationSplit:invalidManifest', ...
         'split manifest has invalid field values: %s', filename);
   end

   % A replay may reorder fields, but it may not overlap, omit, or introduce
   % years relative to the record currently being reconstructed.
   selection = sort(reshape(selection, 1, []));
   evaluation = sort(reshape(evaluation, 1, []));
   pool = unique(years);
   complete = isequal(sort([selection, evaluation]), pool);
   both_sides = numel(pool) < 2 ...
      || (~isempty(selection) && ~isempty(evaluation));
   if ~isempty(intersect(selection, evaluation)) || ~complete || ~both_sides
      error('icemodel:reconstruct:validationSplit:invalidManifest', ...
         ['split manifest must be a disjoint, complete partition of the ' ...
         'current years: %s'], filename);
   end
   split.years_selection = selection;
   split.years_evaluation = evaluation;
end
