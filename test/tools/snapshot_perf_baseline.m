function PerfBaseline = snapshot_perf_baseline(kwargs)
   %SNAPSHOT_PERF_BASELINE Save release snapshot from the rolling perf baseline.
   %
   %  PerfBaseline = snapshot_perf_baseline(baseline_tag="v1.1")
   %  PerfBaseline = snapshot_perf_baseline(baseline_tag="v1.1", simyear=2016)
   %
   % Use this when the current rolling perf baseline should be frozen as a
   % named release baseline. It does not rerun the model. It copies the current
   % rolling baseline into a versioned release file. A custom
   % OUTPUT_FILE is supported only when SMBMODEL resolves to one concrete
   % formal model. Existing release files remain immutable even when the
   % OVERWRITE option is true. snapshotBaseline requires a clean worktree,
   % apart from the release files this sequence writes, for a managed
   % release file, refuses a rolling source that carries the -dirty suffix,
   % and applies the release-only measurement conditions.

   arguments (Input)

      kwargs.baseline_tag (1, :) string

      kwargs.smbmodel (1, :) string ...
         {icemodel.validators.mustBeTestSmbmodelSelector(kwargs.smbmodel)} ...
         = "all"

      kwargs.simyear (1, 1) double {mustBeInteger, mustBePositive} ...
         = 2016

      kwargs.overwrite (1, 1) logical ...
         = false

      kwargs.output_file string = string.empty()
   end

   % Deal out arguments.
   [baseline_tag, smbmodel, simyear, overwrite, output_file] ...
      = deal(kwargs.baseline_tag, kwargs.smbmodel, kwargs.simyear, ...
      kwargs.overwrite, kwargs.output_file);

   % Expand the requested formal model selector once at the entrypoint.
   models = icemodel.test.helpers.resolveRequestedSmbmodels(smbmodel);

   % A custom output file is only coherent for one concrete model snapshot.
   if numel(models) > 1 && ~isblanktext(output_file)
      error(['output_file overrides only one release file. Omit it when ', ...
         'smbmodel expands to more than one formal model.'])
   end

   % Use the caller's custom output only for single-model snapshots.
   if numel(models) > 1
      output_file = "";
   end

   snapshotter = @(kind, tag, model, year) ...
      icemodel.test.helpers.snapshotBaseline( ...
      kind, tag, model, overwrite, output_file, year);
   if isblanktext(output_file)
      PerfBaseline = icemodel.test.helpers.transactionalSnapshotSet( ...
         "perf", baseline_tag, models, simyear, snapshotter);
   else
      % A custom output file needs its own loader and remover, because the
      % defaults resolve the managed baseline path, not this file.
      PerfBaseline = icemodel.test.helpers.transactionalSnapshotSet( ...
         "perf", baseline_tag, models, simyear, snapshotter, ...
         loader=@(kind, tag, model, year) ...
         icemodel.test.helpers.loadBaseline(kind, ...
         filename=output_file, smbmodel=model, baseline_tag=tag, ...
         simyear=year), ...
         remover=@(~, ~, ~, ~) ...
         icemodel.test.helpers.removeReleaseSnapshotArtifacts(output_file));
   end
end
