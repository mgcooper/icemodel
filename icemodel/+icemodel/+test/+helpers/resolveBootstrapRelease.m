function baseline = resolveBootstrapRelease( ...
      kind, baseline_tag, smbmodel, simyear, kwargs)
   %RESOLVEBOOTSTRAPRELEASE Load or create registered release baselines.
   %
   %  baseline = icemodel.test.helpers.resolveBootstrapRelease( ...
   %     kind, baseline_tag, smbmodel, simyear)
   %
   % For an aggregate model request, every immutable model snapshot must be
   % present or absent. A partial release set is rejected.
   %
   % Inputs
   %  kind         - "regression" or "perf"
   %  baseline_tag - registered release tag
   %  smbmodel     - formal model name or "all"
   %  simyear      - performance baseline year
   %
   % Output
   %  baseline - combined baseline table for the requested models
   %
   % Name-value
   %  policy_resolver, loader, regression_snapshotter, perf_snapshotter,
   %  snapshot_remover  Test seams. A test replaces one to drive a failure
   %  path. The first four default to the production call this function makes.
   %  snapshot_remover defaults to blank, which keeps transactionalSnapshotSet's
   %  own removal default.
   %
   % See also: run_test_bootstrap, snapshot_regression_baseline,
   %  snapshot_perf_baseline

   arguments
      kind (1, 1) string {mustBeMember(kind, ["regression", "perf"])}
      baseline_tag (1, 1) string
      smbmodel (1, :) string ...
         {icemodel.validators.mustBeTestSmbmodelSelector(smbmodel)}
      simyear (1, 1) double {mustBeInteger, mustBePositive}
      kwargs.policy_resolver (1, 1) function_handle = ...
         @icemodel.test.helpers.formalBaselinePolicy
      kwargs.loader (1, 1) function_handle = ...
         @icemodel.test.helpers.loadBaseline
      kwargs.regression_snapshotter (1, 1) function_handle = ...
         @snapshot_regression_baseline
      kwargs.perf_snapshotter (1, 1) function_handle = ...
         @snapshot_perf_baseline
      % Blank keeps transactionalSnapshotSet's own removal default, so the
      % removal exists in one place. Tests pass a stub here.
      kwargs.snapshot_remover function_handle = function_handle.empty()
   end

   policy = kwargs.policy_resolver(baseline_tag);
   models = icemodel.test.helpers.resolveRequestedSmbmodels(smbmodel);
   tables = loadModelSet(kwargs.loader, kind, baseline_tag, models, simyear);
   missing = cellfun(@isempty, tables);

   % Reject a partial set because its files can come from different rolling
   % generations.
   if any(missing) && ~all(missing)
      error('icemodel:test:partialReleaseBaseline', ...
         'Release %s has only part of the requested %s baseline set.', ...
         baseline_tag, kind)
   end

   % Create all missing snapshots as one rollback-protected set.
   if any(missing) && policy.snapshot_from_rolling
      if kind == "regression"
         snapshotter = @(~, tag, model, ~) ...
            kwargs.regression_snapshotter( ...
            baseline_tag=tag, smbmodel=model);
      else
         snapshotter = @(~, tag, model, year) ...
            kwargs.perf_snapshotter( ...
            baseline_tag=tag, smbmodel=model, simyear=year);
      end
      loader = @(snapshot_kind, tag, model, year) kwargs.loader( ...
         snapshot_kind, smbmodel=model, baseline_tag=tag, simyear=year);
      % Forward a remover only when the caller supplied one.
      remover_args = {};
      if ~isempty(kwargs.snapshot_remover)
         remover_args = {'remover', kwargs.snapshot_remover};
      end
      baseline = icemodel.test.helpers.transactionalSnapshotSet( ...
         kind, baseline_tag, models(missing), simyear, snapshotter, ...
         'loader', loader, remover_args{:});
      return
   elseif any(missing)
      error('icemodel:test:preservedReleaseMissing', ...
         'Registered immutable %s release %s is missing.', kind, baseline_tag)
   end

   baseline = vertcat(tables{:});
end

function tables = loadModelSet(loader, kind, baseline_tag, models, simyear)
   %LOADMODELSET Load each scalar model table.
   tables = cell(numel(models), 1);
   for k = 1:numel(models)
      tables{k} = loader(kind, smbmodel=models(k), ...
         baseline_tag=baseline_tag, simyear=simyear);
   end
end
