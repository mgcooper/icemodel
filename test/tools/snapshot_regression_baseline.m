function RegressionBaseline = snapshot_regression_baseline(kwargs)
   %SNAPSHOT_REGRESSION_BASELINE Freeze the rolling regression baseline.
   %
   %  RegressionBaseline = snapshot_regression_baseline(baseline_tag="v1.1")
   %
   % Use this to freeze the current rolling regression baseline as a named
   % release baseline. It does not rerun the model. It copies the current
   % rolling baseline into a versioned release file. You can pass a custom
   % OUTPUT_FILE only when SMBMODEL resolves to one concrete formal model.
   % Existing release files stay immutable even when the OVERWRITE option is
   % true.

   arguments (Input)

      kwargs.baseline_tag (1, :) string

      kwargs.smbmodel (1, :) string ...
         {icemodel.validators.mustBeTestSmbmodelSelector(kwargs.smbmodel)} ...
         = "all"

      kwargs.overwrite (1, 1) logical ...
         = false

      kwargs.output_file string = string.empty()
   end

   % Deal out arguments.
   [baseline_tag, smbmodel, overwrite, output_file] = deal( ...
      kwargs.baseline_tag, kwargs.smbmodel, kwargs.overwrite, ...
      kwargs.output_file);

   % Expand the requested formal model selector once at the entrypoint.
   models = icemodel.test.helpers.resolveRequestedSmbmodels(smbmodel);
   if numel(models) > 1 && ~isblanktext(output_file)
      error(['output_file overrides only one release file. Omit it when ', ...
         'smbmodel expands to more than one formal model.'])
   end

   % Use the caller's custom output only for single-model snapshots.
   if numel(models) > 1
      output_file = "";
   end

   % Every model snapshot must come from the same rolling source revision.
   icemodel.test.helpers.assertCommonBaselineRevision( ...
      "regression", "rolling", models, NaN);

   snapshotter = @(kind, tag, model, simyear) ...
      icemodel.test.helpers.snapshotBaseline( ...
      kind, tag, model, overwrite, output_file, simyear);
   if isblanktext(output_file)
      RegressionBaseline = ...
         icemodel.test.helpers.transactionalSnapshotSet( ...
         "regression", baseline_tag, models, NaN, snapshotter, ...
         require_common_revision=true);
   else
      % A custom output file skips the common-revision check, because the
      % revision loader resolves the managed baseline path, not this file.
      RegressionBaseline = ...
         icemodel.test.helpers.transactionalSnapshotSet( ...
         "regression", baseline_tag, models, NaN, snapshotter, ...
         loader=@(kind, tag, model, simyear) ...
         icemodel.test.helpers.loadBaseline(kind, ...
         filename=output_file, smbmodel=model, baseline_tag=tag, ...
         simyear=simyear), ...
         remover=@(~, ~, ~, ~) ...
         icemodel.test.helpers.removeReleaseSnapshotArtifacts(output_file));
   end
end
