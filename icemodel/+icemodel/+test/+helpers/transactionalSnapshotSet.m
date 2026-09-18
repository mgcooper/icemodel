function baseline = transactionalSnapshotSet( ...
      kind, baseline_tag, models, simyear, snapshotter, kwargs)
   %TRANSACTIONALSNAPSHOTSET Create and validate an aggregate release set.
   %
   %  baseline = icemodel.test.helpers.transactionalSnapshotSet( ...
   %     kind, baseline_tag, models, simyear, snapshotter)
   %
   % A release snapshot set covers every formal model. A partial set would let
   % a later rolling generation fill the gap, so any failure removes every
   % snapshot this call created and rethrows.
   %
   % Input
   %  kind           "regression" or "perf".
   %  baseline_tag   Release tag the snapshots belong to.
   %  models         Formal models to snapshot.
   %  simyear        Benchmark year for perf snapshots.
   %  snapshotter    Called as snapshotter(kind, tag, model, simyear).
   %
   % Name-value
   %  loader          Reloads one saved snapshot for validation.
   %  remover         Removes one snapshot this call created. This default is
   %                  the only definition; resolveBootstrapRelease forwards a
   %                  remover only when its caller supplies one.
   %
   % See also: icemodel.test.helpers.resolveBootstrapRelease,
   %  icemodel.test.helpers.removeReleaseSnapshotArtifacts

   arguments
      kind (1, 1) string {mustBeMember(kind, ["regression", "perf"])}
      baseline_tag (1, 1) string
      models string
      simyear (1, 1) double
      snapshotter (1, 1) function_handle
      kwargs.loader (1, 1) function_handle = @loadReleaseSnapshot
      kwargs.remover (1, 1) function_handle = @removeReleaseSnapshot
   end

   % Track what this call created, separately from what it was asked to
   % create, so the rollback below removes only its own snapshots and never a
   % preserved release file that already existed.
   created_models = strings(numel(models), 1);
   n_created = 0;
   try
      for model = reshape(models, 1, [])
         snapshotter(kind, baseline_tag, model, simyear);
         n_created = n_created + 1;
         created_models(n_created) = model;
      end

      % Reload every saved file. A snapshot that writes but reloads without
      % rows would otherwise become an accepted empty release baseline.
      tables = arrayfun(@(model) kwargs.loader( ...
         kind, baseline_tag, model, simyear), models, ...
         'UniformOutput', false);
      if any(cellfun(@isempty, tables))
         error('icemodel:test:emptyReleaseSnapshot', ...
            'A saved %s release snapshot reloaded without rows.', kind)
      end
      baseline = vertcat(tables{:});
   catch err
      % Remove every snapshot this call created, then report the original
      % failure. A removal that also fails is chained onto it, so one bad
      % cleanup cannot hide the cause.
      for model = created_models(1:n_created).'
         try
            kwargs.remover(kind, baseline_tag, model, simyear);
         catch cleanup_err
            err = addCause(err, cleanup_err);
         end
      end
      rethrow(err)
   end
end

function baseline = loadReleaseSnapshot(kind, baseline_tag, smbmodel, simyear)
   %LOADRELEASESNAPSHOT Load one saved snapshot for aggregate validation.
   if kind == "perf"
      baseline = icemodel.test.helpers.loadBaseline(kind, ...
         smbmodel=smbmodel, baseline_tag=baseline_tag, simyear=simyear);
   else
      baseline = icemodel.test.helpers.loadBaseline(kind, ...
         smbmodel=smbmodel, baseline_tag=baseline_tag);
   end
end

function removeReleaseSnapshot(kind, baseline_tag, smbmodel, simyear)
   %REMOVERELEASESNAPSHOT Remove one newly created snapshot and sidecar.
   pathname = icemodel.test.helpers.baselineFilePath(kind, ...
      smbmodel=smbmodel, baseline_tag=baseline_tag, simyear=simyear);
   icemodel.test.helpers.removeReleaseSnapshotArtifacts(pathname)
end
