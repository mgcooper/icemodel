function revision = assertCommonBaselineRevision( ...
      kind, baseline_tag, models, simyear, kwargs)
   %ASSERTCOMMONBASELINEREVISION Require one source revision for a model set.
   %
   %  revision = icemodel.test.helpers.assertCommonBaselineRevision( ...
   %     kind, baseline_tag, models, simyear)
   %
   % See also: run_regression_suite, run_perf_suite,
   %  snapshot_regression_baseline, snapshot_perf_baseline

   arguments
      kind (1, 1) string {mustBeMember(kind, ["regression", "perf"])}
      baseline_tag (1, 1) string
      models string
      simyear (1, 1) double
      kwargs.loader (1, 1) function_handle = ...
         @icemodel.test.helpers.loadBaseline
   end

   revisions = strings(numel(models), 1);
   for k = 1:numel(models)
      [~, meta] = kwargs.loader(kind, smbmodel=models(k), ...
         baseline_tag=baseline_tag, simyear=simyear);
      if isstruct(meta) && isfield(meta, 'git_revision')
         revisions(k) = string(meta.git_revision);
      end
   end

   % Test each entry. A scalar blank test would accept an all-blank set,
   % because isblanktext returns one false for a string array and every entry
   % is then equal. Test ismissing as well, because strlength returns NaN for
   % a missing string and NaN == 0 is false.
   if any(ismissing(revisions) | strlength(revisions) == 0) ...
         || numel(unique(revisions)) ~= 1
      error('icemodel:test:mixedReleaseBaselineProvenance', ...
         ['The %s %s model set does not share one nonblank source ', ...
         'revision.'], baseline_tag, kind)
   end
   revision = revisions(1);
end
