function [ambient_stable, anchor_ratio] = ambientAnchorVerdict( ...
      first_median, anchor_sample_times, anchor_valid, anchor_tol)
   %AMBIENTANCHORVERDICT Decide whether ambient conditions held for a run.
   %
   %  [ambient_stable, anchor_ratio] = ...
   %     icemodel.test.helpers.ambientAnchorVerdict(first_median, ...
   %     anchor_sample_times, anchor_valid, anchor_tol)
   %
   % The ambient anchor re-measures the first executed case at the end of
   % a perf run. The per-case dispersion gate cannot see load or
   % scheduling shifts that are steady WITHIN each case but different
   % ACROSS cases, so the anchor compares the same case against its own
   % first measurement. An invalid anchor sample set cannot certify
   % stability, so ANCHOR_VALID must be true for a stable verdict.
   %
   % Inputs
   %  first_median        - the case's median from its formal measurement [s]
   %  anchor_sample_times - the anchor re-measurement's samples [s]
   %  anchor_valid        - the anchor sample set's validity-gate verdict
   %  anchor_tol          - largest accepted |anchor_ratio - 1|
   %
   % Outputs
   %  ambient_stable - true when the anchor certifies the run's timings
   %  anchor_ratio   - anchor median / first median
   %
   % See also: icemodel.test.helpers.perfMeasurementPolicy,
   %  icemodel.test.helpers.perfSampleValidity, run_perf_suite

   anchor_median = median(anchor_sample_times, 'omitnan');
   anchor_ratio = anchor_median / first_median;
   ambient_stable = anchor_valid && isfinite(anchor_ratio) ...
      && abs(anchor_ratio - 1) <= anchor_tol;
end
