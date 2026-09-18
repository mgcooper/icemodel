function policy = perfMeasurementPolicy()
   %PERFMEASUREMENTPOLICY Formal timing gate thresholds.
   %
   %  policy = icemodel.test.helpers.perfMeasurementPolicy()
   %
   % Every consumer of a formal timing gate reads it from here, so the
   % measurement protocol cannot drift between the per-case validity check
   % and the run-level ambient anchor.
   %
   % Fields:
   %  max_dispersion - largest accepted max/median sample ratio. One slow
   %                   outlier above this gate marks interference.
   %  anchor_tol     - largest accepted |anchor ratio - 1|. The ambient
   %                   anchor re-measures the first executed case at the
   %                   end of the run; a larger drift records
   %                   meta.ambient_stable = false once for the run. Every
   %                   case keeps its own comparison verdict, and the
   %                   run fails the measurement quality conditions of
   %                   perfMeasurementQuality.
   %  tol_perf       - two-sided fractional noise budget for a case
   %                   verdict. run_perf_suite and build_perf_baseline
   %                   default their tol_perf argument from it, and the
   %                   A/A diagnostic (run_aa_acceptance) accepts B/A
   %                   ratios inside [1/(1 + tol_perf), 1 + tol_perf].
   %  load_average_max - largest one-minute load average, sampled at run
   %                   start before the first measurement, that a release
   %                   snapshot accepts in the rolling source's
   %                   attestation. 4.0 is the gate the pre-snow sweep
   %                   applied by hand before each timing run. Later
   %                   samples include the run's own subprocesses and are
   %                   recorded, not gated.
   %
   % See also: icemodel.test.helpers.perfSampleValidity,
   %  icemodel.test.helpers.perfMeasurementQuality,
   %  icemodel.test.helpers.assertReleasePerfBaselineSource, run_perf_suite

   policy = struct( ...
      'max_dispersion', 1.5, ...
      'anchor_tol', 0.15, ...
      'tol_perf', 0.20, ...
      'load_average_max', 4.0);
end
