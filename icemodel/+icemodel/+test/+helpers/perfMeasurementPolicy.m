function policy = perfMeasurementPolicy()
   %PERFMEASUREMENTPOLICY Single source for the formal timing gates.
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
   %                   end of the run; a larger drift marks every verdict
   %                   in the run ambient-invalid.
   %
   % See also: icemodel.test.helpers.perfSampleValidity, run_perf_suite

   policy = struct( ...
      'max_dispersion', 1.5, ...
      'anchor_tol', 0.15);
end
