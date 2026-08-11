function hist = initialize_coupler_history()
   %INITIALIZE_COUPLER_HISTORY Empty iterate history for the coupler accelerators.
   %
   % NaN means "no history yet". aitkenscalar and secantscalar both fall back
   % to the relaxed Picard step until enough iterations have run, so the
   % couplers do not need to special-case the first two passes.
   %
   % See also: icemodel.couplers.accelerate_coupler_iterate
   %
   %#codegen

   hist = struct('Ts_1', nan, 'Ts_2', nan, 'Ts_prev', nan, 'res_prev', nan);
end
