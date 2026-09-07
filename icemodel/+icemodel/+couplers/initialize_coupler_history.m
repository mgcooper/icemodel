function hist = initialize_coupler_history()
   %INITIALIZE_COUPLER_HISTORY Initialize the coupler-iteration history.
   %
   %  hist = icemodel.couplers.initialize_coupler_history()
   %
   % Ts_1 and Ts_2 start as NaN. Aitken acceleration uses the third Picard
   % iterate, after both values are available. Ts_prev and res_prev start as
   % NaN. The secant step uses the second residual when it brackets zero.
   %
   % See also: icemodel.couplers.accelerate_coupler_iterate
   %
   %#codegen

   hist = struct('Ts_1', nan, 'Ts_2', nan, 'Ts_prev', nan, 'res_prev', nan);
end
