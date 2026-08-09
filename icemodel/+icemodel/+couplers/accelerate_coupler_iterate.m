function [Ts_next, hist] = accelerate_coupler_iterate( ...
      hist, Ts_old, Ts_new, cpl_alpha, cpl_jumpmax, cpl_aitken)
   %ACCELERATE_COUPLER_ITERATE Accelerate one surface-temperature Picard step.
   %
   % Every coupler runs the same Picard loop on T_sfc and accelerates it the
   % same way, so the acceleration lives here instead of in each solver.
   %
   % Two stages. Aitken (with relaxation as its fallback) handles the ordinary
   % case. Then, when the last two residuals bracket a root, the secant step
   % replaces a Picard iterate that is not contracting. secantscalar returns
   % the fallback when there is no bracket or no usable history, so the first
   % iterations behave exactly as Aitken alone.
   %
   % Inputs
   %  hist       - iterate history struct from initialize_coupler_history
   %  Ts_old     - surface temperature entering this iteration
   %  Ts_new     - Picard iterate produced by this iteration
   %  cpl_alpha  - relaxation weight on the Picard update
   %  cpl_jumpmax - largest accepted jump from Ts_old
   %  cpl_aitken - enable both acceleration stages, Aitken and secant. False
   %               falls back to relaxation only.
   %
   % Outputs
   %  Ts_next - accelerated surface temperature for the next iteration
   %  hist    - history advanced by one iteration
   %
   % See also: icemodel.numerics.aitkenscalar, icemodel.numerics.secantscalar
   %
   %#codegen

   cpl_res = Ts_new - Ts_old;

   % Relaxation is the innermost fallback: Ts_old + alpha * residual.
   Ts_fallback = icemodel.numerics.aitkenscalar(hist.Ts_2, hist.Ts_1, ...
      Ts_new, Ts_old + cpl_alpha * cpl_res, cpl_jumpmax, cpl_aitken);

   Ts_next = icemodel.numerics.secantscalar(hist.Ts_prev, hist.res_prev, ...
      Ts_old, cpl_res, Ts_fallback, cpl_jumpmax, cpl_aitken);

   % Advance both histories: Aitken needs the last two iterates, the secant
   % needs the last (iterate, residual) pair.
   hist.Ts_2 = hist.Ts_1;
   hist.Ts_1 = Ts_new;
   hist.Ts_prev = Ts_old;
   hist.res_prev = cpl_res;
end
