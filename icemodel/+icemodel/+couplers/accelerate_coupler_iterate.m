function [T_sfc_next, hist] = accelerate_coupler_iterate( ...
      hist, T_sfc_old, T_sfc_new, cpl_alpha, cpl_jumpmax, cpl_aitken)
   %ACCELERATE_COUPLER_ITERATE Accelerate one surface-temperature Picard step.
   %
   % Accelerates the Picard loop on T_sfc that every coupler runs.
   %
   % Two stages. Aitken (with relaxation as its fallback) handles the ordinary
   % case. Then, when the last two residuals bracket a root, the secant step
   % replaces a Picard iterate that is not contracting. secantscalar returns
   % the fallback when there is no bracket or no usable history, so the first
   % iterations behave exactly as Aitken alone.
   %
   % Inputs
   %  hist        - iterate history struct from initialize_coupler_history
   %  T_sfc_old   - surface temperature entering this iteration
   %  T_sfc_new   - Picard iterate produced by this iteration
   %  cpl_alpha   - relaxation weight on the Picard update
   %  cpl_jumpmax - largest accepted jump from T_sfc_old
   %  cpl_aitken  - enable both acceleration stages, Aitken and secant. False
   %                falls back to relaxation only.
   %
   % Outputs
   %  T_sfc_next  - accelerated surface temperature for the next iteration
   %  hist        - history advanced by one iteration
   %
   % See also: icemodel.numerics.aitkenscalar, icemodel.numerics.secantscalar
   %
   %#codegen

   cpl_res = T_sfc_new - T_sfc_old;

   % Relaxation is the innermost fallback: T_sfc_old + alpha * residual.
   T_sfc_fallback = icemodel.numerics.aitkenscalar(hist.T_sfc_2, ...
      hist.T_sfc_1, T_sfc_new, T_sfc_old + cpl_alpha * cpl_res, cpl_jumpmax, ...
      cpl_aitken);

   T_sfc_next = icemodel.numerics.secantscalar(hist.T_sfc_prev, ...
      hist.res_prev, T_sfc_old, cpl_res, T_sfc_fallback, cpl_jumpmax, ...
      cpl_aitken);

   % Advance both histories: Aitken needs the last two iterates, the secant
   % needs the last (iterate, residual) pair.
   hist.T_sfc_2 = hist.T_sfc_1;
   hist.T_sfc_1 = T_sfc_new;
   hist.T_sfc_prev = T_sfc_old;
   hist.res_prev = cpl_res;
end
