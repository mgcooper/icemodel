function [dt_sum, d_liq, d_evp, d_lyr, d_rof, d_vap_liq, ...
      d_vap_ice, diag] = newtimestep(f_liq)
   %NEWTIMESTEP Initialize forcing-step accumulators and diagnostics.
   %
   %  [dt_sum, d_liq, d_evp, d_lyr, d_rof, d_vap_liq, d_vap_ice, diag] = ...
   %     icemodel.timestepping.newtimestep(f_liq)
   %
   % F_LIQ supplies the size and numeric class for the per-cell outputs.
   % DT_SUM is the elapsed substep time [s]. D_LIQ records melt and freeze.
   % D_EVP records the liquid part of surface vapor exchange. D_LYR records
   % remeshing. D_ROF records condensation overflow. D_VAP_LIQ and D_VAP_ICE
   % record all vapor-driven phase changes. DIAG records solver results.
   % All accumulators start at zero for the new forcing step.
   %
   % See also: icemodel.couplers.initialize_solver_diag,
   %  icemodel.timestepping.acceptsubstep
   %
   %#codegen

   % Reset the mass budget delta terms.
   d_liq = 0.0 * f_liq;
   d_evp = 0.0 * f_liq;
   d_lyr = 0.0 * f_liq;
   d_rof = 0.0;
   d_vap_liq = 0.0 * f_liq;
   d_vap_ice = 0.0 * f_liq;

   % Reset the substep time budget.
   dt_sum = 0.0;

   % Reset the solver diagnostics.
   diag = icemodel.couplers.initialize_solver_diag();
end
