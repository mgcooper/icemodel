function balance = evaluate_surface_energy_balance( ...
      chi, albedo, Qsi, Qli, Qle, Qh, Qe, Qc, Qa, Qm)
   %EVALUATE_SURFACE_ENERGY_BALANCE Evaluate the SEB from known flux terms.
   %
   %  balance = icemodel.surface.evaluate_surface_energy_balance(
   %     chi, albedo, Qsi, Qli, Qle, Qh, Qe, Qc, Qa, Qm)
   %
   % All flux terms must already be evaluated at the current surface state
   % before calling this function. The balance is:
   %
   %   balance = chi*Qsi*(1-albedo) + emiss*Qli + Qle + Qh + Qe + Qc + Qa - Qm
   %
   % Pass Qm = 0 to obtain the energy available for melt when no melt is
   % assumed (the pre-melt residual). Pass Qm as the diagnosed melt energy to
   % obtain the closed balance (should be ~0 at steady state).
   %
   % This is the downstream assembler in the diagnostic SEB chain:
   %   surface_energy_balance_terms -> diagnose_melt_freeze_energy -> here
   %
   % See also: icemodel.surface.diagnose_surface_energy_fluxes,
   %           icemodel.surface.surface_energy_balance_terms,
   %           icemodel.surface.diagnose_melt_freeze_energy
   %
   %#codegen
   persistent emiss
   if isempty(emiss)
      emiss = icemodel.parameterLookup('emiss');
   end
   balance = chi * Qsi * (1.0 - albedo) ...
      + emiss * Qli + Qle + Qh + Qe + Qc + Qa - Qm;
end
