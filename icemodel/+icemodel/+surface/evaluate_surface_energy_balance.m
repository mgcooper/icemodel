function balance = evaluate_surface_energy_balance(Qsn, Qln, Qh, Qe, Qc, Qa, Qm)
   %EVALUATE_SURFACE_ENERGY_BALANCE Evaluate the SEB from known flux terms.
   %
   %  balance = icemodel.surface.evaluate_surface_energy_balance(
   %     Qsn, Qln, Qh, Qe, Qc, Qa, Qm)
   %
   % You must evaluate all flux terms at the current surface state before you
   % call this function. The balance is:
   %
   %   balance = Qsn + Qln + Qh + Qe + Qc + Qa - Qm
   %
   % where the pre-computed net fluxes are:
   %   Qsn = chi * Qsi * (1 - albedo)           [W m^-2]  net shortwave
   %   Qln = emiss * Qli - emiss * SB * T_sfc^4 [W m^-2]  net longwave
   %
   % Pass Qm = 0 to get the energy available for melt with no melt (the
   % pre-melt residual). Pass Qm as the diagnosed melt energy to get the
   % closed balance, which is about 0 at steady state.
   %
   % This is the downstream assembler in the diagnostic SEB chain:
   %   surface_energy_balance_terms -> diagnose_melt_freeze_energy -> here
   %
   % See also: icemodel.surface.diagnose_surface_energy_balance,
   %           icemodel.surface.surface_energy_balance_terms,
   %           icemodel.surface.diagnose_melt_freeze_energy
   %
   %#codegen
   balance = Qsn + Qln + Qh + Qe + Qc + Qa - Qm;
end
