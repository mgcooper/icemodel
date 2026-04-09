function residual = surface_energy_balance_residual(T_sfc, tair, Qsi, Qli, ...
      albedo, wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, roL, liqflag, chi, Qc, ...
      ro_sfc, snow_depth, opts)
   %SURFACE_ENERGY_BALANCE_RESIDUAL Return the SEB residual at T_sfc.
   %
   %  residual = icemodel.surface.surface_energy_balance_residual(...)
   %
   % This is the canonical nonlinear surface residual used by the touched
   % SEB solvers. It includes absorbed shortwave, net longwave, conductive
   % heat, sensible heat, latent heat, and precipitation-advected heat.
   %
   %#codegen

   persistent cv_liq emiss SB
   if isempty(cv_liq)
      [cv_liq, SB] = icemodel.physicalConstant('cv_liq', 'SB');
      emiss = icemodel.parameterLookup('emiss');
   end

   % Diagnose all fluxes and compute the residual.
   [Qe, Qh] = icemodel.surface.diagnose_turbulent_heat_fluxes(T_sfc, tair, ...
      wspd, psfc, ea_atm, De, br_coefs, ro_sfc, snow_depth, roL, liqflag, ...
      opts);

   residual = chi * Qsi * (1.0 - albedo) + emiss * (Qli - SB * T_sfc ^ 4) ...
      + Qh + Qe + Qc + icemodel.surface.advective_heat_flux(ppt, tppt, cv_liq);

   % For comparison with evaluate_surface_energy_balance:
   % Qle = LONGOUT(T_sfc, emiss, SB);
   % Qa = icemodel.surface.advective_heat_flux(ppt, tppt, cv_liq);
   %
   % residual = chi * Qsi * (1.0 - albedo) + ...
   %    emiss * Qli + Qle + Qh + Qe + Qc + Qa;
end
