function [T_sfc, Qe, Qh, Qc, Qa, Qle, balance, diag_turbulent] = ...
      surface_energy_balance_terms(T_sfc, tair, Qsi, Qli, albedo, ...
      wspd, ppt, tppt, psfc, De, ea_atm, T, k_eff, dz, roL, chi, br_coefs, ...
      liqflag, ro_sfc, snow_depth, opts)
   %SURFACE_ENERGY_BALANCE_TERMS Evaluate the full SEB term set at T_sfc.
   %
   %  [T_sfc, Qe, Qh, Qc, Qa, Qle, balance] = ...
   %     icemodel.surface.surface_energy_balance_terms(...)
   %
   % The returned terms are:
   %   Qe  turbulent latent heat flux                     [W m^-2]
   %   Qh  turbulent sensible heat flux                  [W m^-2]
   %   Qc  conductive heat flux into the surface         [W m^-2]
   %   Qa  precipitation-advected heat flux              [W m^-2]
   %   Qle emitted longwave heat flux                    [W m^-2]
   %   balance  surface energy-balance residual with Qm=0 [W m^-2]
   %
   % This helper is the canonical "evaluate the surface budget at T_sfc"
   % contract for the touched SEB stack.

   %#codegen

   persistent cv_liq emiss SB
   if isempty(cv_liq)
      [cv_liq, SB] = icemodel.physicalConstant('cv_liq', 'SB');
      emiss = icemodel.parameterLookup('emiss');
   end

   % Turbulent sensible and latent heat exchange.
   if nargout > 7
      [Qe, Qh, diag_turbulent] = ...
         icemodel.surface.diagnose_turbulent_heat_fluxes( ...
         T_sfc, tair, wspd, psfc, ea_atm, De, br_coefs, ro_sfc, snow_depth, ...
         roL, liqflag, opts);
   else
      [Qe, Qh] = icemodel.surface.diagnose_turbulent_heat_fluxes( ...
         T_sfc, tair, wspd, psfc, ea_atm, De, br_coefs, ro_sfc, snow_depth, ...
         roL, liqflag, opts);
      diag_turbulent = struct([]);
   end

   % Surface longwave emission, subsurface conduction, and precipitation heat.
   Qle = LONGOUT(T_sfc, emiss, SB);
   Qc = CONDUCT(k_eff, T, dz, T_sfc);
   Qa = QADVECT(ppt, tppt, cv_liq);

   % The residual is the energy available for melt when Qm is set to zero.
   balance = ENBAL(chi, albedo, Qsi, Qli, Qle, Qh, Qe, Qc, Qa, 0.0);
end
