function [Qe, Qh, Qc, Qa, Qle, Qsn, Qln, diag_turbulent] = ...
      surface_energy_balance_terms(T_sfc, tair, Qsi, Qli, albedo, wspd, ppt, ...
      tppt, psfc, De, ea_atm, br_coefs, roL, liqflag, chi, T_ice, k_eff, dz, ...
      ro_sfc, snow_depth, opts)
   %SURFACE_ENERGY_BALANCE_TERMS Evaluate the SEB term set at T_sfc.
   %
   %  [Qe, Qh, Qc, Qa, Qle] = ...
   %     icemodel.surface.surface_energy_balance_terms(...)
   %  [Qe, Qh, Qc, Qa, Qle, Qsn, Qln] = ...
   %     icemodel.surface.surface_energy_balance_terms(...)
   %  [Qe, Qh, Qc, Qa, Qle, Qsn, Qln, diag_turbulent] = ...
   %     icemodel.surface.surface_energy_balance_terms(...)
   %
   % The returned terms are:
   %   Qe  turbulent latent heat flux                    [W m^-2]
   %   Qh  turbulent sensible heat flux                  [W m^-2]
   %   Qc  conductive heat flux into the surface         [W m^-2]
   %   Qa  precipitation-advected heat flux              [W m^-2]
   %   Qle emitted longwave heat flux                    [W m^-2]
   %
   % Optional outputs (guarded by nargout to avoid extra computation when
   % only the T_sfc-dependent flux terms above are needed):
   %   Qsn net absorbed shortwave flux = chi*Qsi*(1-albedo) [W m^-2]
   %   Qln net longwave flux = emiss*Qli + Qle              [W m^-2]
   %   diag_turbulent struct from the turbulent flux scheme
   %
   % This helper is the canonical "evaluate the surface budget at T_sfc"
   % contract for the SEB diagnostic stack.
   %
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
   Qc = icemodel.surface.conductive_heat_flux(k_eff, T_ice, dz, T_sfc);
   Qa = icemodel.surface.advective_heat_flux(ppt, tppt, cv_liq);

   % Net shortwave and net longwave (optional; guards avoid extra computation).
   if nargout > 5
      Qsn = chi * Qsi * (1.0 - albedo);
      Qln = emiss * Qli + Qle;
   end
end
