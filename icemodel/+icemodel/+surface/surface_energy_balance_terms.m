function [Qe, Qh, Qc, Qa, Qsn, Qln] = surface_energy_balance_terms(T_sfc, ...
      tair, Qsi, Qli, albedo, wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, ...
      ro_air_Lv, liqflag, chi, T_ice, k_eff, dz, ro_sfc, snow_depth, opts)
   %SURFACE_ENERGY_BALANCE_TERMS Evaluate the SEB term set at T_sfc.
   %
   %  [Qe, Qh, Qc, Qa, Qsn, Qln] = ...
   %     icemodel.surface.surface_energy_balance_terms(...)
   %
   % The returned terms are:
   %   Qe  turbulent latent heat flux                    [W m^-2]
   %   Qh  turbulent sensible heat flux                  [W m^-2]
   %   Qc  conductive heat flux into the surface         [W m^-2]
   %   Qa  precipitation-advected heat flux              [W m^-2]
   %   Qsn net absorbed shortwave flux                   [W m^-2]
   %   Qln net longwave flux                             [W m^-2]
   %
   % This helper is the canonical "evaluate the surface budget at T_sfc"
   % contract for the SEB diagnostic stack.
   %
   %#codegen

   % Pass cv_liq to advective_heat_flux. TODO: For snowfall, pass snow density.
   persistent cv_liq
   if isempty(cv_liq)
      cv_liq = icemodel.physicalConstant('cv_liq');
   end

   % Turbulent sensible and latent heat exchange.
   [Qe, Qh] = icemodel.surface.diagnose_turbulent_heat_fluxes( ...
      T_sfc, tair, wspd, psfc, ea_atm, De, br_coefs, ro_sfc, snow_depth, ...
      ro_air_Lv, liqflag, opts);

   % Surface conduction and precipitation heat.
   Qc = icemodel.surface.conductive_heat_flux(k_eff, T_ice, dz, T_sfc);
   Qa = icemodel.surface.advective_heat_flux(ppt, tppt, cv_liq);

   % Net radiative fluxes.
   Qsn = icemodel.surface.net_shortwave_radiation(Qsi, albedo, chi);
   Qln = icemodel.surface.net_longwave_radiation(T_sfc, Qli);
end
