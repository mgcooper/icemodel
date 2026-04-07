function [Qe, Qh, Qc, Qm, Qf, Qbal, diag_turbulent] = ...
      diagnose_surface_energy_fluxes(T, xTs, tair, ...
      Qsi, Qli, albedo, wspd, ppt, tppt, psfc, De, ea_atm, k_eff, dz, roL, chi, ...
      br_coefs, ro_sfc, snow_depth, liqflag, opts)
   %DIAGNOSE_SURFACE_ENERGY_FLUXES Diagnose the component surface fluxes.
   %
   %  [Qe, Qh, Qc, Qm, Qf, Qbal] = ...
   %     icemodel.surface.diagnose_surface_energy_fluxes(...)
   %  [Qe, Qh, Qc, Qm, Qf, Qbal, diag_turbulent] = ...
   %     icemodel.surface.diagnose_surface_energy_fluxes(...)
   %
   %
   % ea_atm  - atmospheric vapor pressure from relative humidity data.
   % es_sfc  - surface saturation water vapor pressure.
   % Qe  - turbulent latent heat flux.
   % Qh  - turbulent sensible heat flux.
   % Qle - surface emitted longwave heat flux.
   % Qc  - conductive heat flux into the surface.
   % Qm  - the energy flux available for melting or freezing.
   % Qbal - the net diagnosed surface energy-balance residual.
   %
   %#codegen

   T_sfc = icemodel.surface.physical_surface_temperature(xTs);

   if nargout > 6
      [T_sfc, Qe, Qh, Qc, Qa, Qle, balance, diag_turbulent] = ...
         icemodel.surface.surface_energy_balance_terms(T_sfc, tair, Qsi, ...
         Qli, albedo, wspd, ppt, tppt, psfc, De, ea_atm, T, k_eff, dz, ...
         roL, chi, br_coefs, liqflag, ro_sfc, snow_depth, opts);
   else
      [T_sfc, Qe, Qh, Qc, Qa, Qle, balance] = ...
         icemodel.surface.surface_energy_balance_terms(T_sfc, tair, Qsi, ...
         Qli, albedo, wspd, ppt, tppt, psfc, De, ea_atm, T, k_eff, dz, ...
         roL, chi, br_coefs, liqflag, ro_sfc, snow_depth, opts);
      diag_turbulent = struct([]);
   end

   [Qm, Qf] = MFENERGY(T_sfc, chi, albedo, Qsi, Qli, Qle, Qh, Qe, Qc, Qa, ...
      tair, wspd, psfc, ppt, tppt, ea_atm, br_coefs, De, T, k_eff, dz, ...
      ro_sfc, snow_depth, roL, opts);

   Qbal = ENBAL(chi, albedo, Qsi, Qli, Qle, Qh, Qe, Qc, Qa, Qm);
end
