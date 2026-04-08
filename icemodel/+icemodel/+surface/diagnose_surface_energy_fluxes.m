function [Qe, Qh, Qc, Qm, Qf, Qbal, diag_turbulent] = ...
      diagnose_surface_energy_fluxes(T_ice, xT_sfc, tair, Qsi, Qli, albedo, ...
      wspd, ppt, tppt, psfc, De, ea_atm, k_eff, dz, roL, chi, br_coefs, ...
      ro_sfc, snow_depth, liqflag, opts)
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

   % Diagnose the physical surface temperature.
   T_sfc_phys = icemodel.surface.physical_surface_temperature(xT_sfc);

   % Diagnose the surface energy fluxes at T_sfc_phys.
   if nargout > 6
      [T_sfc_phys, Qe, Qh, Qc, Qa, Qle, ~, diag_turbulent] = ...
         icemodel.surface.surface_energy_balance_terms(T_sfc_phys, tair, ...
         Qsi, Qli, albedo, wspd, ppt, tppt, psfc, De, ea_atm, T_ice, k_eff, ...
         dz, roL, chi, br_coefs, liqflag, ro_sfc, snow_depth, opts);
   else
      [T_sfc_phys, Qe, Qh, Qc, Qa, Qle, ~] = ...
         icemodel.surface.surface_energy_balance_terms(T_sfc_phys, tair, ...
         Qsi, Qli, albedo, wspd, ppt, tppt, psfc, De, ea_atm, T_ice, k_eff, ...
         dz, roL, chi, br_coefs, liqflag, ro_sfc, snow_depth, opts);
      diag_turbulent = struct([]);
   end

   % Diagnose the energy for melting or freezing. If T_sfc_phys = Tf, Qm is the
   % energy balance residual and Qf = 0. If T_sfc_phys < Tf, Qf is the energy
   % balance residual
   [Qm, Qf] = MFENERGY(T_sfc_phys, chi, albedo, Qsi, Qli, Qle, Qh, Qe, Qc, ...
      Qa, tair, wspd, psfc, ppt, tppt, ea_atm, br_coefs, De, T_ice, k_eff, ...
      dz, ro_sfc, snow_depth, roL, opts);

   % Evaluate the energy balance residual given the diagnosed Qm.
   Qbal = ENBAL(chi, albedo, Qsi, Qli, Qle, Qh, Qe, Qc, Qa, Qm);
end
