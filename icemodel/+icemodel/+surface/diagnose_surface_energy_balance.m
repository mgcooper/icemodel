function [Qe, Qh, Qc, Qsn, Qln, Qa, Qm, Qf, Qbal] = ...
      diagnose_surface_energy_balance(T_sfc, tair, Qsi, Qli, albedo, wspd, ...
      ppt, tppt, psfc, De, ea_atm, T_ice, k_eff, dz, roL, chi, br_coefs, ...
      ro_sfc, snow_depth, liqflag, opts)
   %DIAGNOSE_SURFACE_ENERGY_BALANCE Diagnose the full surface energy budget.
   %
   %  [Qe, Qh, Qc, Qsn, Qln, Qa, Qm, Qf, Qbal] = ...
   %     icemodel.surface.diagnose_surface_energy_balance(...)
   %
   % This function handles the T_sfc_phys contract and returns the full SEB
   % term set plus melt/freezing diagnostics. It wraps the canonical SEB
   % diagnostic stack:
   %  icemodel.surface.surface_energy_balance_terms
   %  icemodel.surface.diagnose_melt_freeze_energy
   %  icemodel.surface.evaluate_surface_energy_balance
   %
   % Diagnostic variable definitions:
   %  ea_atm  - atmospheric vapor pressure from relative humidity data.
   %  Qe  - turbulent latent heat flux.
   %  Qh  - turbulent sensible heat flux.
   %  Qc  - conductive heat flux into the surface.
   %  Qsn - net shortwave radiation.
   %  Qln - net longwave radiation.
   %  Qa  - precipitation-advected heat flux.
   %  Qm  - the energy flux available for melting.
   %  Qf  - the energy flux required to warm the surface to melting.
   %  Qbal - the net diagnosed surface energy-balance residual.
   %
   %#codegen

   % Diagnose the physical surface temperature.
   T_sfc_phys = icemodel.surface.physical_surface_temperature(T_sfc);

   % Diagnose the surface energy terms at T_sfc_phys.
   [Qe, Qh, Qc, Qa, Qsn, Qln] = ...
      icemodel.surface.surface_energy_balance_terms(T_sfc_phys, tair, ...
      Qsi, Qli, albedo, wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, ...
      roL, liqflag, chi, T_ice, k_eff, dz, ro_sfc, snow_depth, opts);

   % Diagnose the energy for melting or freezing. If T_sfc_phys = Tf, Qm is the
   % energy balance residual and Qf = 0. If T_sfc_phys < Tf, Qf is the energy
   % balance residual.
   [Qm, Qf] = icemodel.surface.diagnose_melt_freeze_energy( ...
      T_sfc_phys, Qsn, Qln, Qh, Qe, Qc, Qa, tair, Qsi, Qli, albedo, wspd, ...
      psfc, ppt, tppt, De, ea_atm, chi, br_coefs, roL, T_ice, k_eff, dz, ro_sfc, ...
      snow_depth, opts);

   % Evaluate the energy balance residual given the diagnosed Qm.
   Qbal = icemodel.surface.evaluate_surface_energy_balance( ...
      Qsn, Qln, Qh, Qe, Qc, Qa, Qm);
end
