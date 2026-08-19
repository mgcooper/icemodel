function [d_pevp, pevp, Qe, T_sfc_phys] = ...
      potential_surface_vapor_demand(T_sfc, tair, wspd, psfc, ea_atm, ...
      ro_atm, cv_atm, nu_air, H_h, H_e, hv_atm, br_coefs, liqflag, ...
      f_ice, f_liq, dt, dz, snow_depth, opts)
   %POTENTIAL_SURFACE_VAPOR_DEMAND Diagnose top-cell vapor energy demand.
   %
   %  [d_pevp, pevp, Qe, T_sfc_phys] = ...
   %     icemodel.surface.potential_surface_vapor_demand(...)
   %
   % The model solver can carry a surface temperature above Tf internally.
   % Diagnose the turbulent latent heat flux at the physical surface
   % temperature T_SFC_PHYS = min(T_SFC, Tf), then convert it to the
   % liquid-water-equivalent demand D_PEVP.
   %
   % See also: icemodel.kernels.potential_surface_vapor_demand,
   %  icemodel.surface.potential_surface_vapor_exchange
   %
   %#codegen

   % Surface density belongs to the turbulent heat-flux diagnosis.
   ro_sfc = icemodel.surface.surface_bulk_density(f_ice, f_liq);

   % Use the physical surface temperature for both Qe and its demand.
   T_sfc_phys = icemodel.surface.physical_surface_temperature(T_sfc);
   [Qe, ~] = icemodel.surface.diagnose_turbulent_heat_fluxes( ...
      T_sfc_phys, tair, wspd, psfc, ea_atm, ro_atm, cv_atm, nu_air, ...
      H_h, H_e, hv_atm, br_coefs, liqflag, ro_sfc, snow_depth, opts);

   [d_pevp, pevp] = ...
      icemodel.kernels.potential_surface_vapor_demand(Qe, dt, dz);
end
