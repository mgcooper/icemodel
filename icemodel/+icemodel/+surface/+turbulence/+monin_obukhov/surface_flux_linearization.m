function [Fc, Fp, diag] = surface_flux_linearization(T_sfc, tair, Qsi, ...
      Qli, albedo, wspd, ppt, tppt, psfc, ea_atm, ro_atm, cv_atm, ...
      nu_air, H_h, H_e, hv_atm, liqflag, chi, ro_sfc, snow_depth, opts)
   %SURFACE_FLUX_LINEARIZATION Linearize the bulk-MO surface flux.
   %
   %  [Fc, Fp] = ...
   %     icemodel.surface.turbulence.monin_obukhov.surface_flux_linearization(...)
   %  [Fc, Fp, diag] = ...
   %     icemodel.surface.turbulence.monin_obukhov.surface_flux_linearization(...)
   %
   % The Robin coupler expects a linear boundary-flux form:
   %   Q_sfc(T_sfc) ≈ Fc + Fp * T_sfc
   %
   % This helper evaluates the current surface flux state, then computes the
   % local derivative using a complex-step perturbation. It intentionally
   % linearizes only the atmospheric surface flux (shortwave, longwave,
   % sensible, latent, and precipitation advection) so the subsurface
   % conductive term remains in the `icemodel.column.solve_column_enthalpy`
   % Robin interior solve.
   %
   % cv_atm, hv_atm, ro_atm, nu_air are precomputed per forcing timestep and
   % passed in; they are constant across inner solver iterations.
   %
   %#codegen

   % Set bulk-Richardson stability coefficients empty for Monin-Obukhov scheme.
   % Set T_ice, k_eff, dz to dummy values to ignore the conduction term.
   persistent br_coefs T_ice k_eff dz
   if isempty(br_coefs)
      br_coefs = [];
      T_ice = 0.0;
      k_eff = 0.0;
      dz = 1.0;
   end

   if nargout > 2
      [~, ~, diag_thf] = ...
         icemodel.surface.diagnose_turbulent_heat_fluxes(T_sfc, ...
         tair, wspd, psfc, ea_atm, ro_atm, cv_atm, nu_air, H_h, ...
         H_e, hv_atm, br_coefs, liqflag, ro_sfc, snow_depth, opts);
   else
      diag_thf = struct([]);
   end

   % Fp = dQ_sfc_dT_sfc
   [Q_sfc, Fp] = icemodel.surface.numerical_surface_flux(T_sfc, ...
      tair, Qsi, Qli, albedo, wspd, ppt, tppt, psfc, ea_atm, ...
      ro_atm, cv_atm, nu_air, H_h, H_e, hv_atm, br_coefs, liqflag, ...
      chi, T_ice, k_eff, dz, ro_sfc, snow_depth, opts);

   Fc = Q_sfc - Fp * T_sfc;

   % Parse outputs.
   if nargout > 2
      diag = struct( ...
         'scheme', 'monin_obukhov', ...
         'q_surface', Q_sfc, ...
         'dq_surface_dTs', Fp, ...
         'thf', diag_thf);
   end
end
