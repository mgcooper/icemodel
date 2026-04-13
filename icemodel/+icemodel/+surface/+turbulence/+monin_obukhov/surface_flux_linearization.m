function [Fc, Fp, diag] = surface_flux_linearization(T_sfc, tair, Qsi, Qli, ...
      albedo, wspd, ppt, tppt, psfc, De, ea_atm, chi, roL, liqflag, ro_sfc, ...
      snow_depth, opts)
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
   %#codegen

   if nargout > 2
      [~, ~, diag_thf] = ...
         icemodel.surface.diagnose_turbulent_heat_fluxes(T_sfc, tair, wspd, ...
         psfc, ea_atm, De, [], ro_sfc, snow_depth, roL, liqflag, opts);
   else
      diag_thf = struct([]);
   end

   [Q_sfc, Fp] = icemodel.surface.numerical_surface_flux_linearization( ...
      T_sfc, tair, Qsi, Qli, albedo, wspd, ppt, tppt, psfc, De, ea_atm, ...
      [], roL, liqflag, chi, 0.0, 0.0, 1.0, ro_sfc, snow_depth, opts);
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
