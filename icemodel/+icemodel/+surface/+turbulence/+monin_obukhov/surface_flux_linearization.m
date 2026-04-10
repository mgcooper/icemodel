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

   persistent h
   if isempty(h)
      h = 1e-10;
   end

   % Surface heat flux.
   Q_sfc = surface_flux(T_sfc, tair, Qsi, Qli, albedo, wspd, ppt, tppt, ...
      psfc, De, ea_atm, chi, roL, liqflag, ro_sfc, snow_depth, opts);

   % Surface heat flux complex-step perturbation.
   Q_sfc_step = surface_flux(T_sfc + 1i * h, tair, Qsi, Qli, albedo, wspd, ...
      ppt, tppt, psfc, De, ea_atm, chi, roL, liqflag, ro_sfc, snow_depth, opts);

   % Linearization.
   Fp = imag(Q_sfc_step) / h;
   Fc = Q_sfc - Fp * T_sfc;

   % Parse outputs.
   if nargout > 2
      [~, diag_thf] = surface_flux(T_sfc, tair, Qsi, Qli, albedo, wspd, ppt, ...
         tppt, psfc, De, ea_atm, chi, roL, liqflag, ro_sfc, snow_depth, opts);
      diag = struct( ...
         'scheme', 'monin_obukhov', ...
         'q_surface', Q_sfc, ...
         'dq_surface_dTs', Fp, ...
         'thf', diag_thf);
   end
end

function [Q_sfc, diag_thf] = surface_flux(T_sfc, tair, Qsi, Qli, albedo, ...
      wspd, ppt, tppt, psfc, De, ea_atm, chi, roL, liqflag, ro_sfc, ...
      snow_depth, opts)
   %SURFACE_FLUX Evaluate the non-conductive surface flux at T_sfc.

   % Load physical constants and parameters.
   persistent cv_liq
   if isempty(cv_liq)
      cv_liq = icemodel.physicalConstant('cv_liq');
   end

   % Sensible and latent heat fluxes.
   % monin_obukhov does not use br_coefs; pass [] as placeholder for the dispatcher.
   if nargout > 1
      [Qe, Qh, diag_thf] = ...
         icemodel.surface.diagnose_turbulent_heat_fluxes(T_sfc, tair, wspd, ...
         psfc, ea_atm, De, [], ro_sfc, snow_depth, roL, liqflag, opts);
   else
      [Qe, Qh] = icemodel.surface.diagnose_turbulent_heat_fluxes(T_sfc, ...
         tair, wspd, psfc, ea_atm, De, [], ro_sfc, snow_depth, roL, ...
         liqflag, opts);
      diag_thf = struct([]);
   end

   % Precipitation advected heat flux.
   Qa = icemodel.surface.advective_heat_flux(ppt, tppt, cv_liq);

   % Net non-conductive heat flux (Qc = 0: handled by the Robin interior solve).
   Qsn = icemodel.surface.net_shortwave_radiation(Qsi, albedo, chi);
   Qln = icemodel.surface.net_longwave_radiation(T_sfc, Qli);
   Q_sfc = icemodel.surface.evaluate_surface_energy_balance( ...
      Qsn, Qln, Qh, Qe, 0.0, Qa, 0.0);
end
