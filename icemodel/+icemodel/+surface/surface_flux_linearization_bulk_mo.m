function [Fc, Fp, diag] = surface_flux_linearization_bulk_mo(Ts, Ta, ...
      Qsi, Qli, albedo, wspd, ppt, tppt, Pa, De, ea, chi, roL, liqflag, ...
      ro_sfc, snow_depth, opts)
   %SURFACE_FLUX_LINEARIZATION_BULK_MO Linearize the bulk-MO surface flux.
   %
   %  [Fc, Fp] = icemodel.surface.surface_flux_linearization_bulk_mo(...)
   %  [Fc, Fp, diag] = icemodel.surface.surface_flux_linearization_bulk_mo(...)
   %
   % The Robin coupler expects a linear boundary-flux form:
   %   q_surface(Ts) ≈ Fc + Fp * Ts
   %
   % This helper evaluates the current surface flux state, then computes the
   % local derivative using a complex-step perturbation. It intentionally
   % linearizes only the atmospheric surface flux (shortwave, longwave,
   % sensible, latent, and precipitation advection) so the subsurface
   % conductive term remains in the ICEENBAL / Robin interior solve.
   %
   %#codegen

   persistent h
   if isempty(h)
      h = 1e-10;
   end

   q_surface = surface_flux(Ts, Ta, Qsi, Qli, albedo, wspd, ppt, tppt, ...
      Pa, De, ea, chi, roL, liqflag, ro_sfc, snow_depth, opts);
   q_surface_step = surface_flux(Ts + 1i * h, Ta, Qsi, Qli, albedo, wspd, ...
      ppt, tppt, Pa, De, ea, chi, roL, liqflag, ro_sfc, snow_depth, opts);

   Fp = imag(q_surface_step) / h;
   Fc = q_surface - Fp * Ts;

   if nargout > 2
      [~, diag_thf] = surface_flux(Ts, Ta, Qsi, Qli, albedo, wspd, ppt, ...
         tppt, Pa, De, ea, chi, roL, liqflag, ro_sfc, snow_depth, opts);
      diag = struct( ...
         'scheme', 'bulk_mo', ...
         'q_surface', q_surface, ...
         'dq_surface_dTs', Fp, ...
         'thf', diag_thf);
   end
end

function [q_surface, diag_thf] = surface_flux(Ts, Ta, Qsi, Qli, albedo, wspd, ppt, ...
      tppt, Pa, De, ea, chi, roL, liqflag, ro_sfc, snow_depth, opts)
   %SURFACE_FLUX Evaluate the non-conductive surface flux at Ts.

   persistent cv_air cv_liq emiss SB
   if isempty(cv_liq)
      [cv_air, cv_liq, SB] = icemodel.physicalConstant( ...
         'cv_air', 'cv_liq', 'SB');
      emiss = icemodel.parameterLookup('emiss');
   end

   if nargout > 1
      [Qe, Qh, diag_thf] = icemodel.surface.turbulent_heat_flux(Ta, Ts, ...
         wspd, Pa, De, ea, cv_air, roL, 0.0, liqflag, ro_sfc, ...
         snow_depth, opts);
   else
      [Qe, Qh] = icemodel.surface.turbulent_heat_flux(Ta, Ts, wspd, Pa, ...
         De, ea, cv_air, roL, 0.0, liqflag, ro_sfc, snow_depth, opts);
      diag_thf = struct([]);
   end

   Qle = LONGOUT(Ts, emiss, SB);
   Qa = QADVECT(ppt, tppt, cv_liq);
   q_surface = ENBAL(albedo, emiss, chi, Qsi, Qli, Qle, Qh, Qe, 0.0, Qa, 0.0);
end
