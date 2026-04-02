function [Fc, Fp] = surface_flux_linearization(Ts, Ta, Qsi, Qli, albedo, ...
      wspd, ppt, tppt, Pa, De, ea, cv_air, cv_liq, emiss, SB, roL, ...
      scoef, chi, Tf, liqflag, ro_sfc, snow_depth, opts)
   %SURFACE_FLUX_LINEARIZATION Linearize the non-conductive surface flux.
   %
   %  [Fc, Fp] = icemodel.surface.surface_flux_linearization(...)
   %
   % This is the canonical Robin-boundary linearization contract for the
   % touched SEB stack. The legacy bulk-Richardson path keeps the analytic
   % `SFCFLIN` linearization, while `bulk_mo` uses a dedicated numerical
   % linearization of the nonlinear atmospheric surface-flux closure.
   %
   %#codegen

   if strcmp(opts.turbulent_flux_scheme, 'bulk_mo')
      [Fc, Fp] = icemodel.surface.surface_flux_linearization_bulk_mo( ...
         Ts, Ta, Qsi, Qli, albedo, wspd, ppt, tppt, Pa, De, ea, chi, ...
         roL, liqflag, ro_sfc, snow_depth, opts);
   else
      [Fc, Fp] = SFCFLIN(Ta, Qsi, Qli, albedo, wspd, ppt, tppt, Pa, De, ...
         ea, cv_air, cv_liq, emiss, SB, roL, scoef, chi, Tf, Ts, liqflag);
   end
end
