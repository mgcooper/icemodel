function [Fc, Fp] = surface_flux_linearization(T_sfc, tair, Qsi, Qli, albedo, ...
      wspd, ppt, tppt, psfc, De, ea_atm, roL, br_coefs, chi, liqflag, ...
      ro_sfc, snow_depth, opts)
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
         T_sfc, tair, Qsi, Qli, albedo, wspd, ppt, tppt, psfc, De, ea_atm, chi, ...
         roL, liqflag, ro_sfc, snow_depth, opts);
   else
      [Fc, Fp] = SFCFLIN(tair, Qsi, Qli, albedo, wspd, ppt, tppt, psfc, De, ...
         ea_atm, roL, br_coefs, chi, T_sfc, liqflag);
   end
end
