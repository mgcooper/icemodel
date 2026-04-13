function [Fc, Fp] = surface_flux_linearization(T_sfc, tair, Qsi, Qli, albedo, ...
      wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, ro_air_Lv, liqflag, chi, ...
      ro_sfc, snow_depth, opts)
   %SURFACE_FLUX_LINEARIZATION Linearize the non-conductive surface flux.
   %
   %  [Fc, Fp] = icemodel.surface.surface_flux_linearization(...)
   %
   % This is the canonical Robin-boundary linearization contract for the
   % touched SEB stack. The bulk-Richardson path keeps an analytic
   % linearization, while `monin_obukhov` uses a dedicated numerical
   % linearization of the nonlinear atmospheric surface-flux closure.
   %
   % See also: icemodel.surface.numerical_surface_flux
   %
   %#codegen

   if strcmp(opts.turbulent_flux_scheme, 'monin_obukhov')
      [Fc, Fp] = ...
         icemodel.surface.turbulence.monin_obukhov.surface_flux_linearization( ...
         T_sfc, tair, Qsi, Qli, albedo, wspd, ppt, tppt, psfc, De, ea_atm, chi, ...
         ro_air_Lv, liqflag, ro_sfc, snow_depth, opts);

   elseif strcmp(opts.turbulent_flux_scheme, 'bulk_richardson')
      [Fc, Fp] = ...
         icemodel.surface.turbulence.bulk_richardson.surface_flux_linearization( ...
         T_sfc, tair, Qsi, Qli, albedo, wspd, ppt, tppt, psfc, De, ea_atm, ...
         br_coefs, ro_air_Lv, liqflag, chi);
   else
      error('icemodel:surface_flux_linearization:unknownTurbulentFluxScheme', ...
         'Unrecognized turbulent flux scheme: %s', opts.turbulent_flux_scheme);
   end
end
