function [Qe, Qh, diag] = turbulent_heat_flux(Ta, Ts, wspd, Pa, De, ea, ...
      cv_air, roL, scoef, liqflag, ro_sfc, snow_depth, opts)
   %TURBULENT_HEAT_FLUX Dispatch the configured turbulent-flux scheme.
   %
   %  [Qe, Qh] = icemodel.surface.turbulent_heat_flux(...)
   %  [Qe, Qh, diag] = icemodel.surface.turbulent_heat_flux(...)
   %
   % The production contract is explicit: callers must provide the resolved
   % surface state and opts struct. Optional diagnostics are constructed only
   % when the third output is requested, primarily for the test suite or
   % interactive debugging / inspection.
   %
   %#codegen

   scheme = lower(char(opts.turbulent_flux_scheme));
   switch scheme
      case 'bulk_richardson'
         if nargout > 2
            [Qe, Qh, diag] = ...
               icemodel.surface.turbulent_heat_flux_bulk_richardson( ...
               Ta, Ts, wspd, Pa, De, ea, cv_air, roL, scoef, liqflag, ...
               opts.z0_bulk);
         else
            [Qe, Qh] = icemodel.surface.turbulent_heat_flux_bulk_richardson( ...
               Ta, Ts, wspd, Pa, De, ea, cv_air, roL, scoef, liqflag, ...
               opts.z0_bulk);
         end

      case 'bulk_mo'
         if nargout > 2
            [Qe, Qh, diag] = icemodel.surface.turbulent_heat_flux_bulk_mo( ...
               Ta, Ts, wspd, Pa, ea, cv_air, roL, liqflag, ro_sfc, ...
               snow_depth, opts);
         else
            [Qe, Qh] = icemodel.surface.turbulent_heat_flux_bulk_mo( ...
               Ta, Ts, wspd, Pa, ea, cv_air, roL, liqflag, ro_sfc, ...
               snow_depth, opts);
         end

      otherwise
         error('icemodel:configureRun:unknownTurbulentFluxScheme', ...
            'Unrecognized turbulent flux scheme: %s', opts.turbulent_flux_scheme);
   end
end
