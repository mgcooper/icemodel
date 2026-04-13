function [Qe, Qh, diag] = diagnose_turbulent_heat_fluxes(T_sfc, tair, wspd, ...
      psfc, ea_atm, De, br_coefs, ro_sfc, snow_depth, ro_air_Lv, liqflag, opts)
   %DIAGNOSE_TURBULENT_HEAT_FLUXES Dispatch the configured THF scheme.
   %
   %  [Qe, Qh] = icemodel.surface.diagnose_turbulent_heat_fluxes(...)
   %  [Qe, Qh, diag] = icemodel.surface.diagnose_turbulent_heat_fluxes(...)
   %
   % The production contract requires callers to provide the resolved surface
   % state and opts struct. Optional diagnostics are constructed only when the
   % third output is requested, primarily for the test suite or interactive
   % debugging / inspection.
   %
   %#codegen

   z0_bulk = opts.z0_bulk;

   scheme = lower(char(opts.turbulent_flux_scheme));
   switch scheme
      case 'bulk_richardson'
         if nargout > 2
            [Qe, Qh, diag] = ...
               icemodel.surface.turbulence.bulk_richardson.turbulent_heat_flux( ...
               T_sfc, tair, wspd, psfc, ea_atm, De, br_coefs, ro_air_Lv, ...
               liqflag, z0_bulk);
         else
            [Qe, Qh] = ...
               icemodel.surface.turbulence.bulk_richardson.turbulent_heat_flux( ...
               T_sfc, tair, wspd, psfc, ea_atm, De, br_coefs, ro_air_Lv, ...
               liqflag, z0_bulk);
         end

      case 'monin_obukhov'
         if nargout > 2
            [Qe, Qh, diag] = ...
               icemodel.surface.turbulence.monin_obukhov.turbulent_heat_flux( ...
               T_sfc, tair, wspd, psfc, ea_atm, ro_sfc, snow_depth, ro_air_Lv, ...
               liqflag, opts);
         else
            [Qe, Qh] = ...
               icemodel.surface.turbulence.monin_obukhov.turbulent_heat_flux( ...
               T_sfc, tair, wspd, psfc, ea_atm, ro_sfc, snow_depth, ro_air_Lv, ...
               liqflag, opts);
         end

      otherwise
         error('icemodel:configureRun:unknownTurbulentFluxScheme', ...
            'Unrecognized turbulent flux scheme: %s', ...
            opts.turbulent_flux_scheme);
   end
end
