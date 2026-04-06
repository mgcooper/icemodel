function [Qe, Qh, diag] = turbulent_heat_flux(T_sfc, tair, wspd, psfc, ...
      ea_atm, De, br_coefs, ro_sfc, snow_depth, roL, liqflag, opts)
   %TURBULENT_HEAT_FLUX Dispatch the configured turbulent-flux scheme.
   %
   %  [Qe, Qh] = icemodel.surface.turbulent_heat_flux(...)
   %  [Qe, Qh, diag] = icemodel.surface.turbulent_heat_flux(...)
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
               icemodel.surface.turbulent_heat_flux_bulk_richardson( ...
               T_sfc, tair, wspd, psfc, ea_atm, De, br_coefs, roL, liqflag, ...
               z0_bulk);
         else
            [Qe, Qh] = icemodel.surface.turbulent_heat_flux_bulk_richardson( ...
               T_sfc, tair, wspd, psfc, ea_atm, De, br_coefs, roL, liqflag, ...
               z0_bulk);
         end

      case 'bulk_mo'
         if nargout > 2
            [Qe, Qh, diag] = icemodel.surface.turbulent_heat_flux_bulk_mo( ...
               T_sfc, tair, wspd, psfc, ea_atm, ro_sfc, snow_depth, roL, ...
               liqflag, opts);
         else
            [Qe, Qh] = icemodel.surface.turbulent_heat_flux_bulk_mo( ...
               T_sfc, tair, wspd, psfc, ea_atm, ro_sfc, snow_depth, roL, ...
               liqflag, opts);
         end

      otherwise
         error('icemodel:configureRun:unknownTurbulentFluxScheme', ...
            'Unrecognized turbulent flux scheme: %s', ...
            opts.turbulent_flux_scheme);
   end
end
