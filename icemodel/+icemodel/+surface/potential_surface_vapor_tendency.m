function [d_pevp, pevp, Qe, Ts_phys] = potential_surface_vapor_tendency( ...
      Ts, Tf, Ta, wspd, Pa, De, ea, cv_air, roL, scoef, liqflag, ro_sfc, ...
      snow_depth, opts, Lv, ro_liq, dt, dz)
   %POTENTIAL_SURFACE_VAPOR_TENDENCY Evaluate the top-layer vapor tendency.
   %
   %  [d_pevp, pevp, Qe, Ts_phys] = ...
   %     icemodel.surface.potential_surface_vapor_tendency(...)
   %
   % This helper is the physical-flux contract used by icemodel mass
   % balance updates. Ts may be an internal solver temperature that exceeds
   % Tf; the turbulent latent-heat flux and derived vapor tendency must use
   % the physical diagnosed surface temperature Ts_phys = min(Ts, Tf).
   %
   %#codegen

   Ts_phys = MELTTEMP(Ts, Tf);
   [Qe, ~] = icemodel.surface.turbulent_heat_flux(Ta, Ts_phys, wspd, Pa, ...
      De, ea, cv_air, roL, scoef, liqflag, ro_sfc, snow_depth, opts);
   [d_pevp, pevp] = PEVAP(Qe, Lv, ro_liq, dt, dz);
end
