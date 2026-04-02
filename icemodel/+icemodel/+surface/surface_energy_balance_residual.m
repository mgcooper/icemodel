function residual = surface_energy_balance_residual(Ts, Ta, Qsi, Qli, ...
      albedo, wspd, ppt, tppt, Pa, De, ea, chi, roL, scoef, Qc, ...
      liqflag, ro_sfc, snow_depth, opts)
   %SURFACE_ENERGY_BALANCE_RESIDUAL Return the SEB residual at Ts.
   %
   %  residual = icemodel.surface.surface_energy_balance_residual(...)
   %
   % This is the canonical nonlinear surface residual used by the touched
   % SEB solvers. It includes absorbed shortwave, net longwave, conductive
   % heat, sensible heat, latent heat, and precipitation-advected heat.
   %
   %#codegen

   persistent cv_air cv_liq emiss SB
   if isempty(cv_air)
      [cv_air, cv_liq, SB] = icemodel.physicalConstant( ...
         'cv_air', 'cv_liq', 'SB');
      emiss = icemodel.parameterLookup('emiss');
   end

   [Qe, Qh] = icemodel.surface.turbulent_heat_flux(Ta, Ts, wspd, Pa, De, ...
      ea, cv_air, roL, scoef, liqflag, ro_sfc, snow_depth, opts);

   residual = chi * (1.0 - albedo) * Qsi + emiss * (Qli - SB * Ts ^ 4) ...
      + Qc + Qh + Qe + QADVECT(ppt, tppt, cv_liq);
end
