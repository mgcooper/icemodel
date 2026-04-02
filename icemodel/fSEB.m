function f = fSEB(Ts, Ta, Qsi, Qli, albedo, wspd, ppt, tppt, Pa, De, ...
      ea, chi, roL, scoef, Qc, liqflag, ro_sfc, snow_depth, opts)
   %FSEB Surface energy-balance residual evaluated at Ts.
   %
   % This top-level kernel remains the canonical entry point used by the
   % existing solvers and tests. The touched contract is explicit: callers
   % must provide the surface-state inputs rather than relying on filler
   % defaults for missing roughness or option state.
   %
   %#codegen

   f = icemodel.surface.surface_energy_balance_residual(Ts, Ta, Qsi, Qli, ...
      albedo, wspd, ppt, tppt, Pa, De, ea, chi, roL, scoef, Qc, liqflag, ...
      ro_sfc, snow_depth, opts);
end
