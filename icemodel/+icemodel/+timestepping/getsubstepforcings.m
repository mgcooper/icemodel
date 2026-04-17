function [tair, swd, lwd, albedo, wspd, rh, psfc, ppt, tppt, De, ea] ...
      = getsubstepforcings(opts, metstep, dt_sum, dt_full, liqflag, tair_v, ...
      swd_v, lwd_v, albedo_v, wspd_v, rh_v, psfc_v, ppt_v, tppt_v, De_v)
   %GETSUBSTEPFORCINGS Return scalar forcing values for the current substep.
   %
   % Linearly interpolates between the current (metstep) and next (metstep+1)
   % forcing values based on the fractional substep progress dt_sum/dt_full.
   % When opts.met_substep_interp is false or dt_sum == 0, returns the current
   % timestep value with no interpolation.
   %
   % NOTE: This function was written for experiments with substep-averaged flux
   % diagnostics under sudden forcing transitions (e.g. no-leap calendar
   % boundaries). It is not currently called in production runs. To activate,
   % interpolate _s variants into solve_surface_energy_balance and recompute ea
   % before diagnose_surface_energy_balance.
   %
   %#codegen

   i1 = metstep;
   i2 = min(metstep + 1, numel(tair_v));

   w = 0.0;
   if opts.met_substep_interp && i2 > i1 && dt_full > 0
      w = min(max(dt_sum / dt_full, 0.0), 1.0);
   end

   tair   = (1.0 - w) * tair_v(i1)   + w * tair_v(i2);
   swd    = (1.0 - w) * swd_v(i1)    + w * swd_v(i2);
   lwd    = (1.0 - w) * lwd_v(i1)    + w * lwd_v(i2);
   albedo = (1.0 - w) * albedo_v(i1) + w * albedo_v(i2);
   wspd   = (1.0 - w) * wspd_v(i1)   + w * wspd_v(i2);
   rh     = (1.0 - w) * rh_v(i1)     + w * rh_v(i2);
   psfc   = (1.0 - w) * psfc_v(i1)   + w * psfc_v(i2);
   ppt    = (1.0 - w) * ppt_v(i1)    + w * ppt_v(i2);
   tppt   = (1.0 - w) * tppt_v(i1)   + w * tppt_v(i2);
   De     = (1.0 - w) * De_v(i1)     + w * De_v(i2);

   ea = icemodel.surface.atmospheric_vapor_pressure(tair, rh, liqflag);
end
