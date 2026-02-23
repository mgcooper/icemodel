function [tair, swd, lwd, albedo, wspd, rh, psfc, ppt, tppt, De, ea] ...
      = METSUB(opts, metiter, dt_sum, dt_full, liqflag, Tf, tair_v, ...
      swd_v, lwd_v, albedo_v, wspd_v, rh_v, psfc_v, ppt_v, tppt_v, De_v)
   %METSUB Return scalar forcing values for the current substep.
   %
   %#codegen

   % This was abandoned b/c it would require computing substep-averaged flux
   % diagnostics and mainly helped with sudden forcing changes when using a
   % no-leap calendar which will be removed eventually. If using, add a call
   % like this:
   %
   % while dt_sum + TINY < dt_FULL_STEP
   %
   %    [tair_s, swd_s, lwd_s, albedo_s, wspd_s, ~, psfc_s, ppt_s, ...
   %       tppt_s, De_s, ea] = METSUB(opts, metiter, dt_sum, ...
   %       dt_FULL_STEP, liqflag0, Tf, tair, swd, lwd, albedo, wspd, rh, ...
   %       psfc, ppt, tppt, De);
   %
   % then replace the vars passed to SEBSOLVE with the _s variants, and
   % re-compute ea prior to diagnosing fluxes (SEBFLUX)

   i1 = metiter;
   i2 = min(metiter + 1, numel(tair_v));

   w = 0.0;
   if opts.met_substep_interp && i2 > i1 && dt_full > 0
      w = min(max(dt_sum / dt_full, 0.0), 1.0);
   end

   tair = (1.0 - w) * tair_v(i1) + w * tair_v(i2);
   swd = (1.0 - w) * swd_v(i1) + w * swd_v(i2);
   lwd = (1.0 - w) * lwd_v(i1) + w * lwd_v(i2);
   albedo = (1.0 - w) * albedo_v(i1) + w * albedo_v(i2);
   wspd = (1.0 - w) * wspd_v(i1) + w * wspd_v(i2);
   rh = (1.0 - w) * rh_v(i1) + w * rh_v(i2);
   psfc = (1.0 - w) * psfc_v(i1) + w * psfc_v(i2);
   ppt = (1.0 - w) * ppt_v(i1) + w * ppt_v(i2);
   tppt = (1.0 - w) * tppt_v(i1) + w * tppt_v(i2);
   De = (1.0 - w) * De_v(i1) + w * De_v(i2);

   ea = VAPPRESS(tair, Tf, liqflag) * rh / 100;
end
