function [f_ice, f_liq, x_ice, x_liq] = MAXWATER(f_ice, f_liq, ro_ice, ...
      ro_liq, d_liq, d_ice)

   % NOTE: I did not appreciate that fliqmin / fliqmax depend on f_wat which can
   % be lower than the maximum f_wat due to sublimation. Although I've already
   % learned this lesson when I added them to the beginning of ICEENBAL, it
   % slipped my mind when thinking through all the mass balance assert
   % statements. However, it should't matter because fliqmin is not relevant to
   % those checks, and fliqmax can only get smaller than f_wat_max which is the
   % ro_ice / ro_liq condition, and the checks are all designed to ensure we
   % don't exceeed that. BUT, for numerical error e.g. in MZTRANSFORM, it's
   % necessary to ensure the solution does not skip the melt zone in terms of
   % fliqmin/max rather than ro_ice / ro_liq etc.

   % These must be true:
   %
   % f_ice + f_liq * ro_liq / ro_ice <= 1
   %
   % f_ice * ro_ice / ro_liq + f_liq <= ro_ice / ro_liq
   %
   % f_wat <= ro_ice / ro_liq
   %
   % Thus, if adding new water:
   %
   % f_ice + (f_liq + d_liq) * ro_liq / ro_ice <= 1
   %
   % f_ice * ro_ice / ro_liq + f_liq + d_liq <= ro_ice / ro_liq
   %
   % d_liq <= ro_ice / ro_liq * (1 - f_ice(old)) - f_liq(old)
   %
   % Or:
   %
   % d_liq_max = ro_ice / ro_liq * (1 - f_ice(old)) - f_liq(old)
   %
   % f_ice + f_liq + d_liq <= f_wat
   %
   % If all ice melts:
   %
   % f_wat = f_liq(new)
   %       = f_liq(old) + df_liq
   %
   % Thus the requirement is:
   %
   % f_liq(old) + df_liq <= ro_ice / ro_liq
   % df_liq <= ro_ice / ro_liq - f_liq(old)
   %

   %
   % But the implementation should probably use a small f_ice threshold to
   % disallow complete melting, and/or call combine layers.


   maxfill = ro_ice / ro_liq * (1 - f_ice) - f_liq;

   % therefore:
   f_liq = max()


   assertF(@() all(f_ice + f_liq * ro_wie <= 1 + eps))
end
