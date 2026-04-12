function [TL, TH, f_ell_min, f_ell_max] = meltzone_bounds()
   %MELTZONE_BOUNDS Return the canonical mushy-zone temperature bounds.
   %
   %  [TL, TH, F_ELL_MIN, F_ELL_MAX] = icemodel.column.meltzone_bounds()
   %
   %  TL and TH are the lower and upper melt-zone temperature bounds [K].
   %  F_ELL_MIN and F_ELL_MAX are the corresponding liquid mass fractions on
   %  the Jordan phase-fraction curve at TL and TH.
   %
   %#codegen

   persistent Tf Lf cp_ice cp_liq fcp TL_cached TH_cached f_ell_min_cached ...
      f_ell_max_cached
   if isempty(TL_cached)
      [Tf, Lf, cp_ice, cp_liq] = icemodel.physicalConstant( ...
         'Tf', 'Lf', 'cp_ice', 'cp_liq');
      fcp = icemodel.parameterLookup('fcp');

      TL_cached = Tf - (2.0 * Lf / (fcp ^ 2.0 * cp_ice)) ^ (1.0 / 3.0);
      TH_cached = Tf - cp_liq / (Lf * 2.0 * fcp ^ 2.0);
      f_ell_min_cached = 1.0 / (1.0 + (fcp * (Tf - TL_cached)) ^ 2.0);
      f_ell_max_cached = 1.0 / (1.0 + (fcp * (Tf - TH_cached)) ^ 2.0);
   end

   TL = TL_cached;
   TH = TH_cached;
   f_ell_min = f_ell_min_cached;
   f_ell_max = f_ell_max_cached;
end
