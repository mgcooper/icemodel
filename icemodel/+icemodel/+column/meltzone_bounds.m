function [f_liq_min, f_liq_max] = meltzone_bounds(f_wat)
   %MELTZONE_BOUNDS Return the canonical mushy-zone liquid fraction bounds.
   %
   %  [F_LIQ_MIN, F_LIQ_MAX] = icemodel.column.meltzone_bounds(F_WAT)
   %
   %  TL and TH are the lower and upper melt-zone temperature bounds [K].
   %  F_ELL_MIN and F_ELL_MAX are the corresponding liquid mass fractions on
   %  the Jordan phase-fraction curve at TL and TH given the current F_WAT.
   %
   %  F_LIQ_MIN and F_LIQ_MAX are the corresponding liquid volume fractions.
   %
   %#codegen

   persistent f_ell_min f_ell_max
   if isempty(f_ell_min)
      [f_ell_min, f_ell_max] = icemodel.parameterLookup( ...
         'f_ell_min', 'f_ell_max');
   end

   f_liq_min = f_wat .* f_ell_min;
   f_liq_max = f_wat .* f_ell_max;
end
