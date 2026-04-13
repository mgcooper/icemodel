function f_liq_res = residual_water_fraction(f_ice, f_liq)
   %RESIDUAL_WATER_FRACTION Residual liquid fraction per pore volume.
   %
   %  f_liq_res = icemodel.column.residual_water_fraction(f_ice, f_liq)
   %
   % Returns the residual unfrozen liquid fraction per pore volume for each
   % control volume.  The result is defined by the phase-fraction
   % characteristic function evaluated at the lower melt-zone temperature
   % boundary TL:
   %
   %   f_wat     = f_liq + f_ice * ro_ice / ro_liq
   %   f_liq_min = f_wat * f_ell_min               (meltzone_bounds)
   %   f_liq_res = f_liq_min / (1 - f_ice)         (per pore volume)
   %
   % The function works on scalar inputs or on full column vectors. When
   % called from budget_surface_mass_balance only the top node is passed
   % (f_ice(1), f_liq(1)) so that the returned scalar can be used directly
   % in the evaporation floor logic. When called diagnostically the full
   % column vectors may be passed instead.
   %
   % The residual floor ensures that evaporation does not reduce f_liq below
   % the minimum physically consistent value for a melting node (T > TL).
   %
   % See also: icemodel.column.water_fraction,
   %           icemodel.column.meltzone_bounds,
   %           icemodel.column.budget_surface_mass_balance
   %
   %#codegen

   f_wat     = icemodel.column.water_fraction(f_ice, f_liq);
   f_liq_min = icemodel.column.meltzone_bounds(f_wat);
   f_liq_res = f_liq_min ./ (1 - f_ice);
end
