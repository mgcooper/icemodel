function f_liq_min = residual_water_fraction(f_ice, f_liq)
   %RESIDUAL_WATER_FRACTION Minimum residual liquid fraction (volumetric).
   %
   %  f_liq_min = icemodel.column.residual_water_fraction(f_ice, f_liq)
   %
   % Returns the minimum residual unfrozen liquid fraction in volumetric units
   % [m³_liq m³_CV⁻¹] — the same reference frame as f_liq — for each control
   % volume, defined by the Jordan (1991) phase-fraction characteristic curve
   % evaluated at the lower melt-zone temperature boundary TL:
   %
   %   f_wat     = f_liq + f_ice * ro_ice / ro_liq   (total water fraction)
   %   f_liq_min = f_wat * f_ell_min                 (meltzone_bounds)
   %
   % This is the volumetric companion to
   % icemodel.column.residual_water_pore_fraction, which returns the same
   % quantity normalized to pore volume. The two are related by:
   %
   %   f_liq_min_vol  = f_liq_min                           (this function)
   %   f_liq_min_pore = f_liq_min / (1 - f_ice)             (pore-fraction)
   %
   % Use the pore-fraction form for comparisons against opts.f_liq_resid (which
   % is pore-volume referenced following Jordan 1991).  Use this volumetric form
   % when the result needs to be compared or combined directly with f_liq (e.g.,
   % in infiltration models or as a floor condition on f_liq itself).
   %
   % Inputs
   %   f_ice  - Volumetric ice fraction [-], scalar or any array.
   %   f_liq  - Volumetric liquid-water fraction [-], same shape as f_ice.
   %
   % Output
   %   f_liq_min - Minimum volumetric liquid fraction [-], same shape as inputs.
   %
   % See also: icemodel.column.residual_water_pore_fraction,
   %           icemodel.column.water_fraction,
   %           icemodel.column.meltzone_bounds
   %
   %#codegen

   f_wat     = icemodel.column.water_fraction(f_ice, f_liq);
   f_liq_min = icemodel.column.meltzone_bounds(f_wat);
end
