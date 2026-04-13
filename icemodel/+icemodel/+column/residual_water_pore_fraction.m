function f_liq_res = residual_water_pore_fraction(f_ice, f_liq)
   %RESIDUAL_WATER_PORE_FRACTION Residual liquid fraction per pore volume.
   %
   %  f_liq_res = icemodel.column.residual_water_pore_fraction(f_ice, f_liq)
   %
   % Returns the minimum residual unfrozen liquid fraction per pore volume
   % [m³_liq m³_pore⁻¹] for each control volume, defined by the Jordan (1991)
   % phase-fraction characteristic curve evaluated at the lower melt-zone
   % temperature boundary TL:
   %
   %   f_wat     = f_liq + f_ice * ro_ice / ro_liq   (total water fraction)
   %   f_liq_min = f_wat * f_ell_min                 (minimum vol. liquid, meltzone_bounds)
   %   f_liq_res = f_liq_min / (1 - f_ice)           (convert to per-pore-volume)
   %
   % The result is in the same reference frame as the model's opts.f_liq_resid
   % parameter (pore-volume), so it can be used directly in max(f_liq_res, ...).
   % Multiply by (1 - f_ice) to obtain the volumetric liquid fraction; see
   % icemodel.column.residual_water_fraction for that form.
   %
   % The residual floor ensures that evaporation does not reduce f_liq below
   % the minimum physically consistent value for a melting node (T > TL).
   %
   % Inputs
   %   f_ice  - Volumetric ice fraction [-], scalar or column vector.
   %   f_liq  - Volumetric liquid-water fraction [-], same shape as f_ice.
   %
   % Output
   %   f_liq_res - Residual liquid fraction per pore volume [-].
   %
   % See also: icemodel.column.residual_water_fraction,
   %           icemodel.column.water_fraction,
   %           icemodel.column.meltzone_bounds,
   %           icemodel.column.budget_surface_mass_balance
   %
   %#codegen

   f_wat     = icemodel.column.water_fraction(f_ice, f_liq);
   f_liq_min = icemodel.column.meltzone_bounds(f_wat);
   f_liq_res = f_liq_min ./ (1 - f_ice);
end
