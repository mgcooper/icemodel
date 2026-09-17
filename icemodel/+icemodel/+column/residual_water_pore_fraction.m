function f_res_por = residual_water_pore_fraction(f_ice, f_liq)
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
   %   f_res_por = f_liq_min / (1 - f_ice)           (convert to per-pore-volume)
   %
   % The result uses the same pore-volume reference frame as the model
   % parameters opts.f_res_pore_ice, opts.f_res_pore_snow and
   % opts.f_res_pore_firn, so you can use it directly in max(f_res_por, ...).
   % Multiply by (1 - f_ice) to get the volumetric liquid fraction; see
   % icemodel.column.residual_water_fraction for that form.
   %
   % The residual floor keeps evaporation from reducing f_liq below the
   % minimum physically consistent value for a melting node (T_ice > TL).
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

   f_wat = icemodel.column.water_fraction(f_ice, f_liq);
   f_liq_min = icemodel.column.meltzone_bounds(f_wat);
   f_res_por = f_liq_min ./ (1 - f_ice);
end
