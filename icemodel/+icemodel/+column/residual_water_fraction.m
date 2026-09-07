function [f_res, f_liq_min] = residual_water_fraction(f_ice, f_liq, f_res_por)
   %RESIDUAL_WATER_FRACTION Volumetric residual-water floor for control volumes.
   %
   %  [f_res, f_liq_min] = icemodel.column.residual_water_fraction( ...
   %     f_ice, f_liq, f_res_por)
   %
   % Returns the volumetric residual-water floor f_res (the minimum liquid
   % fraction that must be retained) as the maximum of two contributions:
   %
   %  1. Capillary: f_res_por .* (1 - f_ice)
   %        Capillary/hydraulic irreducible saturation scaled to CV volume.
   %        f_res_por is the per-pore-volume parameter (opts.f_res_pore_ice,
   %        opts.f_res_pore_snow, or opts.f_res_pore_firn). This contribution
   %        dominates for ice and snow because their mushy-zone thermodynamic
   %        minimum liquid water content is extremely small (near zero).
   %
   %  2. Thermodynamic: Thermodynamic minimum from Jordan's (1991)
   %        phase-fraction curve evaluated at the lower melt-zone boundary TL:
   %        f_liq_min = meltzone_bounds(water_fraction(f_ice, f_liq)).
   %        Near-zero for Jp=0 (ice/snow) but non-negligible for soils or
   %        well-saturated firn.
   %
   %  This function returns the maximum of these two minimums:
   %   f_res = max(f_res_por .* (1 - f_ice), f_liq_min)
   %
   % Inputs
   %   f_ice     - Volumetric ice fraction [-], scalar or any array.
   %   f_liq     - Volumetric liquid-water fraction [-], same shape as f_ice.
   %   f_res_por - Capillary residual per pore volume [-], scalar or same shape.
   %               Pass opts.f_res_pore_ice, opts.f_res_pore_snow, or
   %               opts.f_res_pore_firn depending on the surface type.
   %
   % Outputs
   %   f_res     - Volumetric residual-water floor [-], same shape as f_ice.
   %   f_liq_min - Phase fraction curve thermodynamic minimum (volumetric) [-].
   %
   % See also: icemodel.column.residual_water_pore_fraction,
   %           icemodel.column.water_fraction,
   %           icemodel.column.meltzone_bounds
   %
   %#codegen

   f_wat     = icemodel.column.water_fraction(f_ice, f_liq);
   f_liq_min = icemodel.column.meltzone_bounds(f_wat);
   f_res     = max(f_res_por .* (1 - f_ice), f_liq_min);
end
