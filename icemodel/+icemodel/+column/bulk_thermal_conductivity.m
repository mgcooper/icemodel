function [k_eff, k_vap] = bulk_thermal_conductivity(T, f_ice, f_liq, varargin)
   %BULK_THERMAL_CONDUCTIVITY Compute bulk effective thermal conductivity.
   %
   %  [k_eff, k_vap] = bulk_thermal_conductivity(T, f_ice, f_liq)
   %  [k_eff, k_vap] = bulk_thermal_conductivity(T, f_ice, f_liq, k_vap)
   %
   %  Combines three transport mechanisms into a volume-weighted effective
   %  thermal conductivity (denoted gamma or ke in Patankar Eq. 4.9):
   %
   %     k_eff = f_liq * k_liq + f_ice * k_ice + f_air * k_vap
   %
   %  where f_air = 1 - f_liq - f_ice.
   %
   %  The solid-phase conductivity k_ice comes from
   %  firn_thermal_conductivity (Calonne 2019 density- and temperature-
   %  dependent parameterization). See firn_thermal_conductivity for the
   %  full Calonne (2019) Eq. 5 formulation.
   %
   %  The vapor component k_vap supports two calling conventions:
   %     nargin=3: k_vap computed internally via
   %           icemodel.vapor.vapor_thermal_diffusion_coefficient(T, f_liq)
   %     nargin=4: k_vap provided externally (including explicit 0)
   %
   % See also: icemodel.column.firn_thermal_conductivity,
   %  icemodel.vapor.vapor_thermal_diffusion_coefficient,
   %  icemodel.vapor.saturation_vapor_density, icemodel.vapor.vapor_diffusivity
   %
   %#codegen

   persistent k_liq
   if isempty(k_liq)
      k_liq = icemodel.physicalConstant('k_liq');
   end

   % Compute dry snow/firn/ice thermal conductivity (Calonne 2019 Eq. 5)
   k_ice = icemodel.column.firn_thermal_conductivity(T, f_ice);

   % Compute vapor thermal diffusion coefficient
   if nargin < 4
      k_vap = icemodel.vapor.vapor_thermal_diffusion_coefficient(T, f_liq);
   else
      % k_vap provided by an external model.
      k_vap = varargin{1};
   end

   % Volume-weighted bulk effective conductivity:
   k_eff = f_liq .* k_liq + f_ice .* k_ice + (1.0 - f_liq - f_ice) .* k_vap;
end
