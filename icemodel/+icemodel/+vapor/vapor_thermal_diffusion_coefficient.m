function k_vap = vapor_thermal_diffusion_coefficient(T, f_liq, varargin)
   %vapor_thermal_diffusion_coefficient Vapor thermal diffusion coefficient.
   %
   %  k_vap = icemodel.vapor.vapor_thermal_diffusion_coefficient(T, f_liq)
   %  k_vap = ...
   %     icemodel.vapor.vapor_thermal_diffusion_coefficient(T, f_liq, dro_vapdT)
   %
   %  Computes the effective thermal conductivity contribution from vapor
   %  diffusion through porous ice, following Anderson (1976):
   %
   %     k_vap = L * De * dro_vap/dT   [W m-1 K-1]
   %
   %  where De is the effective vapor diffusivity
   %  (from icemodel.vapor.vapor_diffusivity), dro_vap/dT is the temperature
   %  derivative of saturation vapor density
   %  (from icemodel.vapor.saturation_vapor_density), and L is the latent heat.
   %
   %  When dro_vapdT is supplied,
   %  icemodel.vapor.vapor_thermal_diffusion_coefficient reuses it directly so
   %  callers such as `icemodel.column.solve_column_enthalpy` can avoid a second
   %  vapor-density derivative evaluation.
   %
   %  Phase awareness: uses Ls (sublimation) for dry cells and Lv (vaporization)
   %  for wet cells (f_liq > f_liq_phase_switch_threshold), matching the phase
   %  of the saturation vapor pressure computation in
   %  icemodel.vapor.saturation_vapor_density.
   %
   % See also: icemodel.vapor.saturation_vapor_density,
   %           icemodel.vapor.vapor_diffusivity,
   %           icemodel.column.bulk_thermal_conductivity
   %
   %#codegen

   % Vapor density derivative [kg m-3 K-1]
   if nargin == 3
      dro_vapdT = varargin{1};
   else
      [~, dro_vapdT] = icemodel.vapor.saturation_vapor_density(T, f_liq);
   end

   % Vapor diffusivity [m2 s-1]
   De = icemodel.vapor.vapor_diffusivity(T);

   % Phase-aware latent heat: Ls for dry/cold cells, Lv for wet cells.
   Lv = icemodel.vapor.latent_enthalpy_switch(f_liq);

   % Vapor thermal diffusion coefficient [W m-1 K-1]
   k_vap = Lv .* De .* dro_vapdT;
end
