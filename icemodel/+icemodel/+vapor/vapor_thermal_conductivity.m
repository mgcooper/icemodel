function [k_vap, De, dk_vapdT] = vapor_thermal_conductivity(T, f_liq, varargin)
   %VAPOR_THERMAL_CONDUCTIVITY Effective thermal conductivity from vapor
   % diffusion [W m-1 K-1].
   %
   %  k_vap = icemodel.vapor.vapor_thermal_conductivity(T, f_liq)
   %  k_vap = ...
   %     icemodel.vapor.vapor_thermal_conductivity(T, f_liq, dro_vapdT)
   %  [k_vap, De] = ...
   %     icemodel.vapor.vapor_thermal_conductivity(T, f_liq, dro_vapdT)
   %  [k_vap, De, dk_vapdT] = ...
   %     icemodel.vapor.vapor_thermal_conductivity(T, f_liq)
   % 
   %  Computes the effective thermal conductivity contribution from vapor
   %  diffusion through porous ice, following Anderson (1976):
   %
   %     k_vap = L * De * dro_vap/dT   [W m-1 K-1]
   %
   %  where De is the effective vapor diffusivity (from
   %  icemodel.vapor.vapor_diffusivity), dro_vap/dT is the temperature
   %  derivative of saturation vapor density (from
   %  icemodel.vapor.saturation_vapor_density), and L is the latent heat.
   %
   %  When you supply dro_vapdT, this function reuses it directly, so callers
   %  such as `icemodel.column.solve_column_enthalpy` avoid a second
   %  vapor-density derivative evaluation.
   %
   %  The second output returns the effective vapor diffusivity this function
   %  evaluates. The coupled vapor model needs the same De to build its face
   %  quantities. Asking for it here keeps the (T/Tf)^nd power off the hot path
   %  a second time.
   %
   %  The third output returns dk_vap/dT [W m-1 K-2]:
   %
   %     dk_vap/dT = L * (dDe/dT * dro_vap/dT + De * d2ro_vap/dT2)
   %
   %  with dDe/dT = nd/T * De. A Newton linearization of the enthalpy equation
   %  needs this term, whereas the Picard form in
   %  icemodel.column.assemble_enthalpy_system does not. The latent heat
   %  switches on f_liq, not T, so no L derivative is included here.
   %
   %  The third output needs d2ro_vap/dT2, which no caller supplies, so asking
   %  for it evaluates both density derivatives from T and f_liq and ignores a
   %  the optional dro_vapdT argument even if it's supplied. That keeps the two
   %  derivatives consistent. icemodel.vapor.saturation_vapor_density documents
   %  both in the Ambaum (2020) and Romps (2021) coefficient form that the
   %  production model uses.
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

   persistent nd
   if isempty(nd)
      nd = icemodel.parameterLookup('nd');
   end

   % Vapor density derivatives [kg m-3 K-1] and [kg m-3 K-2]. The second
   % derivative comes from the same call, so the pair stays consistent.
   if nargout > 2
      [~, dro_vapdT, d2ro_vapdT2] = ...
         icemodel.vapor.saturation_vapor_density(T, f_liq);
   elseif nargin == 3
      dro_vapdT = varargin{1};
   else
      [~, dro_vapdT] = icemodel.vapor.saturation_vapor_density(T, f_liq);
   end

   % Vapor diffusivity [m2 s-1]
   De = icemodel.vapor.vapor_diffusivity(T);

   % Phase-aware latent heat: Ls for dry/cold cells, Lv for wet cells.
   Lv = icemodel.vapor.latent_enthalpy_switch(f_liq);

   % Vapor thermal conductivity [W m-1 K-1]
   k_vap = Lv .* De .* dro_vapdT;

   % Conductivity derivative [W m-1 K-2]. The product form avoids dividing by
   % dro_vapdT, which becomes small at low temperature.
   if nargout > 2
      dk_vapdT = Lv .* De .* (nd ./ T .* dro_vapdT + d2ro_vapdT2);
   end
end
