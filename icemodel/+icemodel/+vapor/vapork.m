function k_vap = vapork(T, f_liq, varargin)
   %vapork Compute the vapor thermal diffusion coefficient.
   %
   %  k_vap = icemodel.vapor.vapork(T, f_liq)
   %  k_vap = icemodel.vapor.vapork(T, f_liq, dro_vapdT)
   %
   %  Computes the effective thermal conductivity contribution from vapor
   %  diffusion through porous ice, following Anderson (1976):
   %
   %     k_vap = L * De * dro_vap/dT   [W m-1 K-1]
   %
   %  where De is the effective vapor diffusivity (from icemodel.vapor.vapordiffusivity),
   %  dro_vap/dT is the temperature derivative of saturation vapor density
   %  (from icemodel.vapor.vapordensity), and L is the latent heat.
   %
   %  When dro_vapdT is supplied, icemodel.vapor.vapork reuses it directly so callers such
   %  as ICEENBAL can avoid a second vapor-density derivative evaluation.
   %
   %  Phase awareness: uses Ls (sublimation) for dry cells and Lv
   %  (vaporization) for wet cells
   %  (f_liq > f_liq_phase_switch_threshold), matching the phase of the
   %  saturation vapor pressure computation in icemodel.vapor.vapordensity.
   %
   % See also: icemodel.vapor.vapordensity, icemodel.vapor.vapordiffusivity,
   %  BULKTHERMALK
   %
   %#codegen

   persistent Lv Ls f_liq_phase_switch_threshold
   if isempty(Lv)
      [Lv, Ls] = icemodel.physicalConstant('Lv', 'Ls');
      f_liq_phase_switch_threshold = icemodel.parameterLookup( ...
         'f_liq_phase_switch_threshold');
   end

   if nargin > 3
      error('vapork:UnrecognizedNumberOfInputs', 'Unrecognized number of inputs.')
   end

   % Vapor density derivative [kg m-3 K-1]
   if nargin == 3
      dro_vapdT = varargin{1};
   else
      [~, dro_vapdT] = icemodel.vapor.vapordensity(T, f_liq);
   end

   % Vapor diffusivity [m2 s-1]
   De = icemodel.vapor.vapordiffusivity(T);

   % Vapor thermal diffusion coefficient [W m-1 K-1]
   k_vap = Ls * De .* dro_vapdT;

   % Switch to Lv for wet cells
   wet = f_liq > f_liq_phase_switch_threshold;
   if any(wet)
      k_vap(wet) = Lv * De(wet) .* dro_vapdT(wet);
   end
end
