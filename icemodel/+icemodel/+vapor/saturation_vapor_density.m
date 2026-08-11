function [ro_vap, dro_vapdT, d2ro_vapdT2] = saturation_vapor_density(T, f_liq)
   %saturation_vapor_density Compute saturation vapor density in porous ice.
   %
   %  ro_vap = icemodel.vapor.saturation_vapor_density(T, f_liq)
   %  [ro_vap, dro_vapdT, d2ro_vapdT2] = icemodel.vapor.saturation_vapor_density(T, f_liq)
   %
   %  Computes the equilibrium (saturation) water vapor density and its
   %  temperature derivative within the air voids of porous ice. This
   %  function assumes the air voids are saturated over ice or over liquid
   %  water, based on the local liquid fraction.
   %
   %  Phase awareness: this function computes the saturation vapor pressure
   %  over ice by default. It switches to over liquid where
   %  f_liq > f_liq_phase_switch_threshold. That matches the assumption that
   %  the pore air of a wet cell is saturated over liquid water.
   %
   %  ro_vap is vapor mass per vapor volume, analogous to an intrinsic density.
   %
   %  Outputs:
   %     ro_vap    - Saturation vapor density [kg m-3]
   %     dro_vapdT - Derivative of vapor density wrt temperature [kg m-3 K-1]
   %     d2ro_vapdT2 - Second derivative wrt temperature [kg m-3 K-2]
   %
   %  The vapor-density derivatives can be written either in terms of the
   %  pressure derivatives:
   %
   %     dro_vap/dT = (des_dT - es/T) / (Rv*T)
   %     d2ro_vap/dT2 = (d2es_dT2 - 2*des_dT/T + 2*es/T^2) / (Rv*T)
   %
   %  or directly in coefficient form:
   %
   %     dro_vap/dT = ro_vap / T * (c - b/T - 1)
   %     d2ro_vap/dT2 = ro_vap / T^2 * ((c-2) * (c - 2*b/T - 1) + b^2/T^2)
   %
   % See also: icemodel.vapor.vapor_thermal_diffusion_coefficient,
   %  icemodel.vapor.saturation_vapor_pressure,
   %  icemodel.column.vapor_mass_transfer
   %
   %#codegen

   persistent Rv f_liq_phase_switch_threshold
   if isempty(Rv)
      Rv = icemodel.physicalConstant('Rv');
      f_liq_phase_switch_threshold = icemodel.parameterLookup( ...
         'f_liq_phase_switch_threshold');
   end

   if nargout > 2
      [es, des_dT, d2es_dT2] = icemodel.vapor.saturation_vapor_pressure(T, false);
   elseif nargout > 1
      [es, des_dT] = icemodel.vapor.saturation_vapor_pressure(T, false);
   else
      es = icemodel.vapor.saturation_vapor_pressure(T, false);
   end

   % Override with liquid for wet cells
   wet = f_liq > f_liq_phase_switch_threshold;
   if any(wet)
      if nargout > 2
         [es(wet), des_dT(wet), d2es_dT2(wet)] = ...
            icemodel.vapor.saturation_vapor_pressure(T(wet), true);
      elseif nargout > 1
         [es(wet), des_dT(wet)] = ...
            icemodel.vapor.saturation_vapor_pressure(T(wet), true);
      else
         es(wet) = icemodel.vapor.saturation_vapor_pressure(T(wet), true);
      end
   end

   % Equilibrium water vapor density [kg m-3]
   ro_vap = es ./ (Rv * T);

   % Derivatives
   if nargout > 1
      dro_vapdT = (des_dT - es ./ T) ./ (Rv * T); % [kg m-3 K-1]
   end

   if nargout > 2
      d2ro_vapdT2 = (d2es_dT2 - 2 * des_dT ./ T + 2 * es ./ T .^ 2) ...
         ./ (Rv * T);
   end
end
