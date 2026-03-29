function [ro_vap, dro_vapdT, d2ro_vapdT2] = VAPORDENSITY(T, f_liq)
   %VAPORDENSITY Compute saturation vapor density in porous ice.
   %
   %  ro_vap = VAPORDENSITY(T, f_liq)
   %  [ro_vap, dro_vapdT, d2ro_vapdT2] = VAPORDENSITY(T, f_liq)
   %
   %  Computes the equilibrium (saturation) water vapor density and its
   %  temperature derivative within the air voids of porous ice. The air
   %  voids are assumed saturated with respect to either ice or liquid water
   %  depending on the local liquid fraction.
   %
   %  Phase awareness: saturation vapor pressure is computed over ice by
   %  default, switching to over liquid where
   %  f_liq > f_liq_phase_switch_threshold. This is consistent with the
   %  assumption that wet cells have pore air saturated with respect to
   %  liquid water.
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
   % See also: VAPORK, VAPPRESS, VAPORTRANSFER
   %
   %#codegen

   persistent Rv f_liq_phase_switch_threshold
   if isempty(Rv)
      Rv = icemodel.physicalConstant('Rv');
      f_liq_phase_switch_threshold = icemodel.parameterLookup( ...
         'f_liq_phase_switch_threshold');
   end

   if nargout > 2
      [es, des_dT, d2es_dT2] = VAPPRESS(T, false);
   elseif nargout > 1
      [es, des_dT] = VAPPRESS(T, false);
   else
      es = VAPPRESS(T, false);
   end

   % Override with liquid for wet cells
   wet = f_liq > f_liq_phase_switch_threshold;
   if any(wet)
      if nargout > 2
         [es(wet), des_dT(wet), d2es_dT2(wet)] = VAPPRESS(T(wet), true);
      elseif nargout > 1
         [es(wet), des_dT(wet)] = VAPPRESS(T(wet), true);
      else
         es(wet) = VAPPRESS(T(wet), true);
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
