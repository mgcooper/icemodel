function [H, dHdT, dLdT] = total_enthalpy(T, f_ice, f_liq, f_wat, ro_vap)
   %TOTAL_ENTHALPY Compute total enthalpy (J/m3).
   %
   %#codegen

   persistent Tf cv_ice cv_liq roLf Ls Lv f_liq_phase_switch_threshold
   if isempty(Tf)
      [Tf, cv_ice, cv_liq, roLf, Ls, Lv] = icemodel.physicalConstant( ...
         'Tf', 'cv_ice', 'cv_liq', 'roLf', 'Ls', 'Lv');
      f_liq_phase_switch_threshold = icemodel.parameterLookup( ...
         'f_liq_phase_switch_threshold');
   end

   % Option to ignore vapor heat.
   if nargin < 5
      ro_vap = 0;
   end

   % Vapor heat wrt sublimation [J m-3]
   roLv = ro_vap * Ls;

   % Vapor heat wrt evaporation (switch to Lv for wet cells)
   wet = f_liq > f_liq_phase_switch_threshold;
   if any(wet)
      roLv(wet) = ro_vap(wet) * Lv;
   end

   % Total (bulk) enthalpy [J m-3]
   H = ( ...
      (cv_ice * f_ice + cv_liq * f_liq) .* (T - Tf) ...  % specific heat
      + roLf .* f_liq ...                                % latent heat
      + roLv .* (1.0 - f_liq - f_ice) ...                % vapor heat
      );

   % Temperature derivative (specific heat)
   if nargout > 1
      dHdT = cv_ice * f_ice + cv_liq * f_liq; % = cp_bulk * ro_bulk
   end

   % Temperature derivative (latent heat).
   if nargout > 2
      dLdT = icemodel.column.liquid_fraction_derivative(T, f_ice, f_liq, f_wat);
   end
end
