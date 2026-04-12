function [H, dHdT, dLdT, dVdT, k_vap] = total_enthalpy_2(T, f_ice, f_liq)
   %TOTAL_ENTHALPY_2 Compute total enthalpy and its temperature derivative.
   %
   %#codegen

   % Experimental helper. Keep the calculation shape aligned with
   % solve_column_enthalpy.
   persistent cv_ice cv_liq ro_liq Lf Ls Tf
   if isempty(cv_ice)
      [cv_ice, cv_liq, ro_liq, Lf, Ls, Tf] = icemodel.physicalConstant( ...
         'cv_ice', 'cv_liq', 'ro_liq', 'Lf', 'Ls', 'Tf');
   end
   f_wat = icemodel.water_fraction(f_ice, f_liq);

   % Compute the derivative of enthalpy wrt temperature
   dHdT = cv_ice * f_ice + cv_liq * f_liq;
   dLdT = icemodel.column.liquid_fraction_derivative(T, [], [], f_wat);

   % Compute vapor density, its derivative wrt temperature, and k_vap
   [ro_vap, dVdT] = icemodel.vapor.saturation_vapor_density(T, f_liq);
   k_vap = icemodel.vapor.vapor_thermal_diffusion_coefficient(T, f_liq, dVdT);

   % Compute total enthalpy wrt the reference temperature.
   H = icemodel.column.total_enthalpy(T, f_ice, f_liq, f_wat, ro_vap);
end
