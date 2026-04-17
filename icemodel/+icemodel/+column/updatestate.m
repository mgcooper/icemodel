function [H, k_eff, dHdT, dFdT, drovdT, ro_vap, ro_sno, cp_sno] = ...
      updatestate(T, f_ice, f_liq, f_wat)
   %UPDATESTATE Update column thermodynamic state variables.
   %
   % Note: if this is called within solver iterations, ensure f_wat = f_wat_old
   %
   %#codegen

   % Saturation vapor density and temperature derivative [kg m-3, kg m-3 K-1]
   [ro_vap, drovdT] = icemodel.vapor.saturation_vapor_density( ...
      T, f_liq);

   % Vapor thermal diffusion coefficient [W m-1 K-1]
   k_vap = icemodel.vapor.vapor_thermal_diffusion_coefficient( ...
      T, f_liq, drovdT);

   % Effective thermal conductivity [W m-1 K-1]
   k_eff = icemodel.column.bulk_thermal_conductivity( ...
      T, f_ice, f_liq, k_vap);

   % Total enthalpy and temperature derivative [J m-3] and [J m-3 K-1]
   [H, dHdT, dFdT] = icemodel.column.bulk_enthalpy( ...
      T, f_ice, f_liq, f_wat, ro_vap);

   % Optional bulk density (pt) (eq 3) [kg m-3]
   if nargout > 6
      ro_sno = icemodel.column.bulk_density( ...
         f_ice, f_liq);
   end

   % Optional bulk mass-specific heat capacity (ct) (eq 55) [J kg-1 K-1]
   if nargout > 7
      cp_sno = icemodel.column.bulk_specific_heat_capacity( ...
         f_ice, f_liq, ro_sno);
   end

end
