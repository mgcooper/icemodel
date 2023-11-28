function [k_eff, ro_sno, cp_sno, drovdT] = UPDATESTATE(T, f_ice, f_liq, ...
      ro_ice, ro_liq, ro_air, cp_ice, cp_liq, k_liq, Ls, Rv, Tf)
   %UPDATESTATE Update state variables

   % Consituent bulk densities kg/m3
   bd_ice = ro_ice * f_ice;
   bd_liq = ro_liq * f_liq;
   bd_air = ro_air * (1.0 - f_liq - f_ice);

   % Combined bulk density (pt) kg/m3
   ro_sno = (bd_ice + bd_liq + bd_air); % (eq 3)

   % Constituent bulk volume-specific heat capacities J/m3/K
   cv_ice = cp_ice * bd_ice;
   cv_liq = cp_liq * bd_liq;
   % cv_air = cp_air * bd_air;

   % Combined bulk mass-specific heat capacity (ct) J/kg/K
   cp_sno = (cv_ice + cv_liq) ./ ro_sno; % (eq 55)

   % Change in saturation vapor density with temperature kg/m3
   [~, drovdT, k_vap] = VAPORHEAT(T, f_liq, f_ice, Tf, Rv, Ls);
   
   % Effective thermal conductivity
   k_eff = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, k_vap);

   % Energy coefficient W/m2/K
   % aP0 = bd_sno .* cp_sno .* dz ./ dt;
end
