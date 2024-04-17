function [k_eff, ro_sno, cp_sno, dFdT, drovdT, H, ro_vap] = UPDATESTATE( ...
      T, f_ice, f_liq, ro_ice, ro_liq, ro_air, cv_ice, cv_liq, k_liq, roLf, ...
      Ls, Rv, Tf, fcp)
   %UPDATESTATE Update state variables
   %
   % Note: These are "intrinsic" volumetric heat capacities and must be
   % multiplied by the constituent volumetric fractions to get the "bulk"
   % volumetric heat capacities [J/m3/K]
   % cv_ice = cp_ice * ro_ice
   % cv_liq = cp_liq * ro_liq

   % Bulk density (pt) (eq 3) [kg m-3]
   ro_sno = ro_ice * f_ice + ro_liq * f_liq + ro_air * (1.0 - f_liq - f_ice);

   % Bulk mass-specific heat capacity (ct) (eq 55) [J kg-1 K-1]
   cp_sno = (cv_ice * f_ice + cv_liq * f_liq) ./ ro_sno;

   % Change in saturation vapor density with temperature [kg m-3]
   [ro_vap, drovdT, k_vap] = VAPORHEAT(T, f_liq, f_ice, Tf, Rv, Ls);

   % Total enthalpy [J m-3]
   H = TOTALHEAT(T, f_ice, f_liq, cv_ice, cv_liq, roLf, Ls * ro_vap, Tf);

   % Effective thermal conductivity [W m-1 K-1]
   k_eff = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, k_vap);

   % Derivative of liquid water fraction wrt temperature (freezing curve) [K-1]
   dFdT = 2.0 * fcp ^ 2.0 * (Tf - min(T, Tf)) .* f_wat ...
      ./ (1.0 + fcp ^ 2.0 * (Tf - min(T, Tf)) .^ 2.0) .^ 2.0;

   % Energy coefficient [W m-2 K-1]
   % aP0 = (ro_sno .* cp_sno + Lf * ro_liq * dFdT ...
   %    + Ls * (1.0 - f_ice - f_liq) .* drovdT) .* dz / dt;

   % Update the general equation coefficients
   % [aN, aP, aS, b, iM, a1, a2] = GECOEFS(T, ro_sno, cp_sno, f_liq, f_ice, Ls, ...
   %    Lf, ro_liq, dz, dt, dFdT, drovdT, TL, H, H_old, Sc, k_eff, fn, ...
   %    delz, Ts, JJ, Fc, Fp, bc_type);
end
