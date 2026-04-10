function [T_sfc, ok] = solve_surface_temperature(T_sfc, tair, Qsi, Qli, albedo, ...
      wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, roL, liqflag, chi, k_eff, ...
      T_ice, dz)
   %SOLVE_SURFACE_TEMPERATURE Solve the explicit bulk-Richardson SEB for Ts.
   %
   %  [T_sfc, ok] = icemodel.surface.solve_surface_temperature(T_sfc, tair, ...)
   %
   % Solves the nonlinear surface energy balance for T_sfc using a
   % Newton-Raphson iteration with an analytical Jacobian. The residual
   % being minimized is:
   %
   %   F(T_sfc) = chi*Qsi*(1-albedo) + emiss*(Qli - SB*T_sfc^4)
   %              + Qh(T_sfc) + Qe(T_sfc) + Qc(T_sfc) + Qa
   %
   % where Qc(T_sfc) = k_eff(1) * (T_ice(1) - T_sfc) / (dz(1)/2) couples
   % the surface temperature to the top ice layer, and Qh, Qe are evaluated
   % using the bulk-Richardson turbulence closure. The analytical Jacobian is:
   %
   %   dF/dT_sfc = -4 * emiss * SB * T_sfc^3
   %               + d(stability)/dT * (Qh/stability + Qe/stability)
   %               + stability * (dQh/dT|stability + dQe/dT|stability)
   %               + dQc/dT_sfc
   %
   % where dQc/dT_sfc = -k_eff(1) / (dz(1) / 2) is constant.
   %
   % This is the Dirichlet surface solve: Qc and its derivative enter the
   % Newton-Raphson residual and Jacobian directly. In the Robin path,
   % conduction instead enters through the top-node finite-difference
   % equation in `icemodel.column.assemble_enthalpy_system` rather than as
   % an explicit Jacobian term here.
   %
   % T_sfc on input is the outer coupling iterate from
   % solve_surface_energy_balance, used as the initial Newton guess. On
   % convergence T_sfc is the solution.
   %
   % See also: icemodel.surface.solve_surface_energy_balance,
   %           icemodel.surface.conductive_heat_flux,
   %           icemodel.surface.surface_energy_balance_residual
   %
   %#codegen

   % Solver tolerances.
   persistent tol maxiter
   if isempty(tol)
      tol = 1e-3;
      maxiter = 100;
   end

   % Physical constants and parameters for the analytical Jacobian.
   persistent cv_air cv_liq emiss SB epsilon
   if isempty(cv_air)
      [cv_air, cv_liq, SB, epsilon] = icemodel.physicalConstant( ...
         'cv_air', 'cv_liq', 'SB', 'epsilon');
      emiss = icemodel.parameterLookup('emiss');
   end

   % Gather T_sfc-independent terms in the SEB equation.
   Qa = icemodel.surface.advective_heat_flux(ppt, tppt, cv_liq);
   AAA = cv_air * De;                     % [W m-2 K-1]
   CCC = epsilon / psfc;                  % [Pa-1] = [m3 J-1]
   EEE = chi * (1.0 - albedo) * Qsi + emiss * Qli + Qa;
   FFF = roL * De;                        % [W m-2]

   % Initial value for dQc_dT_sfc which is constant over iterations.
   [~, dQc_dT_sfc] = icemodel.surface.conductive_heat_flux( ...
      k_eff, T_ice, dz, T_sfc);

   ok = false;
   old = T_sfc;
   for iter = 1:maxiter

      % Conductive heat flux at the current T_sfc iterate.
      Qc = icemodel.surface.conductive_heat_flux(k_eff, T_ice, dz, old);

      % Surface saturation vapor pressure and its temperature derivative.
      [es_sfc, des_sfc_dT] = icemodel.vapor.saturation_vapor_pressure( ...
         old, liqflag);

      % Bulk-Richardson stability factor and its temperature derivative.
      [stability, dstability] = ...
         icemodel.surface.turbulence.bulk_richardson.stability_factor( ...
         old, tair, wspd, br_coefs);

      % Newton-Raphson residual and analytical Jacobian.
      f = EEE - emiss * SB * old ^ 4 + AAA * (tair - old) * stability ...
         + FFF * CCC * (ea_atm - es_sfc) * stability + Qc;

      dfdT = -4.0 * emiss * SB * old ^ 3 ...
         + stability * -AAA + AAA * (tair - old) * dstability ...
         + stability * -FFF * CCC * des_sfc_dT ...
         + FFF * CCC * (ea_atm - es_sfc) * dstability ...
         + dQc_dT_sfc;

      T_sfc = old - f / dfdT;

      if abs(T_sfc - old) < tol
         ok = true;
         return
      elseif ~isfinite(T_sfc) || iter == maxiter
         ok = false;
         T_sfc = tair;
         return
      end
      old = T_sfc;
   end
end
