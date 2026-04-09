function [T_sfc, ok] = solve_surface_temperature(tair, Qsi, Qli, albedo, ...
      wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, roL, liqflag, chi, varargin)
   %SOLVE_SURFACE_TEMPERATURE Solve the explicit bulk-Richardson SEB for Ts.
   %
   %  [T_sfc, ok] = icemodel.surface.solve_surface_temperature(...)
   %
   % Solves the nonlinear surface energy balance for T_sfc using a
   % Newton-Raphson iteration with an analytical Jacobian. The residual
   % being minimized is:
   %
   %   F(T_sfc) = chi*Qsi*(1-albedo) + emiss*(Qli - SB*T_sfc^4)
   %              + Qh(T_sfc) + Qe(T_sfc) + Qc(T_sfc) + Qa
   %
   % where Qh and Qe are evaluated using the bulk-Richardson turbulence
   % closure (Louis/Liston stability factor). The analytical Jacobian is:
   %
   %   dF/dT_sfc = -4*emiss*SB*T_sfc^3
   %               + d(stability)/dT * (Qh/stability + Qe/stability)
   %               + stability * (dQh/dT|stability + dQe/dT|stability)
   %               - a1   (only when conduction couples to the column)
   %
   % The analytical Jacobian is used in preference to the numerical
   % approach of evaluate_surface_flux for two reasons: (1) efficiency —
   % a single evaluation suffices per Newton step rather than two; and (2)
   % correctness of the conduction Jacobian term (-a1) when T_sfc feeds
   % back into Qc through the k_eff, T, dz column state.
   %
   % This function is specific to the bulk-Richardson turbulence scheme
   % because the analytical Jacobian requires explicit access to the
   % stability factor and its derivative. The outer SEB solver
   % solve_surface_energy_balance dispatches here for solver = 1.
   %
   % Conduction options (varargin):
   %
   %   1) solve_surface_temperature(..., Qc)
   %     Fixed conductive flux Qc. The Jacobian excludes the Qc-T_sfc coupling
   %     term (a1 = 0). Use when Qc is pre-computed and held fixed.
   %
   %   2) solve_surface_temperature(..., k_eff, T, dz)
   %     Conduction is treated as Qc = a1*(T(1) - T_sfc) where
   %     a1 = k_eff(1)/(dz(1)/2). The Jacobian includes the -a1 correction.
   %     Use when the outer coupling loop updates k_eff and T each iteration.
   %
   %   3) solve_surface_temperature(..., Qc, k_eff, T, dz, flag)
   %     Hybrid: when flag = true, uses the k_eff/T/dz form above (a1 != 0);
   %     when flag = false, uses the fixed Qc form (a1 = 0).
   %
   %  The default method is 2.
   %
   % See also: icemodel.surface.solve_surface_energy_balance,
   %           icemodel.surface.evaluate_surface_flux,
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

   % Parse conduction mode from varargin.
   switch numel(varargin)
      case 1
         % Qc included in the SEB but not the Jacobian.
         Qc = varargin{1};
         a1 = 0.0;
      case 3
         % Qc included in the SEB and the Jacobian.
         % Qc = a1 * (T(1) - T_sfc), so EEE += a1*T(1) and d(Qc)/dT_sfc -= a1.
         [k_eff, T, dz] = deal(varargin{:});
         a1 = k_eff(1) / (dz(1) / 2);
         Qc = 0.0;
      case 5
         [Qc, k_eff, T, dz, flag] = deal(varargin{:});
         if flag
            % Equivalent to case 3
            a1 = k_eff(1) / (dz(1) / 2);
            Qc = 0.0;
         else
            % Equivalent to case 1
            a1 = 0.0;
         end
      otherwise
         error('unrecognized number of inputs')
   end

   % Gather T_sfc-independent terms in the SEB equation.
   Qa = icemodel.surface.advective_heat_flux(ppt, tppt, cv_liq);
   AAA = cv_air * De;                     % [W m-2 K-1]
   CCC = epsilon / psfc;                  % [Pa-1] = [m3 J-1]
   EEE = chi * (1.0 - albedo) * Qsi + emiss * Qli + Qc + Qa;
   if a1 ~= 0.0
      EEE = EEE + a1 * T(1); % a1 * T(1) = Qc
   end
   FFF = roL * De;                        % [W m-2]

   T_sfc = nan;
   ok = false;
   old = tair;
   for iter = 1:maxiter

      % Surface saturation vapor pressure and its temperature derivative.
      [es_sfc, des_sfc_dT] = icemodel.vapor.saturation_vapor_pressure( ...
         old, liqflag);

      % Bulk-Richardson stability factor and its temperature derivative.
      [stability, dstability] = ...
         icemodel.surface.turbulence.bulk_richardson.stability_factor( ...
         old, tair, wspd, br_coefs);

      % Newton-Raphson residual and analytical Jacobian.
      f = EEE - emiss * SB * old ^ 4 + AAA * (tair - old) * stability ...
         + FFF * CCC * (ea_atm - es_sfc) * stability - a1 * old;

      dfdT = -4.0 * emiss * SB * old ^ 3 ...
         + stability * -AAA + AAA * (tair - old) * dstability ...
         + stability * -FFF * CCC * des_sfc_dT ...
         + FFF * CCC * (ea_atm - es_sfc) * dstability ...
         - a1;

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
