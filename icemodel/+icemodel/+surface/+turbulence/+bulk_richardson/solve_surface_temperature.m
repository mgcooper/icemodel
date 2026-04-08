function [T_sfc, ok] = solve_surface_temperature(tair, Qsi, Qli, albedo, ...
      wspd, ppt, tppt, psfc, De, ea_atm, chi, roL, br_coefs, liqflag, varargin)
   %SOLVE_SURFACE_TEMPERATURE Solve the explicit bulk-Richardson SEB for Ts.
   %
   % This function uses a Newton-Raphson iteration to solve the explicit
   % bulk-Richardson surface residual for Ts. The residual includes:
   %   absorbed shortwave
   %   net longwave
   %   sensible heat
   %   latent heat
   %   conductive heat
   %   precipitation-advected heat Qa
   %
   %#codegen

   % Default solver options
   persistent tol maxiter
   if isempty(tol)
      tol = 1e-3;
      maxiter = 100;
   end

   % Load physical constants and parameters
   persistent cv_air cv_liq emiss SB epsilon
   if isempty(cv_air)
      [cv_air, cv_liq, SB, epsilon] = icemodel.physicalConstant(...
         'cv_air', 'cv_liq', 'SB', 'epsilon');
      emiss = icemodel.parameterLookup('emiss');
   end

   % Parse inputs
   switch numel(varargin)
      case 1
         Qc = varargin{1};
         a1 = 0.0;
      case 3
         % Qc = a1 * T(1) - a1 * Ts;
         [k_eff, T, dz] = deal(varargin{:});
         a1 = k_eff(1) / (dz(1) / 2);
         Qc = 0.0;
      case 5
         % use this to test w/wo the Qc derivative
         [Qc, k_eff, T, dz, flag] = deal(varargin{:});
         if flag
            a1 = k_eff(1) / (dz(1) / 2);
            Qc = 0.0;
         else
            a1 = 0.0;
         end
      otherwise
         error('unrecognized number of inputs')
   end

   % Gather T_sfc-independent terms in the SEB equation.
   Qa = QADVECT(ppt, tppt, cv_liq); % [W m-2]
   AAA = cv_air * De;   % [W m-2 K-1]
   CCC = epsilon / psfc;  % [Pa-1] = [m3 J-1]
   EEE = chi * (1.0 - albedo) * Qsi + emiss * Qli + Qc + Qa; % [W m-2]
   if a1 ~= 0.0
      EEE = EEE + a1 * T(1);
   end
   FFF = roL * De; % [W m-2]

   T_sfc = nan;
   ok = false;
   old = tair;
   for iter = 1:maxiter

      % Surface saturation vapor pressure and derivative from icemodel.vapor.saturation_vapor_pressure.
      [es_sfc, des_sfc_dT] = icemodel.vapor.saturation_vapor_pressure(old, liqflag);

      % Bulk richardson stability factor and derivative.
      [stability, dstability] = ...
         icemodel.surface.turbulence.bulk_richardson.stability_factor( ...
         old, tair, wspd, br_coefs);

      % Evaluate the flux and derivative
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
