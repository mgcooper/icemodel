function [Ts, ok] = SFCTEMP(Ta, Qsi, Qli, albedo, wspd, ppt, tppt, Pa, ...
      De, ea, chi, roL, scoef, liqflag, varargin)
   %SFCTEMP Solve the explicit bulk-Richardson SEB for surface temperature.
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
   persistent Tf cv_air cv_liq emiss SB epsilon
   if isempty(Tf)
      [Tf, cv_air, cv_liq, SB, epsilon] = icemodel.physicalConstant(...
         "Tf", "cv_air", "cv_liq", "SB", "epsilon");
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

   % Gather Ts-independent terms in the SEB equation.
   Qa = QADVECT(ppt, tppt, cv_liq); % [W m-2]
   AAA = cv_air * De;   % [W m-2 K-1]
   CCC = epsilon / Pa;  % [Pa-1] = [m3 J-1]
   EEE = chi * (1.0 - albedo) * Qsi + emiss * Qli + Qc + Qa; % [W m-2]
   if a1 ~= 0.0
      EEE = EEE + a1 * T(1);
   end
   FFF = roL * De; % [W m-2]

   Ts = nan;
   ok = false;
   old = Ta;
   for iter = 1:maxiter

      % Saturation vapor pressure and derivative from VAPPRESS.
      [es, des_dT] = VAPPRESS(old, liqflag);

      % Bulkk richardson stability factor and derivative.
      [S, dS] = STABLEFN(Ta, old, wspd, scoef);

      % Evaluate the flux and derivative
      f = EEE - emiss * SB * old ^ 4 + AAA * (Ta - old) * S ...
         + FFF * CCC * (ea - es) * S - a1 * old;

      dfdT = -4.0 * emiss * SB * old ^ 3 ...
         + S * -AAA + AAA * (Ta - old) * dS ...
         + S * -FFF * CCC * des_dT ...
         + FFF * CCC * (ea - es) * dS ...
         - a1;

      Ts = old - f / dfdT;

      if abs(Ts - old) < tol
         ok = true;
         return
      elseif isnan(Ts) || iter == maxiter
         ok = false;
         Ts = Ta;
         return
      end
      old = Ts;
   end
end
