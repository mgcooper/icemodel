function [es, des_dT, ro_vap, dro_vapdT, k_vap, T_dew] = buckVaporModel( ...
      T, Tf, liqflag, rh)
   %BUCKVAPORMODEL Buck (1981) vapor model — self-contained reference archive.
   %
   %  Self-contained implementation of the Buck (1981) empirical vapor model,
   %  including saturation vapor pressure, vapor density, vapor thermal
   %  diffusion coefficient, and dew point. This file encodes everything
   %  needed to run the vapor subsystem using Buck instead of Ambaum/Romps.
   %
   %  ES = buckVaporModel(T, TF, LIQFLAG) computes saturation vapor pressure
   %  using the Buck (1981) empirical formula:
   %
   %     ES = A * exp(B * (T - Tf) / (C + T - Tf))   [Pa]
   %
   %  where (al, bl, cl) are over liquid and (ai, bi, ci) are over ice, and
   %  Tf = 273.16 is the triple point temperature.
   %
   %  This is an empirical fit to experimental data, distinct from the
   %  Rankine-Kirchhoff form es = a * exp(b/T) * T^c used in production
   %  (see icemodel.vapor.saturation_vapor_pressure). Buck's formula can equivalently be written in absolute
   %  Kelvin as es = A * exp(B*T/(C'+T)) with C' = C + Tf, but the original
   %  Celsius-offset form is retained here to match the published coefficients.
   %
   %  Note: temperatures for Buck's equations are in deg C, and vapor
   %  pressures are in mb. Here they have been converted so the calculations
   %  are done with temperatures in K, and vapor pressures in Pa. (1 mb =
   %  100 Pa). Curves ew1 and ei2 are adopted for water and ice, respectively
   %  (See Buck Eq. 3,8 and Table 2).
   %
   %  [ES, DES_DT] = buckVaporModel(...) Also computes the derivative of
   %  saturation vapor pressure with respect to temperature.
   %
   %  [ES, DES_DT, RO_VAP, DRO_VAPDT] = buckVaporModel(...) Also computes
   %  saturation vapor density and its temperature derivative.
   %
   %  [ES, DES_DT, RO_VAP, DRO_VAPDT, K_VAP] = buckVaporModel(...) Also
   %  computes the vapor thermal diffusion coefficient:
   %     k_vap = De * Ls * dro_vapdT
   %  where De = De0 * (T/Tf)^nd is the vapor diffusivity.
   %
   %  [..., T_DEW] = buckVaporModel(T, TF, LIQFLAG, RH) Also computes the
   %  dew point temperature in Kelvins from relative humidity [%].
   %
   %  This is the archived Buck (1981) implementation. Production code uses
   %  Ambaum (2020) / Romps (2021) via icemodel.vapor.saturation_vapor_pressure / icemodel.parameterLookup.
   %
   %  Derivative chain summary (Buck):
   %
   %     es       = a * exp(b*(T-Tf) / (c+T-Tf))
   %     des_dT   = b*c*es / (c+T-Tf)^2
   %     ro_vap   = es / (Rv*T)
   %     dro_vapdT = ro_vap * (b*c/(c+T-Tf)^2 - 1/T)
   %              = (des_dT - es/T) / (Rv*T)
   %     k_vap    = De * Ls * dro_vapdT
   %              = De*Ls/(Rv*T) * (des_dT - es/T)
   %
   % See also: icemodel.vapor.saturation_vapor_pressure, icemodel.parameterLookup

   % Buck (1981) coefficients (i = ice, l = liquid)
   al = 611.21;
   bl = 17.502;  % 18.678
   cl = 240.97;
   ai = 611.15;
   bi = 22.452;
   ci = 272.55;

   % Physical constants for vapor density and k_vap
   Rv = 461.5;       % Specific gas constant for water vapor [J kg-1 K-1]
   Ls = 2.834e6;     % Latent heat of sublimation [J kg-1]
   Lv = 2.501e6;     % Latent heat of vaporization [J kg-1]
   De0 = 9e-5;       % Reference vapor diffusivity [m2 s-1]
   nd = 14;          % Temperature exponent (Anderson 1976)

   if nargin < 2
      Tf = 273.16;
   end
   if nargin < 3
      liqflag = false;
   end

   % Select phase coefficients
   if liqflag
      a = al; b = bl; c = cl; L = Lv;
   else
      a = ai; b = bi; c = ci; L = Ls;
   end

   % Saturation vapor pressure [Pa]
   Td = T - Tf;
   es = a * exp(b * Td ./ (c + Td));

   % First derivative [Pa K-1]
   if nargout > 1
      des_dT = b * c * es ./ (c + Td) .^ 2;
   end

   % Vapor density [kg m-3]
   if nargout > 2
      ro_vap = es ./ (Rv * T);
   end

   % Vapor density derivative [kg m-3 K-1]
   if nargout > 3
      dro_vapdT = ro_vap .* (b * c ./ (c + Td) .^ 2 - 1 ./ T);
      % Equivalently: (des_dT - es ./ T) ./ (Rv * T)
   end

   % Vapor thermal diffusion coefficient [W m-1 K-1]
   if nargout > 4
      De = De0 * (T / Tf) .^ nd;
      k_vap = L * De .* dro_vapdT;
      % Equivalently: Ls * De ./ (Rv * T) .* (des_dT - es ./ T);
   end

   % Dew point temperature [K]
   if nargout > 5
      if nargin < 4
         error('buckVaporModel:missingRH', ...
            'Relative humidity (rh) required for dew point computation.');
      end
      e_air = es .* rh / 100;
      T_dew = Tf + c * log(e_air / a) ./ (b - log(e_air / a));
   end
end

%% Notes
%
% To convert from e_sat to T:
%
% T = Tf + c * log(e_sat / a) / (b - log(e_sat / a));
%
% Thus to compute Tdew, substitute e_air for e_sat:
%
% T = Tf + c * log(e_air / a) / (b - log(e_air / a));
%
% Buck suggests ew1, ei2, fw3, fi3. Below are expressions for f, not used
% here, but for reference:
% fw3 = 1.0007 + 3.46e-6 * P; % note: P is in mb, need to convert.
% fi3 = 1.0003 + 4.18e-6 * P;
%
% These enhancement factors are not implemented.
