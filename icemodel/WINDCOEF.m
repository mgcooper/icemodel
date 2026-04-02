function [De, S123, W1] = WINDCOEF(wspd, z0_bulk, z_tair, z_wind)
   %WINDCOEF Compute wind transfer coefficients
   %
   % Exchange coefficient:
   %  D_e = k ^ 2 * u / (log(z / z0)) ^ 2
   %      = c_w * u
   %
   %  where:
   %  c_w = k ^ 2 / (log(z / z0)) ^ 2
   %
   % Thus: c_w = D_e / u
   %
   % Stability coefficient:
   %
   % If Ri < 0 (unstable):
   %     S = 1 - eta * Ri / ( 1 + gamma * abs(Ri) ^ 0.5)
   %
   % If Ri > 0 (stable):
   %     S = 1 / (1 + eta/2 * Ri) ^ 2
   %
   % Define:
   %  T_star = (Ts - Ta) / Ta
   %
   % With some algebra:
   %  If Ri < 0 (unstable):
   %     S = 1 + S2 / u ^ 2 * T_star / (1 + S3 / u * T_star ^ 0.5)
   %
   %  If Ri > 0 (stable):
   %     S = 1 / (1 + S2 / (2 * u ^ 2) * -T_star) ^ 2
   %
   % Where:
   %  S1 = psi * eta * c_w * (z_obs / z0_bulk) ^ 0.5
   %  S2 = eta * g0 * z_obs
   %  S3 = S1 * (g * z_obs) ^ 0.5
   %
   % Note, from Louis, 1979:
   %  c_w = a^2                                    (eq 13)
   %  eta = b                                      (see eq 14)
   %  psi = c_star                                 (see eq 20)
   %  S1 = c = c_star * a^2 * b * (z / z0_bulk)^0.5 (eq 20)
   %  S1 = gamma in Liston et al. 1999, eq. A15
   %
   % Note also that eta_star = eta / 2 = 9.4 / 2 = 4.7 is customarily referred
   % to as "beta" for the stable case, whereas gamma is comparable to the
   % "gamma" normally used for the unstable case.
   %
   %#codegen

   [kappa, gravity] = icemodel.physicalConstant('kappa', 'gravity');
   [eta, psi] = icemodel.parameterLookup('thf_bulk_richardson_eta', ...
      'thf_bulk_richardson_psi');

   if nargin < 4
      z_wind = z_tair;
   end

   z_star = (z_wind ^ 2) / z_tair;

   W1 = (kappa / log(z_tair / z0_bulk)) ^ 2;
   De = wspd .* W1;

   S123 = nan(1, 3);
   S123(1) = psi * eta * W1 * sqrt(z_tair / z0_bulk); % gamma Eq. A15
   S123(2) = eta * gravity * z_star;
   S123(3) = S123(1) * sqrt(gravity * z_star);

   % This follows Glen, and can be compared with STABLEFN
   %
   % B1 = 9.4 * gravity * zobs / (Ta * wspd ^ 2) = scoef(2) / (Ta * wspd ^ 2)
   % B2 = 5.3 * 9.4 * wcoef * sqrt(zobs / z0) * sqrt(gravity * zobs / (Ta * wspd ^ 2));
   %
   % a1 = 5.3 * 9.4;                                    % [-]
   % z1 = z_obs / z0_bulk;                              % [-]
   % C1 = a1 * (kappa / log(z1)) ^ 2 * sqrt(z1);        % [-]
   % C2 = gravity * z_obs / (Ta * wspd ^ 2);            % [K-1]
   % B1 = 9.4 * C2;                                     % [K-1]
   % B2 = C1 * sqrt(C2);                                % [K-1]
   %
   % if (Ts > Ta)                                       % Unstable case.
   %    B3 = 1.0 + B2 * sqrt(Ts - Ta);                  % [-]
   %    S = 1.0 + B1 * (Ts - Ta) / B3;                  % [-]
   % elseif (Ts < Ta)                                   % Stable case.
   %    B8 = B1 / 2.0;                                  % [-]
   %    S = 1.0 / (1.0 + B8 * (Ta - Ts)) ^ 2;           % [-]
   % else                                               % Neutrally stable case.
   %    S = 1.0;                                        % [-]
   % end
end
