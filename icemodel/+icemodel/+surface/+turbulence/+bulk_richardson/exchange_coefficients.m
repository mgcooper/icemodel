function [De_h, S123, W1] = exchange_coefficients(wspd, z0_bulk, z_tair, z_wind)
   %EXCHANGE_COEFFICIENTS Compute bulk-Richardson exchange coefficients.
   %
   % Notation below roughly follows Liston et al. 1999
   %
   % Exchange coefficient [m s-1]:
   %  D_e = k ^ 2 * u / (log(z / z0)) ^ 2
   %      = c_w * u
   %
   % where c_w is the drag coefficient [1]:
   %  c_w = k ^ 2 / (log(z / z0)) ^ 2
   %
   % Thus: c_w = D_e / u
   %
   %  k = Von Karman constant
   %  u = wind speed
   %  z = z_tair = z_wind (assumed equal)
   %  z0 = roughness length for momentum
   %
   % Stability coefficient:
   %
   % If Ri < 0 (unstable):
   %     S = 1 - eta * Ri / ( 1 + gamma * abs(Ri) ^ 0.5)
   %
   % If Ri > 0 (stable):
   %     S = 1 / (1 + eta/2 * Ri) ^ 2
   %
   % Eliminate Ri so S can be computed in terms of Ts, Ta:
   %
   % Define:
   %  T_star = (Ts - Ta) / Ta
   %
   % With some algebra:
   %  if Ri < 0 (unstable):
   %     S = 1 + S2 / u ^ 2 * T_star / (1 + S3 / u * T_star ^ 0.5)
   %
   %  if Ri > 0 (stable):
   %     S = 1 / (1 + S2 / (2 * u ^ 2) * -T_star) ^ 2
   %
   % Where:
   %  S1 = psi * eta * c_w * (z_obs / z0_bulk) ^ 0.5  (eq 20)
   %  S2 = eta * g0 * z_obs
   %  S3 = S1 * (g * z_obs) ^ 0.5
   %
   % Translate to Louis, 1979 notation:
   %  c_w = a^2                                       (eq 13)
   %  eta = b                                         (see eq 14)
   %  psi = c_star                                    (see eq 20)
   %  S1 = c = c_star * a^2 * b * (z / z0_bulk)^0.5   (eq 20)
   %  S1 = gamma in Liston et al. 1999, eq. A15
   %
   % Note also that eta_star = eta / 2 = 9.4 / 2 = 4.7 is customarily referred
   % to as "beta" for the stable case, whereas gamma is comparable to the
   % "gamma" normally used for the unstable case.
   %
   %#codegen

   % Note: vectorized operations are used to support post-run diagnostics
   % where z_tair, z_wind, z0_bulk are allowed to vary in time. They remain
   % scalar-valued for the expected model initialization use-case.

   % Parse the optional z_wind input; set it equal to z_tair if not provided
   if nargin < 4
      z_wind = z_tair;
   end

   % Load physical constants and parameters
   [k_von_karman, gravity] = icemodel.physicalConstant('kappa', 'gravity');
   [eta, psi] = icemodel.parameterLookup('thf_bulk_richardson_eta', ...
      'thf_bulk_richardson_psi');

   % Compute the drag coefficient (a^2 in Louis notation, eq. 13) [1]
   W1 = (k_von_karman ./ log(z_tair ./ z0_bulk)) .^ 2;

   % Compute the aerodynamic exchange coefficient for sensible heat [m s-1]
   De_h = W1 .* wspd;

   % Compute the factor needed to allow z_wind != z_tair (derived by mgc)
   z_star = z_wind .^ 2 ./ z_tair;

   % These coefficients are used in the stability factor calculation and do
   % not depend on dynamic forcing. Pre-compute them here for computational
   % efficiency, to avoid additional operations within the solver.
   S1 = psi * eta .* W1 .* sqrt(z_tair ./ z0_bulk); % gamma Eq. A15
   S2 = eta * gravity .* z_star;
   S3 = S1 .* sqrt(gravity * z_star);
   if isscalar(S1)
      S123 = [S1, S2, S3];
   else
      S123 = [S1(:), S2(:), S3(:)];
   end

   % This follows Glen, and can be compared with stability_factor
   %
   % B1 = 9.4 * gravity * zobs / (Ta * wspd ^ 2) = S123(2) / (Ta * wspd ^ 2)
   % B2 = 5.3 * 9.4 * W1 * sqrt(zobs / z0) * sqrt(gravity * zobs / (Ta * wspd ^ 2));
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
