function [k, dkdT] = thermal_conductivity_air(T, reference)
   %THERMAL_CONDUCTIVITY_AIR Dry-air thermal conductivity and derivative.
   %
   %  [k, dkdT] = icemodel.kernels.thermal_conductivity_air(T)
   %  [k, dkdT] = icemodel.kernels.thermal_conductivity_air(T, reference)
   %
   % Description:
   %  Returns dry-air thermal conductivity k [W m-1 K-1] and its
   %  temperature derivative dkdT [W m-1 K-2] for numeric-array
   %  temperatures T [K]. The available options are temperature-only forms
   %  appropriate for dilute or near-atmospheric dry air, which is the
   %  relevant regime for snow and firn pore space.
   %
   % Inputs:
   %  T         - Temperature [K]
   %  reference - Conductivity formulation (default: "lemmon_jacobsen_2004")
   %
   % Options:
   %  "lemmon_jacobsen_2004" - Lemmon and Jacobsen (2004) dilute-gas air
   %                           formulation. This is the most rigorous
   %                           temperature-only option in the intended
   %                           snow/firn application regime.
   %  "sutherland"           - Sutherland-type engineering fit retained for
   %                           comparison and simple diagnostic use.
   %
   % Notes:
   %  Lemmon and Jacobsen (2004) give a full temperature-density transport
   %  formulation. This helper implements only the dilute-gas contribution
   %  because the current interface is temperature-only and the intended use
   %  is low-density pore air rather than dense or near-critical air.
   %
   % References:
   %
   %   Lemmon, E.W. and Jacobsen, R.T. (2004), Viscosity and Thermal
   %   Conductivity Equations for Nitrogen, Oxygen, Argon, and Air,
   %   International Journal of Thermophysics, 25(1), 21-69.
   %
   % See also: thermal_conductivity_ice, thermal_conductivity_water

   arguments
      T {mustBeNumeric} = 273.15
      reference (1, 1) string {mustBeMember(reference, ...
         ["lemmon_jacobsen_2004", "sutherland"])} = "lemmon_jacobsen_2004"
   end

   switch reference
      case "lemmon_jacobsen_2004"
         [k, dkdT] = lemmon_jacobsen_2004_dilute(T);

      case "sutherland"
         [k, dkdT] = sutherland_fit(T);
   end
end

function [k, dkdT] = lemmon_jacobsen_2004_dilute(T)
   % Lemmon and Jacobsen (2004) dilute-gas coefficients for air.

   collision_coeffs = [0.431, -0.4623, 0.08406, 0.005341, -0.00331];
   molar_mass = 28.9586;   % [g mol^-1]
   sigma = 0.360;          % Lennard-Jones size parameter [nm]
   eps_over_k = 103.3;     % Lennard-Jones energy parameter [K]
   Tc = 132.6312;          % Critical temperature [K]

   N1 = 1.308;
   N2 = 1.405;
   N3 = -1.036;
   t2 = -1.1;
   t3 = -0.3;

   log_Tstar = log(T ./ eps_over_k);
   omega = exp( ...
      collision_coeffs(1) ...
      + collision_coeffs(2) .* log_Tstar ...
      + collision_coeffs(3) .* log_Tstar .^ 2 ...
      + collision_coeffs(4) .* log_Tstar .^ 3 ...
      + collision_coeffs(5) .* log_Tstar .^ 4);

   dlogomega_dT = ( ...
      collision_coeffs(2) ...
      + 2 .* collision_coeffs(3) .* log_Tstar ...
      + 3 .* collision_coeffs(4) .* log_Tstar .^ 2 ...
      + 4 .* collision_coeffs(5) .* log_Tstar .^ 3) ./ T;

   % Dilute-gas viscosity in mPa s.
   eta0 = 0.0266958 .* sqrt(molar_mass .* T) ./ (sigma .^ 2 .* omega);
   deta0dT = eta0 .* (0.5 ./ T - dlogomega_dT);

   tau = Tc ./ T;
   lambda0_mW = N1 .* eta0 + N2 .* tau .^ t2 + N3 .* tau .^ t3;
   dlambda0_mWdT = N1 .* deta0dT ...
      - (N2 .* t2 .* tau .^ t2 + N3 .* t3 .* tau .^ t3) ./ T;

   k = 1e-3 .* lambda0_mW;
   dkdT = 1e-3 .* dlambda0_mWdT;
end

function [k, dkdT] = sutherland_fit(T)
   % Sutherland-type engineering fit.

   k0 = 0.0241;   % Reference conductivity [W m^-1 K^-1]
   T0 = 273.0;    % Reference temperature in published fit [K]
   Sk = 194.0;    % Sutherland constant for air [K]

   k = k0 .* (T ./ T0) .^ (3 / 2) .* ((T0 + Sk) ./ (T + Sk));

   % Using:
   % d/dT ln(k) = 3/(2T) - 1/(T + Sk)
   dkdT = k .* (1.5 ./ T - 1.0 ./ (T + Sk));
end
