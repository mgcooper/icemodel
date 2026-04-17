function diag = bulk_richardson_diagnostics(T_sfc, es_sfc, tair, wspd, ...
      psfc, ea_atm, H_h, cv_atm, hv_atm, stability, z0_bulk)
   %BULK_RICHARDSON_DIAGNOSTICS Assemble the full bulk-Richardson diagnostics.
   %
   % diag = ...
   %    icemodel.surface.turbulence.bulk_richardson.bulk_richardson_diagnostics(...
   %    T_sfc, es_sfc, tair, wspd, psfc, ea_atm, H_h, cv_atm, hv_atm, ...
   %    stability, z0_bulk)
   %
   % Assembles the complete diagnostic struct returned by turbulent_heat_flux
   % when called with nargout > 2. This includes the standard bulk-Richardson
   % fields and a scalar-exchange experiment that adds separate roughness
   % lengths for heat (z0h) and moisture (z0q) via the Andreas (2002)
   % parameterization.
   %
   % De_h (the production neutral exchange coefficient) is reconstructed
   % internally from H_h and cv_atm: De_h = H_h / cv_atm, since
   % H_h = cv_atm * De_h by construction.
   %
   % --------------------------------------------------------------------------
   % SCALAR-EXCHANGE EXPERIMENT
   % --------------------------------------------------------------------------
   % The production bulk-Richardson scheme uses a single aerodynamic roughness
   % length for momentum, heat, and moisture (z0h = z0q = z0m = z0_bulk). Over
   % aerodynamically rough surfaces this overestimates scalar transfer, because
   % the turbulent eddies that dominate momentum exchange are less effective at
   % transporting heat and moisture. The scalar-exchange experiment corrects for
   % this by computing separate z0h and z0q from the friction Reynolds number
   % Re:
   %
   %   Re = u_* z0m / nu_air                 (roughness Reynolds number)
   %   u_* = sqrt(De_h * U)                  (friction velocity, neutral approx)
   %   [z0h, z0q] = scalar_roughness_lengths(z0m, Re)  (Andreas 2002)
   %
   %   De_h = U kappa^2 / (ln(z_obs/z0m) * ln(z_obs/z0h))
   %   De_e = U kappa^2 / (ln(z_obs/z0m) * ln(z_obs/z0q))
   %
   % where z_obs is the implied observation height from the neutral coefficient:
   %   ln(z_obs/z0m) = kappa / sqrt(Cd),   Cd = De_h / U
   %
   % In the aerodynamically rough regime (Re > 2.5) Andreas z0h and z0q are
   % significantly smaller than z0m (0.1-0.5x z0m), so De_h < De and De_e < De
   % thus the scalar-exchange scheme predicts weaker scalar fluxes than the
   % production BR scheme. In the aerodynamically smooth regime (Re < 0.135)
   % z0h and z0q exceed z0m, reversing the sign of the correction. For typical
   % glacier conditions (z0m ~ 1-3 mm, u* ~ 0.1-0.5 m/s), Re ~ 7-100, placing
   % most timesteps in the rough regime and explaining the systematic offset
   % toward smaller |Qe| and |Qh| relative to the production BR scheme.
   %
   % STABILITY NOTE: The Louis/Liston stability factor is re-used unchanged from
   % the production run. This is an approximation: the parameterization was
   % calibrated assuming z0h = z0m, so applying it with distinct scalar
   % roughness lengths is internally inconsistent. The error is second-order
   % relative to the roughness correction itself, but a fully consistent
   % treatment would require separate stability corrections for momentum and
   % scalars (as in the full Monin-Obukhov scheme).
   %
   % --------------------------------------------------------------------------
   %
   %  Inputs:
   %    T_sfc      - surface temperature [K]
   %    es_sfc     - saturation vapor pressure at the surface [Pa]
   %    tair       - air temperature [K]
   %    wspd       - wind speed [m s-1]
   %    psfc       - surface air pressure [Pa]
   %    ea_atm     - atmospheric vapor pressure [Pa]
   %    H_h        - sensible heat transport prefactor [W m-2 K-1]
   %    cv_atm     - volumetric heat capacity of moist air [J m-3 K-1]
   %    hv_atm     - volumetric latent enthalpy of moist air [J m-3]
   %    stability  - Louis/Liston stability factor [-]
   %    z0_bulk    - aerodynamic (momentum) roughness length [m]
   %
   %  Output:
   %    diag - diagnostic struct with fields:
   %      scheme                  - 'bulk_richardson'
   %      stability_factor        - Louis/Liston f(Ri_b) [-]
   %      es_sfc                  - saturation vapor pressure at surface [Pa]
   %      z0m                     - momentum roughness length [m]
   %      z0h                     - scalar heat roughness (Andreas 2002) [m]
   %      z0q                     - scalar moisture roughness (Andreas 2002) [m]
   %      u_star                  - friction velocity (neutral approx) [m s-1]
   %      L                       - Monin-Obukhov length [m] (NaN; BR does
   %                                not converge L)
   %      Re                      - roughness Reynolds number [-]
   %      scalar_exchange_active  - true if Andreas roughness lengths applied
   %      scalar_exchange_z_obs   - implied observation height [m]
   %      scalar_exchange_De_h    - heat exchange coefficient [m s-1]
   %      scalar_exchange_De_e    - moisture exchange coefficient [m s-1]
   %      scalar_exchange_Qh      - sensible heat flux (scalar exchange)
   %                                [W m-2]
   %      scalar_exchange_Qe      - latent heat flux (scalar exchange)
   %                                [W m-2]
   %
   % References:
   %   Andreas, E. L., 2002: Parameterizing scalar transfer over snow and
   %     ice: a review. J. Hydrometeor., 3, 417-432.
   %   Louis, J.-F., 1979: A parametric model of vertical eddy fluxes in
   %     the atmosphere. Boundary-Layer Meteorol., 17, 187-202.
   %
   % See also:
   %   icemodel.surface.turbulence.bulk_richardson.turbulent_heat_flux
   %   icemodel.surface.turbulence.bulk_richardson.stability_factor
   %   icemodel.surface.turbulence.monin_obukhov.turbulent_heat_flux
   %
   %#codegen

   persistent epsilon kappa
   if isempty(epsilon)
      [epsilon, kappa] = icemodel.physicalConstant('epsilon', 'kappa');
   end

   % Reconstruct the production neutral exchange coefficient from the
   % precomputed sensible heat transport prefactor.
   De_h = H_h / cv_atm;

   % Guard: if the production exchange coefficient is zero (calm air or
   % missing aerodynamic signal), the scalar-exchange experiment is
   % undefined. Return production fluxes (also zero) and mark inactive.
   if De_h == 0 || wspd == 0

      es_sfc = icemodel.vapor.saturation_vapor_pressure(T_sfc, false);

      diag = struct( ...
         'scheme',                  'bulk_richardson', ...
         'stability_factor',        stability, ...
         'es_sfc',                  es_sfc, ...
         'z0m',                     z0_bulk, ...
         'z0h',                     z0_bulk, ...
         'z0q',                     z0_bulk, ...
         'u_star',                  0.0, ...
         'L',                       NaN, ...
         'Re',                      0.0, ...
         'scalar_exchange_active',  false, ...
         'scalar_exchange_z_obs',   NaN, ...
         'scalar_exchange_De_h',    0.0, ...
         'scalar_exchange_De_e',    0.0, ...
         'scalar_exchange_Qh',      H_h .* stability .* (tair - T_sfc), ...
         'scalar_exchange_Qe',      hv_atm .* (De_h .* epsilon ./ psfc) ...
            .* stability .* (ea_atm - es_sfc));
      return
   end

   % Local moist-air density and kinematic viscosity.
   ro_atm = icemodel.vapor.moist_air_density(psfc, ea_atm, tair);
   nu_air = icemodel.kernels.air_kinematic_viscosity(tair, ro_atm);

   % Bulk drag coefficient from the neutral exchange coefficient.
   Cd = De_h / wspd;

   % z_obs is not available in the calling function, so recover it from
   % the neutral momentum exchange coefficient. In the production BR
   % scheme, De = U * (kappa / ln(z_obs/z0m))^2, so
   % ln(z_obs/z0m) = kappa/sqrt(Cd).
   log_z_obs_over_z0 = kappa / sqrt(Cd);
   z_obs = z0_bulk * exp(log_z_obs_over_z0);

   % Friction velocity and roughness Reynolds # (neutral approximation).
   u_star = sqrt(De_h * wspd);               % u_* = sqrt(Cd) * U
   Re = u_star * z0_bulk / nu_air;           % Re_* = u_* z0m / nu

   % Scalar roughness lengths from Andreas (2002).
   [z0h, z0q] = scalar_roughness_lengths(z0_bulk, Re);

   % Neutral scalar exchange coefficients with separate roughness
   % lengths for heat and moisture.
   %   De_h = U kappa^2 / (ln(z_obs/z0m) * ln(z_obs/z0h))
   %   De_e = U kappa^2 / (ln(z_obs/z0m) * ln(z_obs/z0q))
   De_h = wspd * kappa ^ 2 / (log(z_obs / z0_bulk) * log(z_obs / z0h));
   De_e = wspd * kappa ^ 2 / (log(z_obs / z0_bulk) * log(z_obs / z0q));

   % Turbulent heat fluxes using the separate exchange coefficients
   % and the production Louis/Liston stability factor. See the
   % stability note in the header for the approximation this entails.
   Qh = cv_atm * De_h * stability * (tair - T_sfc);
   Qe = hv_atm * De_e * stability * (epsilon / psfc * (ea_atm - es_sfc));

   diag = struct( ...
      'scheme', 'bulk_richardson', ...
      'stability_factor', stability, ...
      'es_sfc', es_sfc, ...
      'z0m', z0_bulk, ...
      'z0h', z0h, ...
      'z0q', z0q, ...
      'u_star', u_star, ...
      'L', NaN, ...
      'Re', Re, ...
      'scalar_exchange_active', true, ...
      'scalar_exchange_z_obs', z_obs, ...
      'scalar_exchange_De_h', De_h, ...
      'scalar_exchange_De_e', De_e, ...
      'scalar_exchange_Qh', Qh, ...
      'scalar_exchange_Qe', Qe);
end

function [z0h, z0q] = scalar_roughness_lengths(z0_bulk, Re)
   %SCALAR_ROUGHNESS_LENGTHS Andreas (2002) scalar roughness closures.
   %
   % Returns the scalar roughness lengths for heat (z0h) and moisture
   % (z0q) as a function of the momentum roughness z0m and the roughness
   % Reynolds number Re = u_* z0m / nu. Three regimes:
   %
   %   Smooth    Re < 0.135 : z0h/z0q > z0m
   %   Transition 0.135-2.5 : intermediate behavior
   %   Rough      Re > 2.5  : z0h/z0q < z0m
   %
   % Coefficients from Andreas (2002), Table 1 (snow/ice surface).
   % Reference: Andreas, E. L., 2002, J. Hydrometeor., 3, 417-432.

   if Re < 0.135          % aerodynamically smooth
      z0h = z0_bulk * exp(1.250);
      z0q = z0_bulk * exp(1.610);
   elseif Re < 2.5        % transitional
      log_Re = log(Re);
      z0h = z0_bulk * exp(0.149 - 0.550 * log_Re);
      z0q = z0_bulk * exp(0.351 - 0.628 * log_Re);
   else                   % aerodynamically rough
      log_Re = log(Re);
      z0h = z0_bulk * exp(0.317 - 0.565 * log_Re - 0.183 * log_Re ^ 2);
      z0q = z0_bulk * exp(0.396 - 0.512 * log_Re - 0.180 * log_Re ^ 2);
   end
end
