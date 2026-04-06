function [Qe, Qh, diag] = turbulent_heat_flux_bulk_mo(T_sfc, tair, wspd, ...
      psfc, ea_atm, ro_sfc, snow_depth, roL, liqflag, opts)
   %TURBULENT_HEAT_FLUX_BULK_MO Evaluate the bulk Monin-Obukhov THF scheme.
   %
   %  [Qe, Qh] = icemodel.surface.turbulent_heat_flux_bulk_mo(...)
   %  [Qe, Qh, diag] = icemodel.surface.turbulent_heat_flux_bulk_mo(...)
   %
   % This function implements the glacier-oriented bulk Monin-Obukhov
   % turbulent-flux closure used by the `bulk_mo` scheme (van As et al. 2005).
   % The governing pieces are:
   %
   %  1. Momentum exchange
   %     u_* = kappa * U / (ln(z_u / z0m) - psi_m(z_u/L) + psi_m(z0m/L))
   %
   %  2. Scalar exchange
   %     theta_* = kappa * (theta_air - theta_sfc) / ...
   %        (ln(z_t / z0h) - psi_h(z_t/L) + psi_h(z0h/L))
   %     q_* = kappa * (q_air - q_sfc) / ...
   %        (ln(z_q / z0q) - psi_h(z_q/L) + psi_h(z0q/L))
   %
   %  3. Monin-Obukhov length
   %     L = u_*^2 theta_v / (kappa g theta_v*)
   %
   %  4. Fluxes
   %     Qh = rho_atm c_p u_* theta_*
   %     Qe = rho_atm L_v/s u_* q_*
   %
   % The scalar roughness lengths use Andreas (2002) over snow/firn and
   % Smeets and van den Broeke (2008) over rough bare ice. Stable profile
   % corrections use Holtslag and de Bruin; unstable corrections use the
   % Paulson/Dyer forms. Surface vapor pressure remains phase-aware through
   % the existing liqflag/roL contract.
   %
   % The implementation follows van As et al., 2005:
   %     van As et al., 2005, The Summer Surface Energy Balance of the High
   %     Antarctic Plateau, Boundary-Layer Meteorology.
   %
   %#codegen

   % Load physical constants
   persistent P0 kappa_p ro_air_ref cv_air
   if isempty(P0)
      [P0, kappa_p, ro_air_ref, cv_air] = icemodel.physicalConstant( ...
         'P0', 'kappa_p', 'ro_air', 'cv_air');
   end

   % Load parameters
   persistent wspd_min L_initial L_tol iter_max
   if isempty(P0)
      [wspd_min, L_initial, L_tol, iter_max] = ...
         icemodel.parameterLookup('thf_bulk_ws_min', 'thf_bulk_L_initial', ...
         'thf_bulk_L_tol', 'thf_bulk_iter_max');
   end

   % Observation heights and surface roughness lengths.
   z_q = opts.z_relh;
   z_t = opts.z_tair;
   z_u = opts.z_wind;
   z0m = icemodel.surface.surface_roughness_length(snow_depth, ro_sfc, ...
      opts.z0_ice, opts.z0_snow_low_density, opts.z0_snow_high_density);

   % Atmospheric and surface thermodynamic state. The repo stores air heat
   % capacity and latent energies in volumetric form elsewhere in the SEB
   % stack, so divide by a reference air density to recover the specific
   % quantities c_p [J kg^-1 K^-1] and L_v/s [J kg^-1] needed by the classic
   % flux formulas rho_atm * c_p * u_* * theta_* & rho_atm * L_v/s * u_* * q_*.
   es_sfc = VAPPRESS(T_sfc, liqflag);
   ro_atm = icemodel.kernels.moist_air_density(psfc, ea_atm, tair);
   nu_air = icemodel.kernels.air_kinematic_viscosity(tair, ro_atm);
   cp_air = cv_air / ro_air_ref;
   Le_air = roL / ro_air_ref;

   q_air = icemodel.kernels.specific_humidity_from_vapor_pressure(ea_atm, psfc);
   q_sfc = icemodel.kernels.specific_humidity_from_vapor_pressure(es_sfc, psfc);
   theta_air = tair * (P0 / psfc) ^ kappa_p;
   theta_sfc = T_sfc * (P0 / psfc) ^ kappa_p;

   % Return early if wind speed is below the minimum threshold. This avoids
   % unstable log/roughness diagnostics when the surface layer is effectively
   % decoupled.
   if wspd <= wspd_min
      Qe = 0.0;
      Qh = 0.0;
      if nargout > 2
         diag = bulk_mo_diag(es_sfc, q_air, q_air, theta_air, theta_sfc, ...
            ro_atm, nu_air, z0m, NaN, NaN, 0.0, NaN, 0.0, 0, 0, 0, 0, 0, ...
            0, 0, real(theta_air) >= real(theta_sfc));
      end
      return
   end

   use_snow_closure = snow_depth > 0;
   is_stable = real(theta_air) >= real(theta_sfc);

   % Start from neutral corrections and iterate Monin-Obukhov stability.
   % L_initial is a numerical initialization, not a physical tuning target:
   % it seeds the fixed-point solve from an almost-neutral surface layer and
   % is forgotten once the iteration converges.
   psi_m0 = 0.0;
   psi_mz = 0.0;
   psi_h0 = 0.0;
   psi_hz = 0.0;
   psi_q0 = 0.0;
   psi_qz = 0.0;
   z0h = NaN;
   z0q = NaN;
   u_star = 0.0;
   Re = 0.0;
   L = L_initial;

   for n_iter = 1:iter_max

      % Initial momentum exchange scale (friction velocity) from observed wind
      % speed and the current momentum roughness and stability corrections.


      % Roughness Reynolds number from u_*.


      % Initial scalar roughness lengths from Re.
      [z0h, z0q] = icemodel.surface.scalar_roughness_bulk_mo( ...
         wspd, z0m, psi_m0, psi_mz, nu_air, z_u, use_snow_closure);

      % Profile corrections for Monin-Obukhov stability.
      [psi_m0, psi_mz, psi_h0, psi_hz, psi_q0, psi_qz] = ...
         profile_corrections_bulk_mo(L, z0m, z_u, z0h, z_t, z0q, z_q, is_stable);

      % Recompute scalar roughness, u_* and Reynolds roughness with updated
      % momentum correction.
      [z0h, z0q, u_star, Re] = icemodel.surface.scalar_roughness_bulk_mo( ...
         wspd, z0m, psi_m0, psi_mz, nu_air, z_u, use_snow_closure);

      % The scalar transfer relations define theta_* and q_* once the scalar
      % roughness lengths and profile corrections are known.
      theta_star = exchange_scale(theta_air - theta_sfc, z_t, z0h, psi_hz, psi_h0);
      q_star = exchange_scale(q_air - q_sfc, z_q, z0q, psi_qz, psi_q0);

      % Convert the exchange scales to physical sensible and latent heat
      % fluxes. ro_atm is the local moist-air density, so using ro_atm here
      % is more consistent than assuming the fixed reference density used by
      % the older bulk-Richardson path.
      Qe = ro_atm * Le_air * u_star * q_star;
      if nargout > 1
         Qh = ro_atm * cp_air * u_star * theta_star;
      end

      % Update the Monin-Obukhov stability length from the current exchange
      % state. L depends on the sensible/moisture exchange scales, so this
      % fixed-point loop iterates the roughness, profile corrections, and
      % stability length together until the stability state is self-consistent.
      L_prev = L;
      L = icemodel.kernels.monin_obukhov_length(u_star, theta_air, q_air, ...
         theta_star, q_star);

      % Converge on a relative change in |L| rather than an absolute meter
      % threshold. This avoids temperature-sensitive jumps where one call
      % stops after N iterations and the next stops after N+1.
      L_scale = max([abs(L_prev), abs(L), 1e-12]);
      if abs(L) < 1e-12 || abs(L_prev - L) < L_tol * L_scale
         break
      end
   end

   if nargout > 2
      % Optional diagnostics expose the converged Monin-Obukhov state for
      % tests, interactive inspection, and debug-file dumps.
      diag = bulk_mo_diag(es_sfc, q_air, q_sfc, theta_air, theta_sfc, ...
         ro_atm, nu_air, z0m, z0h, z0q, u_star, L, Re, psi_m0, psi_mz, ...
         psi_h0, psi_hz, psi_q0, psi_qz, n_iter, is_stable);
   end
end

function theta_or_q_star = exchange_scale(delta_value, z_obs, z0_scalar, ...
      psi_scalar_z, psi_scalar_0)
   %EXCHANGE_SCALE Return theta_* or q_* from the scalar transfer relation.
   %
   % The momentum (wind) transfer relation is:
   %   u_* = kappa * delta_wind / ...
   %     (ln(z_wind / z0m) - psi_m(z_wind/L) + psi_m(z0m/L))
   %
   % Here delta_wind is wspd_z_obs - wspd_sfc, with wspd_sfc = 0 and z_wind the
   % wind speed observation height.
   %
   % The scalar transfer relation is:
   %   scalar_* = kappa * delta_scalar / ...
   %      (ln(z_obs / z0_scalar) - psi_scalar(z_obs/L) + psi_scalar(z0/L))
   %
   % Here delta_scalar is either theta_air - theta_sfc or q_air - q_sfc, z_obs
   % is the scalar observation height, z0_scalar is the temperature or
   % humidity roughness length, and the psi terms are the Monin-Obukhov
   % profile corrections for the relevant scalar.

   persistent kappa
   if isempty(kappa)
      kappa = icemodel.physicalConstant('kappa');
   end

   denom = log(z_obs / z0_scalar) - psi_scalar_z + psi_scalar_0;
   if abs(denom) < 1e-12
      denom = denom + icemodel.kernels.sign_or_one(denom) * 1e-12;
   end

   theta_or_q_star = kappa * delta_value / denom;
end

function [psi_m0, psi_mz, psi_h0, psi_hz, psi_q0, psi_qz] = ...
      profile_corrections_bulk_mo(L, z0m, z_u, z0h, z_t, z0q, z_q, is_stable)
   %PROFILE_CORRECTIONS_BULK_MO Return Monin-Obukhov profile corrections.
   %
   % Stable conditions use the Holtslag and de Bruin correction for both
   % momentum and scalars. Unstable conditions use Paulson/Dyer with separate
   % forms for momentum and scalars. The returned values are evaluated both at
   % the observation heights and at the roughness lengths because the flux
   % denominator uses psi(z/L) - psi(z0/L).

   if is_stable
      psi_m0 = icemodel.surface.psi_holtslag(z0m / L);
      psi_mz = icemodel.surface.psi_holtslag(z_u / L);
      psi_h0 = icemodel.surface.psi_holtslag(z0h / L);
      psi_hz = icemodel.surface.psi_holtslag(z_t / L);
      psi_q0 = icemodel.surface.psi_holtslag(z0q / L);
      psi_qz = icemodel.surface.psi_holtslag(z_q / L);
   else
      psi_m0 = icemodel.surface.psi_m_paulson(z0m / L);
      psi_mz = icemodel.surface.psi_m_paulson(z_u / L);
      psi_h0 = icemodel.surface.psi_h_paulson(z0h / L);
      psi_hz = icemodel.surface.psi_h_paulson(z_t / L);
      psi_q0 = icemodel.surface.psi_h_paulson(z0q / L);
      psi_qz = icemodel.surface.psi_h_paulson(z_q / L);
   end
end

function diag = bulk_mo_diag(es_sfc, q_air, q_sfc, theta_air, theta_sfc, ...
      ro_atm, nu_air, z0m, z0h, z0q, u_star, L, Re, psi_m0, psi_mz, ...
      psi_h0, psi_hz, psi_q0, psi_qz, n_iter, is_stable)
   %BULK_MO_DIAG Assemble optional diagnostics for debug/test use.
   %
   % The returned struct mirrors the most important intermediate state of the
   % bulk-MO solve: thermodynamic state, roughness lengths, friction velocity,
   % Monin-Obukhov length, profile corrections, and the number of fixed-point
   % iterations taken to converge.

   diag = struct( ...
      'scheme', 'bulk_mo', ...
      'es_sfc', es_sfc, ...
      'q_air', q_air, ...
      'q_sfc', q_sfc, ...
      'theta_air', theta_air, ...
      'theta_sfc', theta_sfc, ...
      'ro_atm', ro_atm, ...
      'nu_air', nu_air, ...
      'z0m', z0m, ...
      'z0h', z0h, ...
      'z0q', z0q, ...
      'u_star', u_star, ...
      'L', L, ...
      'Re', Re, ...
      'psi_m0', psi_m0, ...
      'psi_mz', psi_mz, ...
      'psi_h0', psi_h0, ...
      'psi_hz', psi_hz, ...
      'psi_q0', psi_q0, ...
      'psi_qz', psi_qz, ...
      'is_stable', is_stable, ...
      'n_iterations', n_iter);
end
