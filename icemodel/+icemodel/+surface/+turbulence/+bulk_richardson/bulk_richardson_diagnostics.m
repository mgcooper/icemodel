function diag = bulk_richardson_diagnostics(T_sfc, es_sfc, tair, wspd, ...
      psfc, De, ea_atm, stability, ro_air_Lv, z0_bulk)
   %BULK_RICHARDSON_DIAGNOSTICS Assemble the full bulk-Richardson diagnostic struct.
   %
   %  diag = icemodel.surface.turbulence.bulk_richardson.bulk_richardson_diagnostics(
   %     T_sfc, es_sfc, tair, wspd, psfc, De, ea_atm, stability, ro_air_Lv, z0_bulk)
   %
   %  Assembles the complete diagnostic struct returned by turbulent_heat_flux
   %  when called with nargout > 2. This includes the standard bulk-Richardson
   %  fields (scheme, stability_factor, es_sfc, z0m, L) and the scalar-exchange
   %  experiment fields (z0h, z0q, u_star, Re, De_h, De_e, scalar_exchange_Qh,
   %  scalar_exchange_Qe) computed by the EXCOEFS_SCALAR roughness closures.
   %
   %  The scalar-exchange experiment uses separate scalar roughness lengths z0h
   %  (heat) and z0q (moisture) derived from the friction Reynolds number Re via
   %  the Andreas (2002) roughness closures. When the experiment is inactive
   %  (guard conditions fail), the scalar fluxes fall back to the standard De.
   %
   %  Inputs:
   %    T_sfc      - surface temperature [K]
   %    es_sfc     - saturation vapor pressure at the surface [Pa]
   %    tair       - air temperature [K]
   %    wspd       - wind speed [m s-1]
   %    psfc       - surface air pressure [Pa]
   %    De         - bulk exchange coefficient [m s-1]
   %    ea_atm     - atmospheric vapor pressure [Pa]
   %    stability  - Louis/Liston stability factor [-]
   %    ro_air_Lv  - air density × latent heat [J m-3]
   %    z0_bulk    - aerodynamic (momentum) roughness length [m]
   %
   %  Output:
   %    diag       - diagnostic struct (see field list below)
   %
   % See also:
   %   icemodel.surface.turbulence.bulk_richardson.turbulent_heat_flux
   %   icemodel.surface.turbulence.bulk_richardson.stability_factor
   %
   %#codegen

   persistent epsilon kappa nu_air cv_air
   if isempty(epsilon)
      [epsilon, kappa, cv_air] = icemodel.physicalConstant( ...
         'epsilon', 'kappa', 'cv_air');
      nu_air = 1.461e-5;
   end

   % Default scalar-exchange fields (inactive / fallback)
   scalar_active = false;
   z_obs   = NaN;
   u_star  = NaN;
   Re      = NaN;
   z0h     = NaN;
   z0q     = NaN;
   De_h    = De;
   De_e    = De;

   if isfinite(wspd) && isfinite(De) && isfinite(z0_bulk) ...
         && wspd > 0 && De > 0 && z0_bulk > 0

      Cd = De / wspd;
      if isfinite(Cd) && Cd > 0

         log_z_obs_over_z0 = kappa / sqrt(Cd);
         z_obs_try  = z0_bulk * exp(log_z_obs_over_z0);
         u_star_try = sqrt(De * wspd);
         Re_try     = u_star_try * z0_bulk / nu_air;
         [z0h_try, z0q_try] = scalar_roughness_lengths(z0_bulk, Re_try);

         if isfinite(z_obs_try) && isfinite(u_star_try) && isfinite(Re_try) ...
               && isfinite(z0h_try) && isfinite(z0q_try) ...
               && z_obs_try > z0h_try && z_obs_try > z0q_try

            scalar_active = true;
            z_obs   = z_obs_try;
            u_star  = u_star_try;
            Re      = Re_try;
            z0h     = z0h_try;
            z0q     = z0q_try;
            De_h = wspd * kappa ^ 2 / (log(z_obs / z0_bulk) * log(z_obs / z0h));
            De_e = wspd * kappa ^ 2 / (log(z_obs / z0_bulk) * log(z_obs / z0q));
         end
      end
   end

   Qh_scalar = cv_air    * De_h * stability * (tair - T_sfc);
   Qe_scalar = ro_air_Lv * De_e * stability * (epsilon / psfc * (ea_atm - es_sfc));

   diag = struct( ...
      'scheme',               'bulk_richardson', ...
      'stability_factor',     stability, ...
      'es_sfc',               es_sfc, ...
      'z0m',                  z0_bulk, ...
      'z0h',                  z0h, ...
      'z0q',                  z0q, ...
      'u_star',               u_star, ...
      'L',                    NaN, ...
      'Re',                   Re, ...
      'scalar_exchange_active', scalar_active, ...
      'scalar_exchange_z_obs',  z_obs, ...
      'scalar_exchange_De_h',   De_h, ...
      'scalar_exchange_De_e',   De_e, ...
      'scalar_exchange_Qh',     Qh_scalar, ...
      'scalar_exchange_Qe',     Qe_scalar);
end

% =========================================================================

function [z0h, z0q] = scalar_roughness_lengths(z0_bulk, Re)
   %SCALAR_ROUGHNESS_LENGTHS Andreas (2002) EXCOEFS_SCALAR roughness closures.

   if Re < 0.135
      z0h = z0_bulk * exp(1.250);
      z0q = z0_bulk * exp(1.610);
   elseif Re < 2.5
      log_Re = log(Re);
      z0h = z0_bulk * exp(0.149 - 0.550 * log_Re);
      z0q = z0_bulk * exp(0.351 - 0.628 * log_Re);
   else
      log_Re = log(Re);
      z0h = z0_bulk * exp(0.317 - 0.565 * log_Re - 0.183 * log_Re ^ 2);
      z0q = z0_bulk * exp(0.396 - 0.512 * log_Re - 0.180 * log_Re ^ 2);
   end
end
