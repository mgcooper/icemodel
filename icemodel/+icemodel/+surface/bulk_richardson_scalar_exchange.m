function diag_scalar = bulk_richardson_scalar_exchange(T_sfc, es_sfc, tair, ...
      wspd, psfc, De, ea_atm, stability, roL, z0_bulk)
   %BULK_RICHARDSON_SCALAR_EXCHANGE Diagnose a scalar-exchange alternative.
   %
   %  diag_scalar = icemodel.surface.bulk_richardson_scalar_exchange(...)
   %
   % This helper implements a recommendation-supporting experiment based on
   % the snowphysics EXCOEFS_SCALAR idea: keep the current bulk-Richardson
   % stability factor, but let heat and moisture use scalar roughness
   % lengths smaller than the momentum roughness over rough surfaces.
   %
   % The production bulk_richardson scheme still uses the legacy single
   % exchange coefficient De. This helper only diagnoses the alternative
   % scalar exchange implied by the current De contract.
   %
   %#codegen

   persistent epsilon kappa nu_air cv_air
   if isempty(epsilon)
      [epsilon, kappa, cv_air] = icemodel.physicalConstant( ...
         'epsilon', 'kappa', 'cv_air');
      nu_air = 1.461e-5;
   end

   diag_scalar = struct( ...
      'active', false, ...
      'z_obs', NaN, ...
      'u_star', NaN, ...
      'Re', NaN, ...
      'z0h', NaN, ...
      'z0q', NaN, ...
      'De_h', De, ...
      'De_e', De, ...
      'Qh', NaN, ...
      'Qe', NaN);

   if ~(isfinite(wspd) && isfinite(De) && isfinite(z0_bulk)) ...
         || wspd <= 0 || De <= 0 || z0_bulk <= 0
      diag_scalar.Qh = cv_air * De * stability * (tair - T_sfc);
      diag_scalar.Qe = roL * De * stability ...
         * (epsilon / psfc * (ea_atm - es_sfc));
      return
   end

   Cd = De / wspd;
   if ~(isfinite(Cd) && Cd > 0)
      diag_scalar.Qh = cv_air * De * stability * (tair - T_sfc);
      diag_scalar.Qe = roL * De * stability ...
         * (epsilon / psfc * (ea_atm - es_sfc));
      return
   end

   log_z_obs_over_z0 = kappa / sqrt(Cd);
   z_obs = z0_bulk * exp(log_z_obs_over_z0);
   u_star = sqrt(De * wspd);
   Re = u_star * z0_bulk / nu_air;
   [z0h, z0q] = scalar_roughness_lengths(z0_bulk, Re);

   if ~(isfinite(z_obs) && isfinite(u_star) && isfinite(Re) && ...
         isfinite(z0h) && isfinite(z0q) && z_obs > z0h && z_obs > z0q)
      diag_scalar.Qh = cv_air * De * stability * (tair - T_sfc);
      diag_scalar.Qe = roL * De * stability ...
         * (epsilon / psfc * (ea_atm - es_sfc));
      return
   end

   De_h = wspd * kappa ^ 2 / (log(z_obs / z0_bulk) * log(z_obs / z0h));
   De_e = wspd * kappa ^ 2 / (log(z_obs / z0_bulk) * log(z_obs / z0q));

   diag_scalar.active = true;
   diag_scalar.z_obs = z_obs;
   diag_scalar.u_star = u_star;
   diag_scalar.Re = Re;
   diag_scalar.z0h = z0h;
   diag_scalar.z0q = z0q;
   diag_scalar.De_h = De_h;
   diag_scalar.De_e = De_e;
   diag_scalar.Qh = cv_air * De_h * stability * (tair - T_sfc);
   diag_scalar.Qe = roL * De_e * stability ...
      * (epsilon / psfc * (ea_atm - es_sfc));
end

function [z0h, z0q] = scalar_roughness_lengths(z0_bulk, Re)
   %SCALAR_ROUGHNESS_LENGTHS Apply the EXCOEFS_SCALAR roughness closures.

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
