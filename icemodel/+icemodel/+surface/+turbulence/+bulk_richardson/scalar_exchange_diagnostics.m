function diag_scalar = scalar_exchange_diagnostics(T_sfc, es_sfc, tair, ...
      wspd, psfc, De, ea_atm, stability, roL, z0_bulk)
   %SCALAR_EXCHANGE_DIAGNOSTICS Diagnose a scalar-exchange alternative.
   %
   %  diag_scalar = ...
   %     icemodel.surface.turbulence.bulk_richardson.scalar_exchange_diagnostics(...)
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

   state_size = size(De + 0 * T_sfc + 0 * tair + 0 * wspd + 0 * psfc ...
      + 0 * ea_atm + 0 * stability + 0 * roL + 0 * z0_bulk);
   T_sfc = expand_input(T_sfc, state_size);
   es_sfc = expand_input(es_sfc, state_size);
   tair = expand_input(tair, state_size);
   wspd = expand_input(wspd, state_size);
   psfc = expand_input(psfc, state_size);
   De = expand_input(De, state_size);
   ea_atm = expand_input(ea_atm, state_size);
   stability = expand_input(stability, state_size);
   roL = expand_input(roL, state_size);
   z0_bulk = expand_input(z0_bulk, state_size);

   diag_scalar = struct( ...
      'active', false(size(De)), ...
      'z_obs', NaN(size(De)), ...
      'u_star', NaN(size(De)), ...
      'Re', NaN(size(De)), ...
      'z0h', NaN(size(De)), ...
      'z0q', NaN(size(De)), ...
      'De_h', De, ...
      'De_e', De, ...
      'Qh', cv_air * (De .* stability .* (tair - T_sfc)), ...
      'Qe', roL .* De .* stability .* (epsilon ./ psfc .* (ea_atm - es_sfc)));

   valid = isfinite(wspd) & isfinite(De) & isfinite(z0_bulk) ...
      & wspd > 0 & De > 0 & z0_bulk > 0;
   if ~any(valid(:))
      return
   end

   Cd = NaN(size(De));
   Cd(valid) = De(valid) ./ wspd(valid);
   valid = valid & isfinite(Cd) & Cd > 0;
   if ~any(valid(:))
      return
   end

   log_z_obs_over_z0 = kappa ./ sqrt(Cd);
   z_obs = z0_bulk .* exp(log_z_obs_over_z0);
   u_star = sqrt(De .* wspd);
   Re = u_star .* z0_bulk / nu_air;
   [z0h, z0q] = scalar_roughness_lengths(z0_bulk, Re);

   valid = valid & isfinite(z_obs) & isfinite(u_star) & isfinite(Re) ...
      & isfinite(z0h) & isfinite(z0q) & z_obs > z0h & z_obs > z0q;
   if ~any(valid(:))
      return
   end

   De_h = NaN(size(De));
   De_e = NaN(size(De));
   De_h(valid) = wspd(valid) * kappa ^ 2 ./ ...
      (log(z_obs(valid) ./ z0_bulk(valid)) .* log(z_obs(valid) ./ z0h(valid)));
   De_e(valid) = wspd(valid) * kappa ^ 2 ./ ...
      (log(z_obs(valid) ./ z0_bulk(valid)) .* log(z_obs(valid) ./ z0q(valid)));

   diag_scalar.active = valid;
   diag_scalar.z_obs = z_obs;
   diag_scalar.u_star = u_star;
   diag_scalar.Re = Re;
   diag_scalar.z0h = z0h;
   diag_scalar.z0q = z0q;
   diag_scalar.De_h(valid) = De_h(valid);
   diag_scalar.De_e(valid) = De_e(valid);
   diag_scalar.Qh(valid) = cv_air * ...
      (De_h(valid) .* stability(valid) .* (tair(valid) - T_sfc(valid)));
   diag_scalar.Qe(valid) = roL(valid) .* De_e(valid) .* stability(valid) ...
      .* (epsilon ./ psfc(valid) .* (ea_atm(valid) - es_sfc(valid)));
end

function [z0h, z0q] = scalar_roughness_lengths(z0_bulk, Re)
   %SCALAR_ROUGHNESS_LENGTHS Apply the EXCOEFS_SCALAR roughness closures.

   z0h = NaN(size(Re));
   z0q = NaN(size(Re));

   low = Re < 0.135;
   mid = Re >= 0.135 & Re < 2.5;
   high = Re >= 2.5;

   z0_bulk_eval = z0_bulk + zeros(size(Re));
   z0h(low) = z0_bulk_eval(low) .* exp(1.250);
   z0q(low) = z0_bulk_eval(low) .* exp(1.610);

   if any(mid(:))
      log_Re = log(Re(mid));
      z0h(mid) = z0_bulk_eval(mid) .* exp(0.149 - 0.550 * log_Re);
      z0q(mid) = z0_bulk_eval(mid) .* exp(0.351 - 0.628 * log_Re);
   end

   if any(high(:))
      log_Re = log(Re(high));
      z0h(high) = z0_bulk_eval(high) .* exp( ...
         0.317 - 0.565 * log_Re - 0.183 * log_Re .^ 2);
      z0q(high) = z0_bulk_eval(high) .* exp( ...
         0.396 - 0.512 * log_Re - 0.180 * log_Re .^ 2);
   end
end

function x = expand_input(x, target_size)
   %EXPAND_INPUT Broadcast scalar diagnostics inputs over vector states.

   if isscalar(x)
      x = repmat(x, target_size);
   end
end
