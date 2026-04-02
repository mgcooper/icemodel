function [z0h, z0q, u_star, Re] = scalar_roughness_bulk_mo(wspd, z0m, ...
      psi_m0, psi_mz, nu_air, z_wind, use_snow_closure)
   %SCALAR_ROUGHNESS_BULK_MO Return scalar roughness lengths for bulk-MO.
   %
   % The bulk-MO scheme carries one momentum roughness length z0m and then
   % diagnoses scalar roughness lengths for heat and moisture:
   %   z0h = temperature roughness length [m]
   %   z0q = humidity roughness length [m]
   %   u_* = friction velocity [m s^-1]
   %   Re_* = u_* z0m / nu_air = roughness Reynolds number [1]
   %
   % The scalar roughness closure depends on surface type:
   %   Andreas (2002) over snow/firn
   %   Smeets and van den Broeke (2008) over rough bare ice
   %
   % Inputs psi_m0 and psi_mz are the Monin-Obukhov momentum profile
   % corrections evaluated at z0m/L and z_wind/L. They enter through the
   % momentum transfer denominator used to diagnose u_*.
   %
   %#codegen

   persistent kappa z0_min z0_fallback
   persistent re1 re2 ch1 ch2 ch3 cq1 cq2 cq3 a0 a1 a2
   if isempty(kappa)
      kappa = icemodel.physicalConstant('kappa');
      [z0_min, z0_fallback, re1, re2, ch1, ch2, ch3, cq1, cq2, cq3, ...
         a0, a1, a2] = icemodel.parameterLookup( ...
         'thf_bulk_scalar_roughness_min', ...
         'thf_bulk_scalar_roughness_fallback', ...
         'thf_bulk_andreas_re1', 'thf_bulk_andreas_re2', ...
         'thf_bulk_andreas_ch1', 'thf_bulk_andreas_ch2', ...
         'thf_bulk_andreas_ch3', 'thf_bulk_andreas_cq1', ...
         'thf_bulk_andreas_cq2', 'thf_bulk_andreas_cq3', ...
         'thf_bulk_smeets_a0', 'thf_bulk_smeets_a1', ...
         'thf_bulk_smeets_a2');
   end

   % A zero-wind limit has no turbulent exchange. Return a tiny fallback
   % scalar roughness so downstream logs stay defined if the diagnostic state
   % is inspected later.
   if wspd <= 0
      z0h = z0_fallback;
      z0q = z0_fallback;
      u_star = 0;
      Re = 0;
      return
   end

   % Convert the observed wind speed to friction velocity using the current
   % momentum roughness and stability corrections.
   denom = log(z_wind / z0m) - psi_mz + psi_m0;
   if abs(denom) < 1e-12
      denom = denom + sign_or_one(denom) * 1e-12;
   end

   % Diagnose friction velocity and the associated roughness Reynolds number.
   % u_* = kappa U / (ln(z/z0m) - psi_m(z/L) + psi_m(z0m/L))
   u_star = kappa * wspd / denom;
   Re = u_star * z0m / nu_air;

   % The Andreas and Smeets closures are functions of the real-valued
   % roughness Reynolds number. Use the real part so complex-step
   % perturbations preserve the same physical closure branch.
   roughness_re = max(real(Re), 1e-12);
   log_Re = log(roughness_re);

   % Snow/firn uses the Andreas piecewise polynomial fits across three
   % roughness-Reynolds regimes.
   if use_snow_closure
      if roughness_re <= re1
         idx = 1;
      elseif roughness_re < re2
         idx = 2;
      else
         idx = 3;
      end

      % Andreas returns distinct scalar roughness lengths for heat and
      % moisture, both scaled from the momentum roughness.
      z0h = z0m * exp(ch1(idx) + ch2(idx) * log_Re + ch3(idx) * log_Re ^ 2);
      z0q = z0m * exp(cq1(idx) + cq2(idx) * log_Re + cq3(idx) * log_Re ^ 2);
   else
      % Rough bare ice uses the Smeets and van den Broeke scalar closure.
      z0h = z0m * exp(a0 + a1 * log_Re + a2 * log_Re ^ 2);
      z0q = z0h;
   end

   % Enforce a small positive lower bound so later logarithms remain defined.
   z0h = max(z0h, z0_min);
   z0q = max(z0q, z0_min);
end
