function [z0h, z0q] = scalar_roughness_lengths(Re, z0m, use_snow_closure)
   %SCALAR_ROUGHNESS_LENGTHS Return scalar roughness lengths for bulk-MO.
   %
   % The bulk-MO scheme carries one momentum roughness length z0m and then
   % diagnoses scalar roughness lengths for heat and moisture:
   %   z0h = temperature roughness length [m]
   %   z0q = humidity roughness length [m]
   %   Re_* = u_* z0m / nu_air = roughness Reynolds number [1]
   %
   % The scalar roughness closure depends on surface type:
   %   Andreas (2002) over snow/firn
   %   Smeets and van den Broeke (2008) over rough bare ice
   %
   % The friction velocity u_* and the associated roughness Reynolds number
   % Re_* are diagnosed by the parent bulk-MO solver before this helper is
   % called. This helper owns only the scalar roughness closure itself.
   %
   %#codegen

   persistent z0_min z0_fallback
   persistent re1 re2 ch1 ch2 ch3 cq1 cq2 cq3 a0 a1 a2
   if isempty(z0_min)
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

   % When the roughness Reynolds number collapses to zero, turbulent exchange
   % also collapses. Return a tiny fallback scalar roughness so downstream
   % logarithms remain defined if the diagnostic state is inspected later.
   z0h = z0_fallback + zeros(size(Re));
   z0q = z0_fallback + zeros(size(Re));

   positive = real(Re) > 0;
   if ~any(positive(:))
      return
   end

   % The Andreas and Smeets closures are functions of the real-valued roughness
   % Reynolds number. Use the real part so complex-step perturbations preserve
   % the same physical closure branch.
   roughness_Re = max(real(Re), 1e-12);
   log_Re = log(roughness_Re);

   z0m_eval = z0m + zeros(size(Re));
   z0h_ice = z0m_eval .* exp(a0 + a1 * log_Re + a2 * log_Re .^ 2);
   z0q_ice = z0h_ice;

   z0h_snow = z0h_ice;
   z0q_snow = z0q_ice;
   low = roughness_Re <= re1;
   mid = roughness_Re > re1 & roughness_Re < re2;
   high = roughness_Re >= re2;
   z0h_snow(low) = z0m_eval(low) .* exp(ch1(1) + ch2(1) * log_Re(low) ...
      + ch3(1) * log_Re(low) .^ 2);
   z0q_snow(low) = z0m_eval(low) .* exp(cq1(1) + cq2(1) * log_Re(low) ...
      + cq3(1) * log_Re(low) .^ 2);
   z0h_snow(mid) = z0m_eval(mid) .* exp(ch1(2) + ch2(2) * log_Re(mid) ...
      + ch3(2) * log_Re(mid) .^ 2);
   z0q_snow(mid) = z0m_eval(mid) .* exp(cq1(2) + cq2(2) * log_Re(mid) ...
      + cq3(2) * log_Re(mid) .^ 2);
   z0h_snow(high) = z0m_eval(high) .* exp(ch1(3) + ch2(3) * log_Re(high) ...
      + ch3(3) * log_Re(high) .^ 2);
   z0q_snow(high) = z0m_eval(high) .* exp(cq1(3) + cq2(3) * log_Re(high) ...
      + cq3(3) * log_Re(high) .^ 2);

   if isscalar(use_snow_closure)
      snow_mask = positive & use_snow_closure;
   else
      snow_mask = positive & use_snow_closure;
   end
   ice_mask = positive & ~snow_mask;

   z0h(ice_mask) = z0h_ice(ice_mask);
   z0q(ice_mask) = z0q_ice(ice_mask);
   z0h(snow_mask) = z0h_snow(snow_mask);
   z0q(snow_mask) = z0q_snow(snow_mask);

   % Enforce a small positive lower bound so later logarithms remain defined.
   z0h = max(z0h, z0_min);
   z0q = max(z0q, z0_min);
end
