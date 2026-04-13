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
   if real(Re) <= 0
      z0h = z0_fallback;
      z0q = z0_fallback;
      return
   end

   % The Andreas and Smeets closures are functions of the real-valued roughness
   % Reynolds number. Use the real part so complex-step perturbations preserve
   % the same physical closure branch.
   roughness_Re = max(real(Re), 1e-12);
   log_Re = log(roughness_Re);

   % Snow/firn uses the Andreas piecewise polynomial fits across three
   % roughness-Reynolds regimes.
   if use_snow_closure
      if roughness_Re <= re1
         idx = 1;
      elseif roughness_Re < re2
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
