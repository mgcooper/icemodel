function psi = psi_m_paulson(zeta)
   %PSI_M_PAULSON Paulson unstable momentum profile correction.
   %
   % x = (1 - gamma zeta)^(1/4)
   % psi_m = ln(((1+x)/2)^2 (1+x^2)/2) - 2 atan(x) + pi/2

   %#codegen

   persistent gamma
   if isempty(gamma)
      gamma = icemodel.parameterLookup('thf_bulk_dyer_gamma');
   end

   x = (1 - gamma * zeta) .^ 0.25;
   psi = log(((1 + x) ./ 2) .^ 2 .* (1 + x .^ 2) ./ 2) - 2 * atan(x) + pi / 2;
end
