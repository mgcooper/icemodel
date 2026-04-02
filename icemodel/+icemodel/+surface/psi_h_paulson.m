function psi = psi_h_paulson(zeta)
   %PSI_H_PAULSON Dyer/Paulson unstable scalar profile correction.
   %
   % y = (1 - gamma zeta)^(1/2)
   % psi_h = ln(((1+y)/2)^2)

   %#codegen

   persistent gamma
   if isempty(gamma)
      gamma = icemodel.parameterLookup('thf_bulk_dyer_gamma');
   end

   y = (1 - gamma * zeta) ^ 0.5;
   psi = log(((1 + y) / 2) ^ 2);
end
