function balance = ENBAL(chi, albedo, Qsi, Qli, Qle, Qh, Qe, Qc, Qa, Qm)
   %ENBAL Compute the energy balance
   %
   %#codegen
   persistent emiss
   if isempty(emiss)
      emiss = icemodel.parameterLookup('emiss');
   end
   balance = chi * Qsi * (1.0 - albedo) ...
      + emiss * Qli + Qle + Qh + Qe + Qc + Qa - Qm;
end
