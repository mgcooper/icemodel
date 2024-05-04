function balance = ENBAL(albedo, emiss, chi, Qsi, Qli, Qle, Qh, Qe, Qc, Qa, Qm)
   %ENBAL Compute the energy balance
   %
   %#codegen
   balance = chi * Qsi * (1.0 - albedo) ...
      + emiss * Qli + Qle + Qh + Qe + Qc + Qa - Qm;
end
