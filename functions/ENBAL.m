function balance = ENBAL(albedo, emiss, chi, Qsi, Qli, Qle, Qh, Qe, Qc, Qm)
   %ENBAL Compute the energy balance
   balance = chi * Qsi * (1.0 - albedo) + emiss * Qli + Qle + Qh + Qe + Qc - Qm;
end
