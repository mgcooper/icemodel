function balance = ENBAL(albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qm,emiss,chi)
   %ENBAL compute the energy balance
   balance = chi*Qsi*(1.0-albedo) + emiss*Qli + Qle + Qh + Qe + Qc - Qm;
end
