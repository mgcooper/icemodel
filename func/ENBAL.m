%--------------------------------------------------------------------------
%   PERFORM AN ENERGY BALANCE CHECK
%--------------------------------------------------------------------------
function balance = ENBAL(albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qm,emiss,chi)
   balance = chi*Qsi*(1.0-albedo) + emiss*Qli + Qle + Qh + Qe + Qc - Qm;