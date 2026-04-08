function balance = ENBAL(chi, albedo, Qsi, Qli, Qle, Qh, Qe, Qc, Qa, Qm)
   %ENBAL Compute the energy balance from known energy fluxes
   %
   %#codegen
   persistent emiss
   if isempty(emiss)
      emiss = icemodel.parameterLookup('emiss');
   end
   balance = chi * Qsi * (1.0 - albedo) ...
      + emiss * Qli + Qle + Qh + Qe + Qc + Qa - Qm;


   % what surface_energy_balance_residual does:
   % [Qe, Qh] = icemodel.surface.diagnose_turbulent_heat_fluxes(...)
   % Qle = - emiss * SB * T_sfc ^ 4;
   % Qa = QADVECT(ppt, tppt, cv_liq);
   %
   % residual = chi * Qsi * (1.0 - albedo) ...
   %    + emiss * Qli + Qle + Qh + Qe + Qc + Qa;

   % so this function assumes all fluxes are known, whereas
   % surface_energy_balance_residual does not.
end
