function [Qm, Qf, Qh, Qe, Qc, Qle, balance, T_sfc, ea_atm, ok] = ENBALANCE( ...
      T_sfc, tair, Qsi, Qli, albedo, wspd, rh, ppt, tppt, psfc, De, T, k_eff, dz, ...
      chi, roL, br_coefs, liqflag, Tflag, solver, ro_sfc, snow_depth, opts)
   %ENBALANCE compute the surface energy-balance and solve for tsfc
   %
   %#codegen

   % Atmospheric vapor pressure from relative humidity data.
   ea_atm = icemodel.surface.atmospheric_vapor_pressure(tair, rh, liqflag);

   % Incoming longwave if not provided
   % Qli = LONGIN(tair, ea_atm, SB);

   % Compute the average station pressure.
   % psfc = PRESSURE(topo);

   % Compute the turbulent exchange coefficients.
   % De = EXCOEFS(wspd, wcoef);

   % Solve the energy balance for the surface temperature.
   [T_sfc, ok] = SEBSOLVE(T_sfc, tair, Qsi, Qli, albedo, wspd, ppt, tppt, psfc, De, ...
      ea_atm, chi, roL, br_coefs, liqflag, T, k_eff, dz, ro_sfc, snow_depth, ...
      solver, false, opts);

   % Ensure T_sfc_phys <= Tf for surface flux calculations.
   % Let T_sfc remain > Tf for the upper boundary condition on ICEENBAL.
   T_sfc_phys = icemodel.kernels.physical_surface_temperature(T_sfc);

   % Gather the terms in the surface energy balance.
   [~, Qe, Qh, Qc, Qa, Qle] = ...
      icemodel.surface.surface_energy_balance_terms(T_sfc_phys, tair, Qsi, ...
      Qli, albedo, wspd, ppt, tppt, psfc, De, ea_atm, T, k_eff, dz, roL, chi, ...
      br_coefs, liqflag, ro_sfc, snow_depth, opts);

   % Compute the energy flux available for melting or freezing.
   [Qm, Qf] = MFENERGY(T_sfc_phys, chi, albedo, Qsi, Qli, Qle, Qh, Qe, Qc, Qa, ...
      tair, wspd, psfc, ppt, tppt, ea_atm, br_coefs, De, T, k_eff, dz, ...
      ro_sfc, snow_depth, roL, opts);

   % Perform an energy balance check. When T_sfc == Tf, balance = 0.0.
   balance = ENBAL(chi, albedo, Qsi, Qli, Qle, Qh, Qe, Qc, Qa, Qm);

   % For a 'skin' surface energy balance model, reset T_sfc
   if Tflag
      T_sfc = T_sfc_phys;
   end
end
