function [Qm, Qf, Qh, Qe, Qc, Qle, balance, Ts, ea, ok] = ENBALANCE( ...
      Ts, Ta, Qsi, Qli, albedo, wspd, rh, ppt, tppt, Pa, De, T, k_eff, dz, ...
      chi, roL, scoef, liqflag, Tflag, solver, ro_sfc, snow_depth, opts)
   %ENBALANCE compute the surface energy-balance and solve for tsfc
   %
   %#codegen

   % Load physical constants and parameters
   persistent Tf emiss
   if isempty(Tf)
      Tf = icemodel.physicalConstant("Tf");
      emiss = icemodel.parameterLookup('emiss');
   end

   % Atmospheric vapor pressure from relative humidity data.
   ea = icemodel.surface.atmospheric_vapor_pressure(Ta, rh, liqflag);

   % Incoming longwave if not provided
   % Qli = LONGIN(Ta, ea, SB);

   % Compute the average station pressure.
   % Pa = PRESSURE(topo);

   % Compute the turbulent exchange coefficients.
   % De = EXCOEFS(wspd, wcoef);

   % Solve the energy balance for the surface temperature.
   [Ts, ok] = SEBSOLVE(Ts, Ta, Qsi, Qli, albedo, wspd, ppt, tppt, Pa, De, ...
      ea, chi, roL, scoef, liqflag, T, k_eff, dz, ro_sfc, snow_depth, ...
      solver, false, opts);

   % Ensure Ts_0 <= 0 C for surface flux calculations.
   % Let Ts remain > Tf for the upper boundary condition on ICEENBAL'
   Ts0 = MELTTEMP(Ts, Tf);

   % Gather the terms in the surface energy balance.
   [~, Qe, Qh, Qc, Qa, Qle] = ...
      icemodel.surface.surface_energy_balance_terms(Ts0, Ta, Qsi, ...
      Qli, albedo, wspd, ppt, tppt, Pa, De, ea, T, k_eff, dz, roL, chi, ...
      scoef, liqflag, ro_sfc, snow_depth, opts);

   % Compute the energy flux available for melting or freezing.
   [Qm, Qf] = MFENERGY(albedo, Qsi, Qli, Qle, Qh, Qe, Qc, Qa, Ts, ...
      Tf, Ta, wspd, ppt, tppt, De, ea, roL, Pa, k_eff, T, dz, scoef, chi, ...
      ro_sfc, snow_depth, opts);

   % Perform an energy balance check. When Ts == Tf, balance = 0.0.
   balance = ENBAL(albedo, emiss, chi, Qsi, Qli, Qle, Qh, Qe, Qc, Qa, Qm);

   % For a 'skin' surface energy balance model, reset Ts
   if Tflag
      Ts = Ts0;
   end
end
