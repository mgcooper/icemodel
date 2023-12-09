function [Qm, Qf, Qh, Qe, Qc, Qle, balance, Ts, ea] = ENBALANCE( ...
      Ta, Qsi, Qli, albedo, wspd, rh, Pa, De, T, k_eff, Tf, dz, chi, Ts, ...
      cv_air, emiss, SB, roL, scoef, epsilon, liqflag, isicemodel)
   %ENBALANCE compute the surface energy-balance and solve for tsfc

   % Atmospheric vapor pressure from relative humidity data.
   ea = VAPPRESS(Ta, Tf, liqflag) * rh / 100;

   % Incoming longwave if not provided
   % Qli = LONGIN(Ta, ea, SB);

   % Compute the average station pressure.
   % Pa = PRESSURE(topo);

   % Compute the turbulent exchange coefficients.
   % De = EXCOEFS(wspd, wcoef);

   % Compute the flux contribution due to conduction.
   Qc = CONDUCT(k_eff, T, dz, Ts);

   % Solve the energy balance for the surface temperature.
   [Ts, ~] = SFCTEMP(Ta, Qsi, Qli, albedo, wspd, Pa, De, ...
      ea, cv_air, emiss, SB, Tf, chi, roL, scoef, Qc, liqflag);
   
   % Make the Ts_0 <= 0 C for surface flux calculations.
   % Let Ts remain > Tf for the upper boundary condition on ICEENBAL
   Ts0 = MELTTEMP(Ts, Tf);

   % Compute the stability function.
   S = STABLEFN(Ta, Ts0, wspd, scoef);

   % Compute saturation water vapor pressure at the surface.
   es = VAPPRESS(Ts0, Tf, liqflag);

   % Compute the latent heat flux.
   Qe = LATENT(De, S, ea, es, roL, epsilon, Pa);

   % Compute the sensible heat flux.
   Qh = SENSIBLE(De, S, Ta, Ts0, cv_air);

   % Compute the longwave flux emitted by the surface.
   Qle = LONGOUT(Ts0, emiss, SB);

   % Compute the energy flux available for melting or freezing.
   [Qm, Qf] = MFENERGY(albedo, Qsi, Qli, Qle, Qh, Qe, Qc, Ts, ...
      Tf, Ta, wspd, De, ea, roL, Pa, cv_air, emiss, SB, k_eff, T, dz, ...
      epsilon, scoef, chi, ctype);

   % Perform an energy balance check. When Ts == Tf, balance = 0.0.
   balance = ENBAL(albedo, emiss, chi, Qsi, Qli, Qle, Qh, Qe, Qc, Qm);

   % For a 'skin' surface energy balance model, reset Ts
   if Tflag
      Ts = Ts0;
   end
end
