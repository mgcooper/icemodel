function [Qm, Qf, Qh, Qe, Qc, Qle, balance, Tsfc, ea] = ENBALANCE( ...
      Tair, wspd, rh, Qsi, Qli, albedo, Pa, De, T, k_eff, Tf, dz, chi, Tsfc, ...
      cv_air, emiss, SB, roL, scoef, epsilon, fopts, liqflag, isicemodel)
   %ENBALANCE compute the surface energy-balance and solve for tsfc

   % Atmospheric vapor pressure from relative humidity data.
   ea = VAPPRESS(Tair, Tf, liqflag) * rh / 100;

   % Incoming longwave if not provided
   % Qli = LONGIN(ea, Tair, SB);

   % Compute the average station pressure.
   % Pa = PRESSURE(topo);

   % Compute the turbulent exchange coefficients.
   % De = EXCOEFS(wspd, wcoef);

   % Compute the flux contribution due to conduction.
   Qc = CONDUCT(k_eff, T, dz, Tsfc);

   % Solve the energy balance for the surface temperature.
   [Tsfc, ~] = SFCTEMP(Tair, Qsi, Qli, ea, albedo, De, Pa, wspd, ...
      cv_air, emiss, SB, Tf, Qc, Tsfc, chi, roL, scoef, fopts, liqflag);

   % Make the Tsfc_0 <= 0 C for surface flux calculations.
   % Let Tsfc remain > Tf for the upper boundary condition on ICEENBAL
   Tsfc0 = MELTTEMP(Tsfc, Tf);

   % Compute the stability function.
   S = STABLEFN(Tair, Tsfc0, wspd, scoef);

   % Compute saturation water vapor pressure at the surface.
   es = VAPPRESS(Tsfc0, Tf, liqflag);

   % Compute the latent heat flux.
   Qe = LATENT(De, S, ea, es, roL, epsilon, Pa);

   % Compute the sensible heat flux.
   Qh = SENSIBLE(De, S, Tair, Tsfc0, cv_air);

   % Compute the longwave flux emitted by the surface.
   Qle = LONGOUT(Tsfc0, emiss, SB);

   % Compute the energy flux available for melting or freezing.
   Qm = 0.0;
   Qf = 0.0;

   % if melting, compute melt energy
   if Tsfc >= Tf
      Qm = chi * Qsi * (1.0 - albedo) + emiss * Qli + Qle + Qh + Qe + Qc;
   else
      % else compute energy needed to reach melt temp
      Qf = -(chi * (1.0 - albedo) * Qsi ...
         + emiss * Qli ...
         + LONGOUT(Tf, emiss, SB) ...
         + LATENT(De, STABLEFN(Tair, Tf, wspd, scoef), ea, ...
            VAPPRESS(Tf, Tf, true), roL, epsilon, Pa) ...
         + SENSIBLE(De, STABLEFN(Tair, Tf, wspd, scoef), Tair, Tf, cv_air) ...
         + CONDUCT(k_eff, T, dz, Tf));
   end

   % Perform an energy balance check.
   balance = chi * Qsi * (1.0 - albedo) + emiss * Qli + Qle + Qh + Qe + Qc - Qm;
   
   % For a 'skin' surface energy balance model, reset Tsfc
   if ~isicemodel
      Tsfc = Tsfc0;
   end
end
