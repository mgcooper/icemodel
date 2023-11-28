function [Qm, Qf] = MFENERGY(albedo, Qsi, Qli, Qle, Qh, Qe, Qc, Tsfc, ...
      Tf, Tair, wspd, De, ea, roL, Pa, cv_air, emiss, SB, k_eff, T, dz, ...
      epsilon, scoef,chi)
   %MFENERGY Compute the energy for melting/freezing
   %
   % If Qm is > 0, this is the energy available for melting.
   %   If Qm is < 0, this is the energy available for freezing
   %   For the case of melting conditions, the fluxes that are passed into
   %   this function are computed using Tsfc0, which is set to Tf, meaning the
   %   residual Qm is correct. For the case of freezing, Tsfc is unaltered, so
   %   it can be passed into the Qf block below to get the correct fluxes.
   %
   % See also:

   Qm = 0.0;
   Qf = 0.0;
   
   if Tsfc >= Tf
      % Compute melt energy
      Qm = chi * Qsi * (1.0 - albedo) + emiss * Qli + Qle + Qh + Qe + Qc;
      
   else
      % Compute energy needed to reach melt temp (freeze energy)
      Qf = -(chi * (1.0 - albedo) * Qsi ...
         + emiss * Qli + LONGOUT(Tf, emiss, SB) ...
         + LATENT(De, STABLEFN(Tair, Tf, wspd, scoef), ea, ...
         VAPPRESS(Tf, Tf, true), roL, epsilon, Pa) ...
         + SENSIBLE(De, STABLEFN(Tair, Tf, wspd, scoef), Tair, Tf, cv_air) ...
         + CONDUCT(k_eff, T, dz, Tf));

      % xS = STABLEFN(Tair,Tf,wspd,scoef);
      % xes = VAPPRESS(Tf,Tf,true); % assume melting
      % xQe = LATENT(De,xS,ea,xes,roL,epsilon,Pa);
      % xQh = SENSIBLE(De,xS,Tair,Tf,cv_air);
      % xQle = LONGOUT(Tf,emiss,SB);
      % xQc = CONDUCT(k_eff,T,dz,Tf);
      % xQn = chi*Qsi*(1.0-albedo) + emiss*Qli;
      % Qf = -(xQn + xQle + xQh + xQe + xQc);
   end
end
