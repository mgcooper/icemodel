function [Qm, Qf] = MFENERGY(albedo, Qsi, Qli, Qle, Qh, Qe, Qc, Qa, Ts, ...
      Tf, Ta, wspd, ppt, tppt, De, ea, roL, Pa, cv_liq, cv_air, emiss, SB, ...
      k_eff, T, dz, epsilon, scoef, chi)
   %MFENERGY Compute the energy for melting/freezing
   %
   %  Qm > 0 = the energy surplus available for melting ice.
   %  Qm < 0 = the energy deficit required to refreeze liquid water.
   %
   %  The fluxes passed into this function must be computed using min(Ts, Tf).
   %  Thus for melting (Ts >= Tf), the fluxes are computed using Ts = Tf, and
   %  the ENBAL residual is Qm. For freezing (Ts < Tf), the fluxes passed into
   %  this function are computed using Ts < Tf, and the energy deficit Qf is
   %  computed below using Ts = Tf. Noticing that the fluxes used here are
   %  computed with Ts = Tf in both cases, it is tempting to do so here or
   %  outside this function, but remember that when Ts < Tf, the fluxes must be
   %  computed using Ts unaltered, so it is better to compute them once,
   %  correctly, outside this function, then re-compute them, if needed, in the
   %  freezing case below.
   %
   % See also:
   %
   %#codegen

   Qm = 0.0;
   Qf = 0.0;

   if Ts >= Tf
      % Compute melt energy
      Qm = ENBAL(albedo, emiss, chi, Qsi, Qli, Qle, Qh, Qe, Qc, Qa, 0.0);
   else
      % Compute energy needed to reach melt temp (energy deficit)
      Qf = -(chi * (1.0 - albedo) * Qsi ...
         + emiss * Qli + LONGOUT(Tf, emiss, SB) ...
         + LATENT(De, STABLEFN(Ta, Tf, wspd, scoef), ea, ...
         VAPPRESS(Tf, Tf, true), roL, epsilon, Pa) ...
         + SENSIBLE(De, STABLEFN(Ta, Tf, wspd, scoef), Ta, Tf, cv_air) ...
         + QADVECT(ppt, tppt, cv_liq) ...
         + CONDUCT(k_eff, T, dz, Tf) ...
         );
   end
end
