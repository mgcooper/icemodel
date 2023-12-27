function [Sc, Sp] = SFCFLIN(Ta, Qsi, Qli, albedo, wspd, Pa, De, ...
      ea, cv_air, emiss, SB, roL, scoef, chi, Tf, Ts, liqflag)
   %SFCFLIN Linearize the surface energy balance equation
   %
   % The linearization is of the form: F = Fc + Fp * T
   %
   % To construct the linearization, the SEB is separated into terms that are
   % independent of T (call them F0) and terms that depend on T. The latter
   % are linearized using truncated Taylor expansions. Terms which are
   % independent of the current T (i.e., those that depend on past values of
   % T or the derivative dF/dT) are grouped with the F0 terms to form Fc:
   %
   %  SEB = F0 + F(T)
   %      = F0 + F' + (dF/dT)' * (T - T')
   %      = F0 + F' + (dF/dT)' * T - (dF/dT)' * T'
   %      = F0 + F' - (dF/dT)' * T' + (dF/dT)' * T
   %      = Fc + Fp * T
   %
   % where
   %
   %  Fc = F0 + F' - (dF/dT)' * T'
   %  Fp = (dF/dT)'
   %
   % Note that the outgoing longwave and the saturation vapor pressure are
   % linearized around T_old but the stability function is not linearized,
   % thus T_old should be used to compute the stability function.
   %
   % All terms passed here are 'old', meaning the linearizations are computed
   % once at the start of the timestep, and on iterations updated as
   % F = Fc + Fp * T_new
   %
   % See also: SFCFLUX, SEBFLUX, SEBSOLVE

   if liqflag == true
      % Over water.
      A = 611.210;
      B = 17.502;
      C = 240.97;
   else
      % Over ice.
      A = 611.150;
      B = 22.452;
      C = 272.55;
   end

   % Compute the constants used in the stability coefficient computations
   B1 = scoef(2) / (Ta * wspd ^ 2);
   B2 = scoef(3) / (sqrt(Ta) * wspd);

   % Compute saturation vapor pressure
   es = A * exp(B * (Ts - Tf) / (C + Ts - Tf));

   % This accounts for an increase in turbulent fluxes under unstable conditions
   if Ts > Ta
      B3 = 1 + B2 * sqrt(Ts - Ta);
      S = 1 + B1 * (Ts - Ta) / B3;
   elseif Ts < Ta
      S = 1 / (1 + B1 / 2 * (Ta - Ts)) ^ 2;
   else
      S = 1.0;
   end

   % linearizations

   % outgoing longwave (note: -Qle = -emiss * SB * Ts ^ 4):
   Sc_Qle = 3 * emiss * SB * Ts ^ 4;
   Sp_Qle = -4 / 3 * Sc_Qle / Ts;

   % sensible heat flux:
   Sc_Qh = cv_air * De * S * Ta;
   Sp_Qh = -cv_air * De * S;

   % latent heat flux:
   Sc_Qe = roL * De * 0.622 / Pa * S ...
      * (ea - es * (1 - B * C * Ts / (C + Ts - Tf) ^ 2));
   Sp_Qe = -roL * De * 0.622 / Pa * S ...
      * es * B * C / (C + Ts - Tf) ^ 2;

   % combine net sw, incoming lw, conduction, and snow/rain heat flux:
   Sc = Sc_Qle + Sc_Qh + Sc_Qe + emiss * Qli + chi * Qsi * (1.0 - albedo);
   Sp = Sp_Qle + Sp_Qh + Sp_Qe;
end
