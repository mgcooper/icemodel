function [Sc, Sp] = SFCFLIN(Ta, Qsi, Qli, albedo, wspd, ppt, tppt, Pa, De, ...
      ea, cv_air, cv_liq, emiss, SB, roL, scoef, chi, Tf, Ts, liqflag)
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
   %
   %#codegen

   persistent epsilon
   if isempty(epsilon)
      epsilon = icemodel.physicalConstant('epsilon');
   end

   % Saturation vapor pressure and derivative from VAPPRESS
   [es, des_dT] = VAPPRESS(Ts, liqflag);

   % Keep the linearization based on the bulk richardson scheme.
   S = STABLEFN(Ta, Ts, wspd, scoef);

   % linearizations

   % outgoing longwave (note: -Qle = -emiss * SB * Ts ^ 4):
   Sc_Qle = 3 * emiss * SB * Ts ^ 4;
   Sp_Qle = -4 / 3 * Sc_Qle / Ts;

   % sensible heat flux:
   Sc_Qh = cv_air * De * S * Ta;
   Sp_Qh = -cv_air * De * S;

   % latent heat flux (linearization of es around Ts):
   % es(T) ≈ es(Ts) + des_dT * (T - Ts) = (es - des_dT * Ts) + des_dT * T
   Sc_Qe = roL * De * epsilon / Pa * S ...
      * (ea - es + des_dT * Ts);
   Sp_Qe = -roL * De * epsilon / Pa * S * des_dT;

   % precipitation-advected heat:
   Sc_Qa = QADVECT(ppt, tppt, cv_liq);
   Sp_Qa = 0;

   % combine net sw, incoming lw, and precipitation-advected heat:
   Sc = Sc_Qle + Sc_Qh + Sc_Qe + emiss * Qli + chi * Qsi * (1.0 - albedo) ...
      + Sc_Qa;
   Sp = Sp_Qle + Sp_Qh + Sp_Qe + Sp_Qa;
end
