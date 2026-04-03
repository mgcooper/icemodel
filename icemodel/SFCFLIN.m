function [Sc, Sp] = SFCFLIN(tair, Qsi, Qli, albedo, wspd, ppt, tppt, psfc, De, ...
      ea_atm, roL, br_coefs, chi, T_sfc, liqflag)
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

   persistent cv_air cv_liq emiss SB epsilon
   if isempty(cv_air)
      [cv_air, cv_liq, SB, epsilon] = icemodel.physicalConstant( ...
         'cv_air', 'cv_liq', 'SB', 'epsilon');
      emiss = icemodel.parameterLookup('emiss');
   end

   % Saturation vapor pressure and derivative from VAPPRESS
   [es_sfc, des_dT] = VAPPRESS(T_sfc, liqflag);

   % Keep the linearization based on the bulk richardson scheme.
   stability = STABLEFN(T_sfc, tair, wspd, br_coefs);

   % linearizations

   % outgoing longwave (note: -Qle = -emiss * SB * T_sfc ^ 4):
   Sc_Qle = 3 * emiss * SB * T_sfc ^ 4;
   Sp_Qle = -4 / 3 * Sc_Qle / T_sfc;

   % sensible heat flux:
   Sc_Qh = cv_air * De * stability * tair;
   Sp_Qh = -cv_air * De * stability;

   % latent heat flux (linearization of es_sfc around T_sfc):
   % es(T) ≈ es_sfc + des_dT * (T - T_sfc) = (es_sfc - des_dT * T_sfc) + des_dT * T
   Sc_Qe = roL * De * epsilon / psfc * stability ...
      * (ea_atm - es_sfc + des_dT * T_sfc);
   Sp_Qe = -roL * De * epsilon / psfc * stability * des_dT;

   % precipitation-advected heat:
   Sc_Qa = QADVECT(ppt, tppt, cv_liq);
   Sp_Qa = 0;

   % combine net sw, incoming lw, and precipitation-advected heat:
   Sc = Sc_Qle + Sc_Qh + Sc_Qe + emiss * Qli + chi * Qsi * (1.0 - albedo) ...
      + Sc_Qa;
   Sp = Sp_Qle + Sp_Qh + Sp_Qe + Sp_Qa;
end
