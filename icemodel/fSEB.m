function f = fSEB(Ts, Ta, Qsi, Qli, albedo, wspd, ppt, tppt, Pa, De, ...
      ea, cv_air, cv_liq, emiss, SB, Tf, chi, roL, scoef, Qc, liqflag)
   %FSEB Surface energy balance equation
   %
   %#codegen

   persistent epsilon
   if isempty(epsilon)
      epsilon = icemodel.physicalConstant('epsilon');
   end

   f = chi * (1.0 - albedo) * Qsi + emiss * (Qli - SB * Ts ^ 4) ...
      + Qc + cv_air * De * (Ta - Ts) * STABLEFN(Ta, Ts, wspd, scoef) ...
      + roL * De * epsilon / Pa * (ea - VAPPRESS(Ts, liqflag)) ...
      * STABLEFN(Ta, Ts, wspd, scoef) ...
      + QADVECT(ppt, tppt, cv_liq);
end
