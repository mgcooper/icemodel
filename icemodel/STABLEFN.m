function S = STABLEFN(Ta, Ts, wspd, scoef)
   %STABLEFN Compute the stability function
   %
   % S = STABLEFN(Ta, Ts, wspd, scoef)
   %
   % Following Appendix 1 of Liston et al. 1999:
   %
   %    S = 1.0;
   %    Ri = 9.81 * z_obs / (Ta * wspd ^ 2) * (Ta - Ts);
   %    De = wcoef * wspd;
   %    gamma = 5.3 * 9.4 * De / wspd * sqrt((z_obs / z_0));
   %
   %    if Ri<0
   %       S = 1 - 9.4 * Ri / (1 + gamma * sqrt(abs(Ri)));
   %    elseif Ri > 0
   %       S = 1 / (1 + 4.7 * Ri) ^ 2;
   %    end
   %
   % For unstable, the absolute value of the Richardson number is implicit
   % because Ts>Ta implies Ts-Ta>0, see Appendix of Liston et al 1999
   %
   % Units: [1]
   %
   % See also: WINDCOEF
   %
   %#codegen

   % neutral blending width for neutral-transition behavior [K].
   dT0 = 0.05;

   if (Ts < Ta - dT0) % Stable case (Ri > 0).

      S = 1 / (1 + scoef(2) / (2 * wspd ^ 2) * (Ta - Ts) / Ta) ^ 2;

   elseif (Ts > Ta + dT0) % Unstable case (Ri < 0).

      S = 1 + scoef(2) / wspd ^ 2 * (Ts - Ta) / Ta ...
         / (1 + scoef(3) / wspd * sqrt((Ts - Ta) / Ta));

   else % Near-neutral transition, blend stable/unstable branches smoothly.
      Ts_lo = Ta - dT0;
      Ts_hi = Ta + dT0;
      S_stable = 1 / (1 + scoef(2) / (2 * wspd ^ 2) * (Ta - Ts_lo) / Ta) ^ 2;
      S_unstable = 1 + scoef(2) / wspd ^ 2 * (Ts_hi - Ta) / Ta ...
         / (1 + scoef(3) / wspd * sqrt((Ts_hi - Ta) / Ta));
      w = (Ts - Ts_lo) / (Ts_hi - Ts_lo);
      S = (1 - w) * S_stable + w * S_unstable;
   end
end
