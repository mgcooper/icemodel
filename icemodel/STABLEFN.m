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

   if (Ts < Ta) % Stable case (Ri > 0).

      S = 1 / (1 + scoef(2) / (2 * wspd ^ 2) * (Ta - Ts) / Ta) ^ 2;

   elseif (Ts > Ta) % Unstable case (Ri < 0).

      S = 1 + scoef(2) / wspd ^ 2 * (Ts - Ta) / Ta ...
         / (1 + scoef(3) / wspd * sqrt((Ts - Ta) / Ta));

   else % Neutrally stable case. (Ts == Ta)
      S = 1.0;
   end
end
