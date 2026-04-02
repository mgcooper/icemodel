function [S, dSdTs] = STABLEFN(Ta, Ts, wspd, scoef)
   %STABLEFN Compute the stability function
   %
   % S = STABLEFN(Ta, Ts, wspd, scoef)
   %
   % Following the Louis/Liston bulk-Richardson parameterization
   % (see Appendix 1 of Liston et al. 1999):
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
   % Legacy implementation, for comparison with notes.
   %
   %  if Ts < Ta % Stable case
   %
   %     S = 1 / (1 + B1 / 2 * (Ta - Ts)) ^ 2;
   %     dS = 2 * B1 / 2 / (1 + B1 / 2 * (Ta - Ts)) ^ 3;
   %
   %  elseif Ts > Ta % Unstable case
   %
   %     S = 1 + B1 * (Ts - Ta) / (1 + B2 * sqrt(Ts - Ta));
   %     dS = B1 / (1 + B2 * sqrt(Ts - Ta)) - (B1 * B2 * (Ts - Ta)) ...
   %        / (2 * (1 + B2 * sqrt(Ts - Ta)) ^ 2 * sqrt(Ts - Ta));
   %
   %  else % Neutrally stable case.
   %     S = 1.0;
   %     dS = 0;
   %  end
   %
   % Units: [1]
   %
   % See also: WINDCOEF
   %
   %#codegen

   % Neutral blending width for neutral-transition behavior [K].
   persistent dT0
   if isempty(dT0)
      dT0 = icemodel.parameterLookup( ...
         'thf_bulk_richardson_neutral_transition_width');
   end

   [B1, B2] = bulk_stability_terms(Ta, wspd, scoef);

   % Branch on the real temperature only. This keeps the stable/unstable
   % switch fixed during complex-step differentiation while still letting the
   % downstream branch formulas operate on the possibly complex Ts value.
   Ts_real = real(Ts);

   if Ts_real < Ta - dT0

      [S, dSdTs] = stableBranch(Ta, Ts, B1);

   elseif Ts_real > Ta + dT0

      [S, dSdTs] = unstableBranch(Ta, Ts, B1, B2);

   else
      [S, dSdTs] = neutralBranch(Ta, Ts, B1, B2, dT0);
   end
end

function [B1, B2] = bulk_stability_terms(Ta, wspd, scoef)
   %BULK_STABILITY_TERMS Return the bulk-Richardson branch coefficients.

   B1 = scoef(2) / (Ta * wspd ^ 2);
   B2 = scoef(3) / (sqrt(Ta) * wspd);
end

function [S, dSdTs] = stableBranch(Ta, Ts, B1)
   %STABLEBRANCH Evaluate the stable Louis/Liston branch and derivative.
   %
   % S = (1 + B1/2 * (Ta - Ts)) ^ -2;
   % dSdTs = B1 * S ^ 3/2;

   B3 = 1 + 0.5 * B1 * (Ta - Ts);
   S = B3 ^ -2;

   if nargout > 1
      dSdTs = B1 * B3 ^ -3;
   end
end

function [S, dSdTs] = unstableBranch(Ta, Ts, B1, B2)
   %UNSTABLEBRANCH Evaluate the unstable Louis/Liston branch and derivative.
   %
   % B3 = 1 + B2 * sqrt(Ts - Ta)
   % S = 1 + B1 * (Ts - Ta) / B3
   % dSdTs = B1 / B3 - (B1 * B2 * (Ts - Ta)) / (2 * B3 ^ 2 * sqrt(Ts - Ta))

   dT = Ts - Ta;
   sqrt_dT = sqrt(dT);
   B3 = 1 + B2 * sqrt_dT;

   S = 1 + B1 * dT / B3;

   if nargout > 1
      dSdTs = B1 / B3 - 0.5 * B1 * B2 * sqrt_dT / B3^2;
   end
end

function [S, dSdTs] = neutralBranch(Ta, Ts, B1, B2, dT0)
   % Near-neutral transition, blend stable/unstable branches smoothly.

   Ts_lo = Ta - dT0;
   Ts_hi = Ta + dT0;

   S_stable = stableBranch(Ta, Ts_lo, B1);
   S_unstable = unstableBranch(Ta, Ts_hi, B1, B2);

   w = (Ts - Ts_lo) / (Ts_hi - Ts_lo);
   S = (1 - w) * S_stable + w * S_unstable;

   if nargout > 1
      dSdTs = (S_unstable - S_stable) / (Ts_hi - Ts_lo);
   end
end
