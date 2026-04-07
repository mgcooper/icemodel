function [stability, dstability] = stability_factor(T_sfc, tair, wspd, br_coefs)
   %STABILITY_FACTOR Compute the stability function and derivative wrt T_sfc.
   %
   %  stability = ...
   %     icemodel.surface.turbulence.bulk_richardson.stability_factor(...)
   %
   % Following the Louis/Liston bulk-Richardson parameterization
   % (see Appendix 1 of Liston et al. 1999):
   %
   %    S = 1.0;
   %    Ri = 9.81 * z_obs / (tair * wspd ^ 2) * (tair - T_sfc);
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
   % because T_sfc>tair implies T_sfc-tair>0, see Appendix of Liston et al 1999
   %
   % Legacy implementation, for comparison with notes.
   %
   %  if T_sfc < tair % Stable case
   %
   %     S = 1 / (1 + B1 / 2 * (tair - T_sfc)) ^ 2;
   %     dS = 2 * B1 / 2 / (1 + B1 / 2 * (tair - T_sfc)) ^ 3;
   %
   %  elseif T_sfc > tair % Unstable case
   %
   %     S = 1 + B1 * (T_sfc - tair) / (1 + B2 * sqrt(T_sfc - tair));
   %     dS = B1 / (1 + B2 * sqrt(T_sfc - tair)) - (B1 * B2 * (T_sfc - tair)) ...
   %        / (2 * (1 + B2 * sqrt(T_sfc - tair)) ^ 2 * sqrt(T_sfc - tair));
   %
   %  else % Neutrally stable case.
   %     S = 1.0;
   %     dS = 0;
   %  end
   %
   % Units: [1]
   %
   % See also: icemodel.surface.turbulence.bulk_richardson.exchange_coefficients
   %
   %#codegen

   % Neutral blending width for neutral-transition behavior [K].
   persistent dT0
   if isempty(dT0)
      dT0 = icemodel.parameterLookup( ...
         'thf_bulk_richardson_neutral_transition_width');
   end

   [B1, B2] = bulk_stability_terms(tair, wspd, br_coefs);

   % Branch on the real temperature only. This keeps the stable/unstable
   % switch fixed during complex-step differentiation while still letting the
   % downstream branch formulas operate on the possibly complex T_sfc value.
   T_sfc_real = real(T_sfc);

   if T_sfc_real < tair - dT0

      [stability, dstability] = stableBranch(T_sfc, tair, B1);

   elseif T_sfc_real > tair + dT0

      [stability, dstability] = unstableBranch(T_sfc, tair, B1, B2);

   else
      [stability, dstability] = neutralBranch(T_sfc, tair, B1, B2, dT0);
   end
end

function [B1, B2] = bulk_stability_terms(tair, wspd, br_coefs)
   %BULK_STABILITY_TERMS Return the bulk-Richardson branch coefficients.

   B1 = br_coefs(2) / (tair * wspd ^ 2);
   B2 = br_coefs(3) / (sqrt(tair) * wspd);
end

function [S, dSdT] = stableBranch(T_sfc, tair, B1)
   %STABLEBRANCH Evaluate the stable Louis/Liston branch and derivative.
   %
   % S = (1 + B1/2 * (tair - T_sfc)) ^ -2;
   % dSdT = B1 * S ^ 3/2;

   B3 = 1 + 0.5 * B1 * (tair - T_sfc);
   S = B3 ^ -2;

   if nargout > 1
      dSdT = B1 * B3 ^ -3;
   end
end

function [S, dSdT] = unstableBranch(T_sfc, tair, B1, B2)
   %UNSTABLEBRANCH Evaluate the unstable Louis/Liston branch and derivative.
   %
   % B3 = 1 + B2 * sqrt(T_sfc - tair)
   % S = 1 + B1 * (T_sfc - tair) / B3
   % dSdT = B1/B3 - (B1*B2 * (T_sfc - tair)) / (2 * B3^2 * sqrt(T_sfc - tair))

   dT = T_sfc - tair;
   sqrt_dT = sqrt(dT);
   B3 = 1 + B2 * sqrt_dT;

   S = 1 + B1 * dT / B3;

   if nargout > 1
      dSdT = B1 / B3 - 0.5 * B1 * B2 * sqrt_dT / B3^2;
   end
end

function [S, dSdT] = neutralBranch(T_sfc, tair, B1, B2, dT0)
   % Near-neutral transition, blend stable/unstable branches smoothly.

   T_sfc_lo = tair - dT0;
   T_sfc_hi = tair + dT0;

   S_stable = stableBranch(T_sfc_lo, tair, B1);
   S_unstable = unstableBranch(T_sfc_hi, tair, B1, B2);

   w = (T_sfc - T_sfc_lo) / (T_sfc_hi - T_sfc_lo);
   S = (1 - w) * S_stable + w * S_unstable;

   if nargout > 1
      dSdT = (S_unstable - S_stable) / (T_sfc_hi - T_sfc_lo);
   end
end
