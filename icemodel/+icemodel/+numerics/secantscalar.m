function [x_next, ok] = secantscalar(x_prev, r_prev, x, r, ...
      x_fallback, jumpmax, use_secant)
   %SECANTSCALAR Apply a safeguarded scalar secant step.
   %
   %  [x_next, ok] = icemodel.numerics.secantscalar(x_prev, r_prev, x, r, ...
   %     x_fallback, jumpmax, use_secant)
   %
   % Acts only when the last two residuals bracket a root, meaning both are
   % finite and nonzero with opposite signs. Otherwise the caller's fallback
   % is returned unchanged, so early iterations proceed on the fallback alone.
   % A step further than jumpmax from x is clamped to x +/- jumpmax. The
   % fallback is returned only if the clamped step leaves the bracket.
   %
   % Inputs
   %  x_prev      - iterate from the previous call
   %  r_prev      - residual at x_prev
   %  x           - current iterate
   %  r           - residual at x
   %  x_fallback  - value to return when no safeguarded step is available
   %  jumpmax     - largest accepted move away from x
   %  use_secant  - false returns x_fallback without attempting a step
   %
   % Outputs
   %  x_next - the secant or bisection step, or x_fallback
   %  ok     - true when a safeguarded step was taken
   %
   % See also: icemodel.numerics.aitkenscalar
   %
   %#codegen

   % Return the caller's fallback unless a finite pair brackets a root.
   x_next = x_fallback;
   ok = false;
   if ~use_secant || ~isfinite(x_prev) || ~isfinite(r_prev) || ...
         ~isfinite(x) || ~isfinite(r) || r_prev == 0 || r == 0 || ...
         sign(r_prev) == sign(r)
      return
   end

   % Use the bracket midpoint whenever interpolation is ill-conditioned or
   % lands indistinguishably close to an endpoint.
   x_lo = min(x_prev, x);
   x_hi = max(x_prev, x);
   x_mid = 0.5 * x_lo + 0.5 * x_hi;
   edge_tol = 1e-12 * max([1.0, abs(x_lo), abs(x_hi)]);
   den = r - r_prev;
   den_tol = 1e-12 * max([1.0, abs(r), abs(r_prev)]);
   if abs(den) <= den_tol
      x_sec = x_mid;
   else
      x_sec = x - r * (x - x_prev) / den;
      if ~isfinite(x_sec) || x_sec <= x_lo + edge_tol || ...
            x_sec >= x_hi - edge_tol
         x_sec = x_mid;
      end
   end

   % Bound the bracketed step from the current evaluated endpoint. The
   % caller's fallback can lie outside the bracket, so it cannot define the
   % safety distance once a sign-changing residual pair is available.
   d_sec = x_sec - x;
   if abs(d_sec) > jumpmax
      x_sec = x + sign(d_sec) * jumpmax;
   end

   % Retain the caller's fallback only if the bounded step is unusable.
   if x_sec <= x_lo || x_sec >= x_hi
      return
   end

   x_next = x_sec;
   ok = true;
end
