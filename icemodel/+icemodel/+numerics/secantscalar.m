function [x_next, ok] = secantscalar(x_prev, r_prev, x, r, ...
      x_fallback, jumpmax, use_secant)
   %SECANTSCALAR Apply a safeguarded scalar secant step.
   %
   %#codegen

   % The caller's fallback is returned unless a finite pair brackets a fixed point.
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
