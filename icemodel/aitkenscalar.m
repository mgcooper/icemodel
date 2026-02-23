function [x_next, ok] = aitkenscalar( ...
      xkm2, xkm1, xk, x_fallback, jumpmax, varargin)
   %AITKENSCALAR Scalar Aitken Delta-squared acceleration with safeguards.
   %
   %#codegen

   % If not provided, assume use_aitken = true
   if nargin == 5
      use_aitken = true;
   else
      use_aitken = varargin{1};
   end

   x_next = x_fallback;
   if ~use_aitken
      ok = true;
      return
   end
   ok = false;

   % Need three Picard iterates before Aitken can be used.
   if isnan(xkm2) || isnan(xkm1)
      return
   end

   den = xk - 2.0 * xkm1 + xkm2;
   tol = 1e-12 * max([1.0, abs(xk), abs(xkm1), abs(xkm2)]);
   if abs(den) <= tol
      return
   end

   x_acc = xkm2 - (xkm1 - xkm2) ^ 2 / den;
   if ~isfinite(x_acc)
      return
   end

   % Reject extreme accelerated jumps and fallback to relaxed Picard.
   if abs(x_acc - x_fallback) > jumpmax
      return
   end

   x_next = x_acc;
   ok = true;
end
