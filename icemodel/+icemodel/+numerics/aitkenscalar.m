function [x_next, ok] = aitkenscalar(xk2, xk1, xk, x_fallback, jumpmax, varargin)
   %AITKENSCALAR Apply scalar Aitken Delta-squared acceleration with safeguards.
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
   if isnan(xk2) || isnan(xk1)
      return
   end

   den = xk - 2.0 * xk1 + xk2;
   tol = 1e-12 * max([1.0, abs(xk), abs(xk1), abs(xk2)]);
   if abs(den) <= tol
      return
   end

   x_acc = xk2 - (xk1 - xk2) ^ 2 / den;
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
