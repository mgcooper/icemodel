function x = TRISOLVE(low, mid, upp, rhs)
   %TRISOLVE Solve tridiagonal matrix equation Ax = b for x.

   % Initialize
   x = zeros(size(rhs));
   n = numel(mid);

   % Forward elimination
   for k = 2:n
      mid(k) = mid(k) - low(k) / mid(k-1) * upp(k-1);
      rhs(k) = rhs(k) - low(k) / mid(k-1) * rhs(k-1);
   end

   % Back substitution
   x(n) = rhs(n) / mid(n);
   for k = n-1:-1:1
      x(k) = (rhs(k) - upp(k) * x(k+1)) / mid(k);
   end
end
