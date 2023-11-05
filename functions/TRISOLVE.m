function x = TRISOLVE(lower, middle, upper, rhs)
   %TRISOLVE Solve tridiagonal matrix

   % Initialize
   x = zeros(size(rhs));
   n = numel(middle);

   % Forward elimination
   for k = 2:n
      xfactor     =  lower(k)  / middle(k-1);
      middle(k)   =  middle(k) - xfactor * upper(k-1);
      rhs(k)      =  rhs(k)    - xfactor * rhs(k-1);
   end

   % Back substitution
   x(n) = rhs(n) / middle(n);
   for k = n-1:-1:1
      x(k) = (rhs(k) - upper(k) * x(k+1)) / middle(k);
   end
end
