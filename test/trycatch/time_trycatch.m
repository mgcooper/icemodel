function t = time_trycatch(a, b, option)
   
   if nargin == 0 || isempty(a) && isempty(b)
      a(1) = 2;
      a(2) = 2;
      b(1) = 1;
      b(2) = "dog";
      option(1) = 1;
      option(2) = 2;
   elseif nargin == 2
      option = 1;
      option(1) = 1;
      option(2) = 2;
   end

   t = zeros(2, 2);
   
   for n = 1:2
      opt = option(n);
      for m = 1:2
         f = @() call_trycatch(a(m), b(m), opt);
         t(n, m) = timeit(f);
         
         if opt == 1
            fprintf('%.10f seconds: Time for trycatch(%d, %d) with handling\n', ...
               t(n, m), a(m), b(m));
         elseif opt == 2
            fprintf('%.10f seconds: Time for trycatch(%d, %d) without handling\n', ...
               t(n, m), a(m), b(m));
         end
      end
   end
end
