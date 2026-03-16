classdef TryCatchTest < matlab.perftest.TestCase
   %TRYCATCHTEST Microbenchmark for `catch ME` vs bare `catch`.
   %
   % This is not a core-model hotspot benchmark. Preserve it as a language
   % overhead experiment in case try/catch becomes relevant in future code.

   properties (TestParameter)
      handling_mode = {1, 2}
      throws_error = {false, true}
   end

   methods(Test)
      function testTryCatch(testCase, handling_mode, throws_error)
         if throws_error
            a = 2;
            b = "dog";
            batch_size = 1000;
         else
            a = 2;
            b = 1;
            batch_size = 10000;
         end

         if handling_mode == 1
            trycatch_impl = @trycatchWithHandling;
         else
            trycatch_impl = @trycatchWithoutHandling;
         end

         % Preflight the selected branch before measuring it repeatedly.
         trycatch_impl(a, b);

         while testCase.keepMeasuring
            for n = 1:batch_size
               trycatch_impl(a, b);
            end
         end
      end
   end
end

function c = trycatchWithHandling(a, b)

   try
      c = a - b;
   catch ME
      c = [];
   end
end

function c = trycatchWithoutHandling(a, b)

   try
      c = a - b;
   catch
      c = [];
   end
end
