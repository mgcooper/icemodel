classdef InplaceTest < matlab.perftest.TestCase

   % Test speed impact of in-place operations
   properties
      A
      B
   end
   methods(TestMethodSetup)
      function generateTestData(testCase)
         testCase.A = rand(1000000, 1);
         testCase.B = rand(1000000, 1);
      end
   end

   methods(Test)
      function test_inplace(testCase)

         A = testCase.A;
         B = testCase.B;
         
         while(testCase.keepMeasuring)
            A(:) = A(:) + B(:);
         end
      end
      
      function test_notinplace(testCase)

         A = testCase.A;
         B = testCase.B;
         
         while(testCase.keepMeasuring)
            A = A + B;
         end
      end

      % function testForLoop(testCase)
      % 
      %    while(testCase.keepMeasuring)
      %       for i=1:1e5
      %          x(i) = 1;
      %       end
      %    end
      %    testCase.verifyNumElements(x,1e5)
      % end
   end
end
