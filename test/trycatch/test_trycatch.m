classdef test_trycatch < matlab.perftest.TestCase

   properties
      a
      b
   end

   properties (TestParameter)
      
   end
   
%    methods(TestMethodSetup)
%       function setup(testCase)
% 
%          testCase.a = a;
%          testCase.b = b;
%       end
%    end

   methods(Test)
      function testTryCatch(testCase)

         while(testCase.keepMeasuring)
            call_trycatch(testCase.a, testCase.b);
         end
      end
   end
end
