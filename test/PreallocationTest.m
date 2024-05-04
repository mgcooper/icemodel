classdef PreallocationTest < matlab.perftest.TestCase

   % This is mathworks example for class-based test of fast-executing code

   methods(Test)
      function testOnes(testCase)

         % standard approach (for "slow" code)
         % testCase.startMeasuring
         % x = ones(1,1e5);
         % testCase.stopMeasuring

         % "keepMeasuring" method for fast code
         while(testCase.keepMeasuring)
            x = ones(1,1e5);
         end

         testCase.verifyEqual(size(x),[1 1e5])
      end

      function testIndexingWithVariable(testCase)
         import matlab.unittest.constraints.IsSameSetAs

         % standard approach (for "slow" code)
         % testCase.startMeasuring
         % id = 1:1e5;
         % x(id) = 1;
         % testCase.stopMeasuring

         % "keepMeasuring" method for fast code
         while(testCase.keepMeasuring)
            id = 1:1e5;
            x(id) = 1;
         end

         testCase.verifyThat(x,IsSameSetAs(1))
      end

      function testIndexingOnLHS(testCase)
         import matlab.unittest.constraints.EveryElementOf
         import matlab.unittest.constraints.IsEqualTo

         % standard approach (for "slow" code)
         % testCase.startMeasuring
         % x(1:1e5) = 1;
         % testCase.stopMeasuring

         % "keepMeasuring" method for fast code
         while(testCase.keepMeasuring)
            x(1:1e5) = 1;
         end

         testCase.verifyThat(EveryElementOf(x),IsEqualTo(1))
      end

      function testForLoop(testCase)

         % standard approach (for "slow" code)
         % testCase.startMeasuring
         % for i=1:1e5
         %    x(i) = 1;
         % end
         % testCase.stopMeasuring

         % "keepMeasuring" method for fast code
         while(testCase.keepMeasuring)
            for i=1:1e5
               x(i) = 1;
            end
         end
         testCase.verifyNumElements(x,1e5)
      end
   end
end
