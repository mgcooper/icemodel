classdef IcemodelPerfTest < matlab.perftest.TestCase

   properties
      opts
      savedata = false;
      sitename = 'behar';
      forcings = 'kanm';
      userdata = 'modis';
      uservars = 'albedo';
      simmodel = 'icemodel';
      simyears = 2016;
   end

   properties (TestParameter)
      solver = {1, 2, 3}
   end
   
   methods(TestMethodSetup)
      function setOpts(testCase)

         testCase.opts = icemodel.setopts( ...
            testCase.simmodel, ...
            testCase.sitename, ...
            testCase.simyears, ...
            testCase.forcings, ...
            testCase.userdata, ...
            testCase.uservars, ...
            testCase.savedata);
      end
   end

   methods(Test)
      function testSolver(testCase, solver)

         testCase.opts.solver = solver;
         
         % standard approach (for "slow" code)
         testCase.startMeasuring
         icemodel(testCase.opts);
         testCase.stopMeasuring
      end
   end
end