classdef SebSolverTest < matlab.perftest.TestCase

   properties
      T
      tair
      swd
      lwd
      albedo
      psfc
      wspd
      rh
      ea
      k_eff
      liqflag
      De
      wcoef
      scoef
      dz
      Qc
      chi
   end

   properties (TestParameter)
      solver = {1, 2, 3}
   end

   methods(TestMethodSetup)
      function generateTestData(testCase)

         testCase.liqflag = false;
         testCase.tair = 274;
         testCase.swd = 500;
         testCase.lwd = 300;
         testCase.albedo = 0.5;
         testCase.psfc = 8600;
         testCase.wspd = 5;
         testCase.rh = 80;
         testCase.ea = VAPPRESS(testCase.tair, 273.16, testCase.liqflag) ...
            .* testCase.rh/100;
         testCase.k_eff = 2;

         z_0 = 1e-3;
         z_obs = 3;
         [De, wcoef, scoef] = WINDCOEF(testCase.wspd, z_0, z_obs);

         testCase.De = De;
         testCase.wcoef = wcoef;
         testCase.scoef = scoef;

         testCase.dz = 0.04 * ones(500, 1);

         testCase.Qc = 0;
         testCase.chi = 0.5;

         Z = cumsum(testCase.dz);
         testCase.T = testCase.tair * ones(500, 1);
         testCase.T = testCase.T .* exp(-0.1 * Z./Z(end));
      end
   end

   methods(Test)
      function testSolver(testCase, solver)

         [cv_air, emiss, SB, Tf, roL] = icemodel.physicalConstant(...
            "cv_air", "emiss", "SB", "Tf", "roLv");

         % standard approach (for "slow" code)
         while(testCase.keepMeasuring)
            Tsfc = SEBSOLVE(testCase.tair, testCase.swd, testCase.lwd, ...
               testCase.albedo, testCase.De, testCase.psfc, testCase.wspd, ...
               testCase.ea, cv_air, emiss, SB, Tf, testCase.chi, roL, ...
               testCase.scoef, testCase.Qc, testCase.liqflag, testCase.tair, ...
               testCase.T, testCase.k_eff, testCase.dz, solver);
         end
      end
   end
end