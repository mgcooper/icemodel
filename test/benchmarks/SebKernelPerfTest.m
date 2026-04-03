classdef SebKernelPerfTest < matlab.perftest.TestCase
   %SEBKERNELPERFTEST Benchmark the standalone SEB surface-solve kernel.

   properties
      T
      Ta
      Qsi
      Qli
      albedo
      Pa
      wspd
      ppt
      tppt
      rh
      ea_atm
      k_eff
      liqflag
      De
      br_coefs
      dz
      chi
      Ts0
      cv_air
      cv_liq
      emiss
      SB
      Tf
      roL
   end

   properties (TestParameter)
      seb_solver = {0, 1, 2}
   end

   methods (TestMethodSetup)
      function generateTestData(testCase)
         % Build one representative SEB state and reuse it across samples.
         testCase.liqflag = false;
         testCase.Ta = 274;
         testCase.Qsi = 500;
         testCase.Qli = 300;
         testCase.albedo = 0.5;
         testCase.Pa = 86000;
         testCase.wspd = 5;
         testCase.ppt = 0;
         testCase.tppt = 274;
         testCase.rh = 80;
         testCase.k_eff = 2;
         testCase.chi = 0.5;
         testCase.Ts0 = testCase.Ta;
         testCase.ea_atm = VAPPRESS(testCase.Ta, testCase.liqflag) ...
            .* testCase.rh / 100;

         [testCase.cv_air, testCase.cv_liq, ...
            testCase.SB, testCase.Tf, testCase.roL] = ...
            icemodel.physicalConstant( ...
            "cv_air", "cv_liq", "SB", "Tf", "roLv");
         testCase.emiss = icemodel.parameterLookup('emiss');

         z0_bulk = 1e-3;
         z_tair = 3;
         [testCase.De, testCase.br_coefs] = WINDCOEF(testCase.wspd, z0_bulk, z_tair);

         testCase.dz = 0.04 * ones(500, 1);
         z = cumsum(testCase.dz);
         testCase.T = testCase.Ta * ones(500, 1);
         testCase.T = testCase.T .* exp(-0.1 * z ./ z(end));
      end
   end

   methods (Test)
      function testSolver(testCase, seb_solver)
         % Benchmark one SEBSOLVE root-finder variant on a fixed surface
         % state so only solver-path overhead changes between cases.
         batch_size = sebBenchmarkBatchSize(seb_solver);

         % Benchmark only the SEBSOLVE call for each available root finder.
         [Ts, ok] = SEBSOLVE(testCase.Ts0, testCase.Ta, testCase.Qsi, ...
            testCase.Qli, testCase.albedo, testCase.wspd, testCase.ppt, ...
            testCase.tppt, testCase.Pa, testCase.De, testCase.ea_atm, ...
            testCase.chi, testCase.roL, testCase.br_coefs, ...
            testCase.liqflag, testCase.T, testCase.k_eff, testCase.dz, ...
            [], [], seb_solver, false, []);
         testCase.assertTrue(ok, 'SEBSOLVE setup did not converge before timing')
         testCase.assertThat(Ts, matlab.unittest.constraints.IsFinite)

         while testCase.keepMeasuring
            for n = 1:batch_size
               [Ts, ok] = SEBSOLVE(testCase.Ts0, testCase.Ta, ...
                  testCase.Qsi, testCase.Qli, testCase.albedo, ...
                  testCase.wspd, testCase.ppt, testCase.tppt, testCase.Pa, ...
                  testCase.De, testCase.ea_atm, testCase.chi, testCase.roL, ...
                  testCase.br_coefs, testCase.liqflag, testCase.T, ...
                  testCase.k_eff, testCase.dz, [], [], seb_solver, false, ...
                  []);
            end

            if ~ok || ~isfinite(Ts)
               error('SEBSOLVE benchmark failed for seb_solver=%d', seb_solver)
            end
         end
      end
   end
end

function batch_size = sebBenchmarkBatchSize(seb_solver)
   %SEBBENCHMARKBATCHSIZE Return calibrated batch sizes for SEB perf cases.

   switch seb_solver
      case 1
         % The secant/newton variant stays quick enough that it still needs
         % a large batch to stay above framework precision.
         batch_size = 2048;
      case 2
         % The constrained solve is similarly quick once the setup state is
         % fixed, so batch it aggressively as well.
         batch_size = 2048;
      otherwise
         % The fallback solver is also fast relative to the framework
         % clock, so use the same scale as the other SEB variants.
         batch_size = 2048;
   end
end
