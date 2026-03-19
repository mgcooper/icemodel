classdef SpectralKernelPerfTest < matlab.perftest.TestCase
   %SPECTRALKERNELPERFTEST Benchmark the spectral-radiative kernel path.

   properties
      workspace
      state
      bulkcoefs
      a
      r
      z_walls
      ro_sno
   end

   methods (TestClassSetup)
      function buildSpectralState(testCase)
         % Build one synthetic spectral state and cache the derived density
         % transform used by the two benchmarked spectral paths.
         testCase.workspace = icemodel.test.fixtures.makeSyntheticWorkspace( ...
            2016, configure=true, nsteps=24, dt_seconds=900);
         testCase.state = icemodel.test.fixtures.makeSyntheticColumnState( ...
            testCase.workspace, 'icemodel', solver=3, ...
            include_spectral=true, testname='spectral_perf_kernel');

         s = testCase.state;

         % Recreate the current UPDATEEXTCOEFS density-to-bulkcoefs path so
         % the microbenchmarks exercise the same spectral bottleneck.
         testCase.ro_sno = s.ro_ice * s.f_ice + s.ro_liq * s.f_liq + ...
            s.ro_air * (1 - s.f_ice - s.f_liq);
         ro_sno_spect = max(interp1(cumsum(s.dz) - s.dz / 2, ...
            testCase.ro_sno, cumsum(s.dz_spect) - s.dz_spect / 2, ...
            'nearest', 'extrap'), 300);
         testCase.bulkcoefs = -log((sum(s.solardwavl .* ...
            exp(s.spect_S .* ro_sno_spect), 2)) ./ ...
            (sum(s.solardwavl .* exp(s.spect_N .* ro_sno_spect), 2))) ...
            ./ s.dz_spect(1);
         testCase.bulkcoefs = vertcat(testCase.bulkcoefs, ...
            testCase.bulkcoefs(end), testCase.bulkcoefs(end));
         [testCase.a, testCase.r] = GETAANDR(testCase.bulkcoefs, s.albedo);
         [~, ~, ~, testCase.z_walls] = CVMESH(s.opts.z0_spectral, ...
            s.opts.dz_spectral);
      end
   end

   methods (TestClassTeardown)
      function cleanupSpectralState(testCase)
         % Tear down the synthetic workspace after the class finishes.
         icemodel.test.fixtures.cleanupSyntheticWorkspace(testCase.workspace);
      end
   end

   methods (Test)
      function testTwoStreamSolve(testCase)
         % Benchmark just the two-stream linear solve on fixed coefficients.
         s = testCase.state;
         % The bare two-stream solve is fast enough that a wide batch is
         % needed for stable error-controlled timing.
         batch_size = 384;
         xynet = SOLVETWOSTREAM(testCase.a, testCase.r, ...
            testCase.bulkcoefs, s.Q0, s.albedo, testCase.z_walls);
         testCase.assertTrue(all(isfinite(xynet)));

         while testCase.keepMeasuring
            for n = 1:batch_size
               xynet = SOLVETWOSTREAM(testCase.a, testCase.r, ...
                  testCase.bulkcoefs, s.Q0, s.albedo, testCase.z_walls);
            end
            if ~all(isfinite(xynet))
               error('SOLVETWOSTREAM benchmark produced a non-finite result')
            end
         end
      end

      function testUpdateExtCoefs(testCase)
         % Benchmark the inline spectral coefficient update path used during
         % production runs.
         s = testCase.state;
         % UPDATEEXTCOEFS is already dominated by the exponential spectral
         % transform, so keep the batch large but bounded to preserve a
         % short overall suite runtime.
         batch_size = 8192;
         [Sc_new, chi] = UPDATEEXTCOEFS(s.swd, s.albedo, s.Q0, ...
            s.dz_spect, s.spect_N, s.spect_S, s.solardwavl, s.Sc, s.dz, ...
            testCase.ro_sno);
         testCase.assertTrue(all(isfinite(Sc_new)));
         testCase.assertTrue(isfinite(chi));

         while testCase.keepMeasuring
            for n = 1:batch_size
               [Sc_new, chi] = UPDATEEXTCOEFS(s.swd, s.albedo, s.Q0, ...
                  s.dz_spect, s.spect_N, s.spect_S, s.solardwavl, s.Sc, ...
                  s.dz, testCase.ro_sno);
            end
            if ~all(isfinite(Sc_new)) || ~isfinite(chi)
               error('UPDATEEXTCOEFS benchmark produced a non-finite result')
            end
         end
      end
   end
end
