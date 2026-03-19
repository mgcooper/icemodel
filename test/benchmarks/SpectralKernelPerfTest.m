classdef SpectralKernelPerfTest < matlab.perftest.TestCase
   %SPECTRALKERNELPERFTEST Benchmark the spectral-radiative kernel path.

   properties
      workspace
      state
      ro_sno
      ro_sno_spect
      bulkcoefs
      a
      r
      z_walls
      grid_thermal
      grid_spectral
      bulk_lookup
   end

   methods (TestClassSetup)
      function buildSpectralState(testCase)
         % Build one synthetic spectral state and cache the exact and
         % approximate transforms used by the benchmarked paths.
         testCase.workspace = icemodel.test.fixtures.makeSyntheticWorkspace( ...
            2016, configure=true, nsteps=96, dt_seconds=900);
         testCase.state = icemodel.test.fixtures.makeSyntheticColumnState( ...
            testCase.workspace, 'icemodel', solver=3, ...
            include_spectral=true, metstep=49, ...
            testname='spectral_perf_kernel');

         s = testCase.state;

         % Cache the thermal-grid density and the spectral-grid remap used
         % by every UPDATEEXTCOEFS variant below.
         testCase.ro_sno = s.ro_ice * s.f_ice + s.ro_liq * s.f_liq + ...
            s.ro_air * (1 - s.f_ice - s.f_liq);
         testCase.grid_thermal = cumsum(s.dz) - s.dz / 2;
         testCase.grid_spectral = cumsum(s.dz_spect) - s.dz_spect / 2;
         testCase.ro_sno_spect = max(GRIDFORWARD(testCase.ro_sno, ...
            testCase.grid_thermal, testCase.grid_spectral), 300);

         % Build the exact bulk coefficients and two-stream coefficients once
         % for the narrower solve-only benchmark below.
         testCase.bulkcoefs = BULKEXTCOEFS(s.dz_spect, ...
            testCase.ro_sno_spect, s.spect_N, s.spect_S, s.solardwavl);
         [testCase.a, testCase.r] = GETAANDR(testCase.bulkcoefs, s.albedo);
         [~, ~, ~, testCase.z_walls] = CVMESH(s.opts.z0_spectral, ...
            s.opts.dz_spectral);

         % Precompute the approximate density lookup used by the fast-path
         % exploratory benchmarks.
         testCase.bulk_lookup = MAKEBULKEXTCOEFSLOOKUP(s.dz_spect, ...
            s.spect_N, s.spect_S, s.solardwavl);

         % Verify that the exact function-call decomposition reproduces the
         % current inlined UPDATEEXTCOEFS path before timing either branch.
         [Sc_inline, chi_inline] = UPDATEEXTCOEFSINLINELEGACY(s.swd, ...
            s.albedo, s.Q0, s.dz_spect, s.spect_N, s.spect_S, ...
            s.solardwavl, s.Sc, s.dz, testCase.ro_sno);
         [Sc_function, chi_function] = UPDATEEXTCOEFSDECOMPOSED(s.swd, ...
            s.albedo, s.Q0, s.dz_spect, s.spect_N, s.spect_S, ...
            s.solardwavl, s.Sc, s.dz, testCase.ro_sno);
         assert(max(abs(Sc_inline - Sc_function)) < 1e-10 * ...
            max(1, max(abs(Sc_inline))));
         assert(abs(chi_inline - chi_function) < 1e-12);
      end
   end

   methods (TestClassTeardown)
      function cleanupSpectralState(testCase)
         % Tear down the synthetic workspace after the class finishes.
         icemodel.test.fixtures.cleanupSyntheticWorkspace(testCase.workspace);
      end
   end

   methods (Test)
      function testBulkExtCoefsExact(testCase)
         % Benchmark the exact exponential/log bulk-extinction transform.
         s = testCase.state;
         batch_size = 32;
         bulkcoefs_out = BULKEXTCOEFS(s.dz_spect, testCase.ro_sno_spect, ...
            s.spect_N, s.spect_S, s.solardwavl);
         testCase.assertTrue(all(isfinite(bulkcoefs_out)));

         while testCase.keepMeasuring
            for n = 1:batch_size
               bulkcoefs_out = BULKEXTCOEFS(s.dz_spect, testCase.ro_sno_spect, ...
                  s.spect_N, s.spect_S, s.solardwavl);
            end
            if ~all(isfinite(bulkcoefs_out))
               error('BULKEXTCOEFS benchmark produced a non-finite result')
            end
         end
      end

      function testBulkExtCoefsLookup(testCase)
         % Benchmark the approximate lookup-table bulk-extinction path.
         batch_size = 1024;
         bulkcoefs_out = BULKEXTCOEFSLOOKUP(testCase.ro_sno_spect, ...
            testCase.bulk_lookup);
         testCase.assertTrue(all(isfinite(bulkcoefs_out)));

         while testCase.keepMeasuring
            for n = 1:batch_size
               bulkcoefs_out = BULKEXTCOEFSLOOKUP(testCase.ro_sno_spect, ...
                  testCase.bulk_lookup);
            end
            if ~all(isfinite(bulkcoefs_out))
               error(['BULKEXTCOEFSLOOKUP benchmark produced a non-finite ', ...
                  'result'])
            end
         end
      end

      function testTwoStreamSolve(testCase)
         % Benchmark just the two-stream linear solve on fixed coefficients.
         s = testCase.state;
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

      function testUpdateExtCoefsInlineLegacy(testCase)
         % Benchmark the historical inlined UPDATEEXTCOEFS implementation.
         s = testCase.state;
         batch_size = 256;
         [Sc_new, chi] = UPDATEEXTCOEFSINLINELEGACY(s.swd, s.albedo, ...
            s.Q0, s.dz_spect, s.spect_N, s.spect_S, s.solardwavl, s.Sc, ...
            s.dz, testCase.ro_sno);
         testCase.assertTrue(all(isfinite(Sc_new)));
         testCase.assertTrue(isfinite(chi));

         while testCase.keepMeasuring
            for n = 1:batch_size
               [Sc_new, chi] = UPDATEEXTCOEFSINLINELEGACY(s.swd, s.albedo, ...
                  s.Q0, s.dz_spect, s.spect_N, s.spect_S, s.solardwavl, ...
                  s.Sc, s.dz, testCase.ro_sno);
            end
            if ~all(isfinite(Sc_new)) || ~isfinite(chi)
               error(['inline UPDATEEXTCOEFS benchmark produced a ', ...
                  'non-finite result'])
            end
         end
      end

      function testUpdateExtCoefsFunctionExact(testCase)
         % Benchmark the exact drop-in helper-call replacement for UPDATEEXTCOEFS.
         s = testCase.state;
         batch_size = 256;
         [Sc_new, chi] = UPDATEEXTCOEFSDECOMPOSED(s.swd, s.albedo, ...
            s.Q0, s.dz_spect, s.spect_N, s.spect_S, s.solardwavl, ...
            s.Sc, s.dz, testCase.ro_sno);
         testCase.assertTrue(all(isfinite(Sc_new)));
         testCase.assertTrue(isfinite(chi));

         while testCase.keepMeasuring
            for n = 1:batch_size
               [Sc_new, chi] = UPDATEEXTCOEFSDECOMPOSED(s.swd, s.albedo, ...
                  s.Q0, s.dz_spect, s.spect_N, s.spect_S, ...
                  s.solardwavl, s.Sc, s.dz, testCase.ro_sno);
            end
            if ~all(isfinite(Sc_new)) || ~isfinite(chi)
               error(['function UPDATEEXTCOEFS benchmark produced a ', ...
                  'non-finite result'])
            end
         end
      end

      function testUpdateExtCoefsFunctionExactCached(testCase)
         % Benchmark the exact helper-call path with cached spectral geometry.
         s = testCase.state;
         batch_size = 256;
         [Sc_new, chi] = UPDATEEXTCOEFSDECOMPOSEDCACHED(s.swd, s.albedo, ...
            s.Q0, s.dz_spect, s.spect_N, s.spect_S, s.solardwavl, s.Sc, ...
            s.dz, testCase.ro_sno, testCase.grid_thermal, ...
            testCase.grid_spectral, testCase.z_walls);
         testCase.assertTrue(all(isfinite(Sc_new)));
         testCase.assertTrue(isfinite(chi));

         while testCase.keepMeasuring
            for n = 1:batch_size
               [Sc_new, chi] = UPDATEEXTCOEFSDECOMPOSEDCACHED(s.swd, ...
                  s.albedo, s.Q0, s.dz_spect, s.spect_N, s.spect_S, ...
                  s.solardwavl, s.Sc, s.dz, testCase.ro_sno, ...
                  testCase.grid_thermal, testCase.grid_spectral, ...
                  testCase.z_walls);
            end
            if ~all(isfinite(Sc_new)) || ~isfinite(chi)
               error(['cached function UPDATEEXTCOEFS benchmark produced a ', ...
                  'non-finite result'])
            end
         end
      end

      function testUpdateExtCoefsFunctionLookup(testCase)
         % Benchmark the exploratory function-call path with lookup bulkcoefs.
         s = testCase.state;
         batch_size = 256;
         [Sc_new, chi] = UPDATEEXTCOEFSLOOKUP(s.swd, s.albedo, ...
            s.Q0, s.dz_spect, s.Sc, s.dz, testCase.ro_sno, ...
            testCase.grid_thermal, testCase.grid_spectral, ...
            testCase.z_walls, testCase.bulk_lookup);
         testCase.assertTrue(all(isfinite(Sc_new)));
         testCase.assertTrue(isfinite(chi));

         while testCase.keepMeasuring
            for n = 1:batch_size
               [Sc_new, chi] = UPDATEEXTCOEFSLOOKUP(s.swd, s.albedo, ...
                  s.Q0, s.dz_spect, s.Sc, s.dz, testCase.ro_sno, ...
                  testCase.grid_thermal, testCase.grid_spectral, ...
                  testCase.z_walls, testCase.bulk_lookup);
            end
            if ~all(isfinite(Sc_new)) || ~isfinite(chi)
               error(['lookup UPDATEEXTCOEFS benchmark produced a ', ...
                  'non-finite result'])
            end
         end
      end
   end
end
