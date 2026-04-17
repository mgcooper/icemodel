classdef SpectralKernelPerfTest < matlab.perftest.TestCase
   %SPECTRALKERNELPERFTEST Benchmark the spectral-radiative kernel path.

   properties
      workspace
      state
      ro_sno
      ro_sno_spect
      k_bulk
      z_nodes
      z_nodes_spect
      z_edges_spect
      k_bulk_lookup
      k_bulk_lookup_empty
   end

   methods (TestClassSetup)
      function buildSpectralState(testCase)

         % Configure the workspace
         testCase.workspace = icemodel.test.fixtures.makeSyntheticWorkspace( ...
            2016, configure=true, nsteps=96, dt_seconds=900);

         % Build one synthetic spectral state and cache the exact and
         % approximate transforms used by the benchmarked paths.
         testCase.state = icemodel.test.fixtures.makeSyntheticColumnState( ...
            testCase.workspace, 'icemodel', solver=3, include_spectral=true, ...
            metstep=49, testname='spectral_perf_kernel');
         s = testCase.state;

         % Cache the thermal-grid density and the spectral-grid remap used by
         % every spectral source-term variant below.
         testCase.ro_sno = s.ro_ice * s.f_ice + s.ro_liq * s.f_liq + ...
            s.ro_air * (1 - s.f_ice - s.f_liq);
         testCase.z_nodes = s.z_nodes;
         testCase.z_nodes_spect = s.z_nodes_spect;
         testCase.z_edges_spect = s.z_edges_spect;
         testCase.ro_sno_spect = max(interp1(testCase.z_nodes, ...
            testCase.ro_sno, testCase.z_nodes_spect, ...
            'nearest', 'extrap'), 300);

         % Build the exact bulk extinction coefficients once for the narrow
         % bulk and two-stream benchmarks below.
         testCase.k_bulk = icemodel.radiation.bulk_extinction_coefficients( ...
            s.dz_spect, testCase.ro_sno_spect, s.tau_N, s.tau_S, s.solar_dwavel);

         % Precompute the lookup used by the fast-path exploratory benchmarks.
         testCase.k_bulk_lookup = icemodel.radiation.make_bulk_extinction_lookup( ...
            s.dz_spect, s.tau_N, s.tau_S, s.solar_dwavel);

         % Create an empty lookup table for the exact path.
         testCase.k_bulk_lookup_empty = struct([]);

         % Verify that the exact functions path reproduces the historical
         % inlined path before timing either branch.
         [Sc_inlined, chi_inlined] = SPECTRALSOURCETERM_INLINE(s.swd, ...
            s.albedo, s.I0, s.dz_spect, s.tau_N, s.tau_S, s.solar_dwavel, ...
            s.dz, testCase.ro_sno, testCase.z_nodes, testCase.z_nodes_spect);

         [Sc_functions, chi_functions] = icemodel.column.shortwave_source_term( ...
            s.swd, s.albedo, s.I0, s.dz_spect, s.tau_N, s.tau_S, s.solar_dwavel, ...
            s.dz, testCase.z_nodes, testCase.z_nodes_spect, testCase.z_edges_spect, ...
            testCase.ro_sno, testCase.k_bulk_lookup_empty);

         assert(max(abs(Sc_inlined - Sc_functions)) < 1e-10 * ...
            max(1, max(abs(Sc_inlined))));

         assert(abs(chi_inlined - chi_functions) < 1e-12);
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

         % Compute the bulk extinction coefficients and verify they're finite
         k_bulk_out = icemodel.radiation.bulk_extinction_coefficients( ...
            s.dz_spect, testCase.ro_sno_spect, s.tau_N, s.tau_S, s.solar_dwavel);
         testCase.assertTrue(all(isfinite(k_bulk_out)));

         while testCase.keepMeasuring
            for n = 1:batch_size
               k_bulk_out = icemodel.radiation.bulk_extinction_coefficients( ...
                  s.dz_spect, testCase.ro_sno_spect, s.tau_N, s.tau_S, ...
                  s.solar_dwavel);
            end
            if ~all(isfinite(k_bulk_out))
               error(['icemodel.radiation.bulk_extinction_coefficients ' ...
                  'benchmark produced a non-finite result'])
            end
         end
      end

      function testBulkExtCoefsLookup(testCase)
         % Benchmark the approximate lookup-table bulk-extinction path.
         s = testCase.state;
         batch_size = 1024;

         k_bulk_out = icemodel.radiation.bulk_extinction_coefficients( ...
            s.dz_spect, testCase.ro_sno_spect, s.tau_N, s.tau_S, ...
            s.solar_dwavel, testCase.k_bulk_lookup);

         testCase.assertTrue(all(isfinite(k_bulk_out)));

         while testCase.keepMeasuring
            for n = 1:batch_size
               k_bulk_out = icemodel.radiation.bulk_extinction_coefficients( ...
                  s.dz_spect, testCase.ro_sno_spect, s.tau_N, s.tau_S, ...
                  s.solar_dwavel, testCase.k_bulk_lookup);
            end
            if ~all(isfinite(k_bulk_out))
               error(['icemodel.radiation.bulk_extinction_coefficients ' ...
                  'lookup benchmark produced a non-finite result'])
            end
         end
      end

      function testTwoStreamSolve(testCase)
         % Benchmark just the two-stream linear solve on fixed coefficients.
         s = testCase.state;
         batch_size = 384;
         xynet = icemodel.radiation.solvetwostream(s.I0, s.albedo, ...
            testCase.k_bulk, testCase.z_edges_spect);
         testCase.assertTrue(all(isfinite(xynet)));

         while testCase.keepMeasuring
            for n = 1:batch_size
               xynet = icemodel.radiation.solvetwostream(s.I0, s.albedo, ...
                  testCase.k_bulk, testCase.z_edges_spect);
            end
            if ~all(isfinite(xynet))
               error(['icemodel.radiation.solvetwostream benchmark produced ' ...
                  'a non-finite result'])
            end
         end
      end

      function testSourceTermInlined(testCase)
         % Benchmark the historical inlined spectral source-term implementation.
         s = testCase.state;
         batch_size = 256;

         [Sc_new, chi] = SPECTRALSOURCETERM_INLINE(s.swd, s.albedo, s.I0, ...
            s.dz_spect, s.tau_N, s.tau_S, s.solar_dwavel, s.dz, ...
            testCase.ro_sno, testCase.z_nodes, testCase.z_nodes_spect);

         testCase.assertTrue(all(isfinite(Sc_new)));
         testCase.assertTrue(isfinite(chi));

         while testCase.keepMeasuring
            for n = 1:batch_size
               [Sc_new, chi] = SPECTRALSOURCETERM_INLINE(s.swd, s.albedo, ...
                  s.I0, s.dz_spect, s.tau_N, s.tau_S, s.solar_dwavel, s.dz, ...
                  testCase.ro_sno, testCase.z_nodes, testCase.z_nodes_spect);
            end
            if ~all(isfinite(Sc_new)) || ~isfinite(chi)
               error(['inlined spectral source-term benchmark produced a ', ...
                  'non-finite result'])
            end
         end
      end

      function testSourceTermFunctions(testCase)
         % Benchmark the organized exact spectral source-term path.
         s = testCase.state;
         batch_size = 256;

         [Sc_new, chi] = icemodel.column.shortwave_source_term(s.swd, s.albedo, ...
            s.I0, s.dz_spect, s.tau_N, s.tau_S, s.solar_dwavel, s.dz, ...
            testCase.z_nodes, testCase.z_nodes_spect, testCase.z_edges_spect, ...
            testCase.ro_sno, testCase.k_bulk_lookup_empty);

         testCase.assertTrue(all(isfinite(Sc_new)));
         testCase.assertTrue(isfinite(chi));

         while testCase.keepMeasuring
            for n = 1:batch_size
               [Sc_new, chi] = icemodel.column.shortwave_source_term(s.swd, ...
                  s.albedo, s.I0, s.dz_spect, s.tau_N, s.tau_S, s.solar_dwavel, ...
                  s.dz, testCase.z_nodes, testCase.z_nodes_spect, ...
                  testCase.z_edges_spect, testCase.ro_sno, ...
                  testCase.k_bulk_lookup_empty);
            end
            if ~all(isfinite(Sc_new)) || ~isfinite(chi)
               error(['functions spectral source-term benchmark produced a ', ...
                  'non-finite result'])
            end
         end
      end

      function testSourceTermLookup(testCase)
         % Benchmark the shared source-term path with lookup bulk coefficients.
         s = testCase.state;
         batch_size = 256;

         [Sc_new, chi] = icemodel.column.shortwave_source_term(s.swd, s.albedo, ...
            s.I0, s.dz_spect, s.tau_N, s.tau_S, s.solar_dwavel, s.dz, ...
            testCase.z_nodes, testCase.z_nodes_spect, testCase.z_edges_spect, ...
            testCase.ro_sno, testCase.k_bulk_lookup);

         testCase.assertTrue(all(isfinite(Sc_new)));
         testCase.assertTrue(isfinite(chi));

         while testCase.keepMeasuring
            for n = 1:batch_size
               [Sc_new, chi] = icemodel.column.shortwave_source_term(s.swd, ...
                  s.albedo, s.I0, s.dz_spect, s.tau_N, s.tau_S, s.solar_dwavel, ...
                  s.dz, testCase.z_nodes, testCase.z_nodes_spect, ...
                  testCase.z_edges_spect, testCase.ro_sno, testCase.k_bulk_lookup);
            end
            if ~all(isfinite(Sc_new)) || ~isfinite(chi)
               error(['lookup spectral source-term benchmark produced a ', ...
                  'non-finite result'])
            end
         end
      end
   end
end
