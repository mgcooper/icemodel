classdef ColumnKernelPerfTest < matlab.perftest.TestCase
   %COLUMNKERNELPERFTEST Benchmark representative column-solver kernels.

   properties
      workspace
      skin
      ice
   end

   methods (TestClassSetup)
      function buildSyntheticColumns(testCase)
         % Build one stable skin column and one stable ice column so each
         % benchmark reuses the same resolved kernel inputs.
         testCase.workspace = icemodel.test.fixtures.makeSyntheticWorkspace( ...
            2016, configure=true, nsteps=24, dt_seconds=900);
         testCase.skin = icemodel.test.fixtures.makeSyntheticColumnState( ...
            testCase.workspace, 'skinmodel', solver=1, ...
            testname='skin_perf_kernel');
         testCase.ice = icemodel.test.fixtures.makeSyntheticColumnState( ...
            testCase.workspace, 'icemodel', solver=3, ...
            testname='ice_perf_kernel');
      end
   end

   methods (TestClassTeardown)
      function cleanupSyntheticColumns(testCase)
         % Tear down the synthetic workspace after the class finishes.
         icemodel.test.fixtures.cleanupSyntheticWorkspace(testCase.workspace);
      end
   end

   methods (Test)
      function testUpdateState(testCase)
         % Benchmark the state update used by the enthalpy solve.
         s = testCase.ice;

         % Batch this fast calculation to reduce measurement noise.
         batch_size = 4096;

         [H, k_eff, dHdT, dFdT, drovdT, ro_vap] = ...
            icemodel.column.updatestate(s.T_ice, s.f_ice, s.f_liq, s.f_wat);
         testCase.assertTrue(all(isfinite([H; k_eff; dHdT; dFdT; ...
            drovdT; ro_vap])));

         while testCase.keepMeasuring
            for n = 1:batch_size
               [H, k_eff, dHdT, dFdT, drovdT, ro_vap] = ...
                  icemodel.column.updatestate( ...
                  s.T_ice, s.f_ice, s.f_liq, s.f_wat);
            end
            if ~all(isfinite([H; k_eff; dHdT; dFdT; drovdT; ro_vap]))
               error('updatestate benchmark produced a non-finite result')
            end
         end
      end

      function testSkinSolve(testCase)
         % Benchmark the standalone skin-column solve on a fixed state.
         s = testCase.skin;

         % icemodel.column.solve_column_temperature is already coarse enough
         % that a modest batch gives a stable sample without stretching the
         % suite runtime.
         batch_size = 64;

         % Keep failure dumps off so the timing measures only the solve.
         settings = s.settings;
         settings.debug = false;

         [T_ice, f_ice, f_liq, k_eff, ok] = ...
            icemodel.column.solve_column_temperature(s.T_sfc, s.T_ice, ...
            s.f_ice, s.f_liq, s.dz, s.delz, s.fn, s.opts.dt, settings);

         testCase.assertTrue(ok);
         testCase.assertTrue(all(isfinite([T_ice; f_ice; f_liq; k_eff])));

         while testCase.keepMeasuring
            for n = 1:batch_size
               [T_ice, f_ice, f_liq, k_eff, ok] = ...
                  icemodel.column.solve_column_temperature(s.T_sfc, s.T_ice, ...
                  s.f_ice, s.f_liq, s.dz, s.delz, s.fn, s.opts.dt, ...
                  settings);
            end
            if ~ok || ~all(isfinite([T_ice; f_ice; f_liq; k_eff]))
               error('solve_column_temperature benchmark failed to converge')
            end
         end
      end

      function testSkinEbSolve(testCase)
         % Benchmark the coupled skin energy-balance solve.
         s = testCase.skin;

         % The coupled skin solve has moderate jitter, so use a wider batch
         % than the standalone solve to keep strict-profile runs stable.
         batch_size = 64;

         [T_sfc, T_ice, f_ice, f_liq, k_eff, diag] = ...
            icemodel.couplers.solve_skin_surface_column( ...
            s.T_sfc, s.T_ice, s.f_ice, s.f_liq, s.dz, s.delz, s.fn, ...
            s.opts.dt, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, ...
            s.tppt, s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, ...
            s.H_e, s.hv_atm, s.br_coefs, s.liqflag, s.chi, s.ro_sfc, ...
            s.snow_depth, s.settings, s.opts);

         testCase.assertTrue(diag.ok_seb && diag.ok_ieb && diag.ok_cpl);
         testCase.assertTrue( ...
            all(isfinite([T_sfc; T_ice; f_ice; f_liq; k_eff])));

         while testCase.keepMeasuring
            for n = 1:batch_size
               [T_sfc, T_ice, f_ice, f_liq, k_eff, diag] = ...
                  icemodel.couplers.solve_skin_surface_column( ...
                  s.T_sfc, s.T_ice, s.f_ice, s.f_liq, s.dz, s.delz, s.fn, ...
                  s.opts.dt, s.tair, s.swd, s.lwd, s.albedo, s.wspd, ...
                  s.ppt, s.tppt, s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, ...
                  s.nu_air, s.H_h, s.H_e, s.hv_atm, s.br_coefs, s.liqflag, ...
                  s.chi, s.ro_sfc, s.snow_depth, s.settings, s.opts);
            end
            if ~(diag.ok_seb && diag.ok_ieb && diag.ok_cpl) || ...
                  ~all(isfinite([T_sfc; T_ice; f_ice; f_liq; k_eff]))
               error('solve_skin_surface_column benchmark failed to converge')
            end
         end
      end

      function testIceEnbal(testCase)
         % Benchmark the interior enthalpy/phase solve for one full column.
         s = testCase.ice;

         % icemodel.column.solve_column_enthalpy is one of the heavier
         % column kernels but still benefits from a wider batch to hit the
         % runner's error target reliably.
         batch_size = 256;

         % Solve with a Dirichlet upper boundary. Keep failure dumps off so
         % the timing measures only the solve.
         settings = s.settings;
         settings.solver = 1;
         settings.debug = false;

         [T_ice, f_ice, f_liq, k_eff, ~, ~, ok] = ...
            icemodel.column.solve_column_enthalpy( ...
            s.T_sfc, s.T_ice, s.f_ice, s.f_liq, s.Fc, s.Fp, s.Sc, s.Sp, ...
            s.dz, s.delz, s.fn, s.opts.dt, s.opts.f_res_pore_ice, settings);

         testCase.assertTrue(ok);
         testCase.assertTrue(all(isfinite([T_ice; f_ice; f_liq; k_eff])));

         while testCase.keepMeasuring
            for n = 1:batch_size
               [T_ice, f_ice, f_liq, k_eff, ~, ~, ok] = ...
                  icemodel.column.solve_column_enthalpy(s.T_sfc, s.T_ice, ...
                  s.f_ice, s.f_liq, s.Fc, s.Fp, s.Sc, s.Sp, s.dz, ...
                  s.delz, s.fn, s.opts.dt, s.opts.f_res_pore_ice, settings);
            end
            if ~ok || ~all(isfinite([T_ice; f_ice; f_liq; k_eff]))
               error('solve_column_enthalpy benchmark failed to converge')
            end
         end
      end

      function testIceEbSolve(testCase)
         % Benchmark the coupled ice-column + SEB solve.
         s = testCase.ice;

         % The coupled Robin icemodel solve carries both column and SEB
         % work, so scale the batch enough to keep variance low without
         % making the suite too slow.
         batch_size = 128;

         [T_sfc, T_ice, f_ice, f_liq, k_eff, ~, ~, diag] = ...
            icemodel.couplers.solve_surface_column_robin(s.T_sfc, s.T_ice, ...
            s.f_ice, s.f_liq, s.Sc, s.Sp, s.dz, s.delz, s.fn, s.opts.dt, ...
            s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, ...
            s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, s.hv_atm, ...
            s.br_coefs, s.liqflag, s.chi, s.ro_sfc, s.snow_depth, ...
            s.opts.f_res_pore_ice, s.settings, s.opts);

         testCase.assertTrue(diag.ok_ieb);
         testCase.assertTrue( ...
            all(isfinite([T_sfc; T_ice; f_ice; f_liq; k_eff])));

         while testCase.keepMeasuring
            for n = 1:batch_size
               [T_sfc, T_ice, f_ice, f_liq, k_eff, ~, ~, diag] = ...
                  icemodel.couplers.solve_surface_column_robin(s.T_sfc, ...
                  s.T_ice, s.f_ice, s.f_liq, s.Sc, s.Sp, s.dz, s.delz, s.fn, ...
                  s.opts.dt, s.tair, s.swd, s.lwd, s.albedo, s.wspd, ...
                  s.ppt, s.tppt, s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, ...
                  s.nu_air, s.H_h, s.H_e, s.hv_atm, s.br_coefs, ...
                  s.liqflag, s.chi, s.ro_sfc, s.snow_depth, ...
                  s.opts.f_res_pore_ice, s.settings, s.opts);
            end
            if ~diag.ok_ieb || ...
                  ~all(isfinite([T_sfc; T_ice; f_ice; f_liq; k_eff]))
               error('surface_column_robin benchmark failed to converge')
            end
         end
      end
   end
end
