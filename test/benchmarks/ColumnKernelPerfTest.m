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
         % Benchmark the state-update algebra used by the enthalpy solve.
         s = testCase.ice;

         % icemodel.column.updatestate is extremely fast, so each timed sample
         % needs a large batch to avoid framework noise at the default and
         % strict sampling targets.
         batch_size = 4096;

         [H, k_eff, dHdT, dLdT, drovdT, ro_vap] = ...
            icemodel.column.updatestate(s.T, s.f_ice, s.f_liq, s.f_wat);

         testCase.assertTrue(all(isfinite([H; k_eff; dHdT; dLdT; drovdT; ...
            ro_vap])));

         while testCase.keepMeasuring
            for n = 1:batch_size
               [H, k_eff, dHdT, dLdT, drovdT, ro_vap] = ...
                  icemodel.column.updatestate(s.T, s.f_ice, s.f_liq, ...
                  s.f_wat);
            end
            if ~all(isfinite([H; k_eff; dHdT; dLdT; drovdT; ro_vap]))
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

         [T, f_ice, f_liq, k_eff, ok] = ...
            icemodel.column.solve_column_temperature(s.T, s.f_ice, ...
            s.f_liq, s.dz, s.delz, s.fn, s.opts.dt, s.Ts, s.tol, ...
            s.maxiter, s.alpha, false);

         testCase.assertTrue(ok);
         testCase.assertTrue(all(isfinite([T; f_ice; f_liq; k_eff])));

         while testCase.keepMeasuring
            for n = 1:batch_size
               [T, f_ice, f_liq, k_eff, ok] = ...
                  icemodel.column.solve_column_temperature(s.T, s.f_ice, ...
                  s.f_liq, s.dz, s.delz, s.fn, s.opts.dt, s.Ts, s.tol, ...
                  s.maxiter, s.alpha, false);
            end
            if ~ok || ~all(isfinite([T; f_ice; f_liq; k_eff]))
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

         [Ts, T, f_ice, f_liq, k_eff, ok_seb, ok_ieb, ok] = SKINEBSOLVE( ...
            s.Ts, s.T, s.f_ice, s.f_liq, s.dz, s.delz, s.fn, s.opts.dt, ...
            s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, ...
            s.psfc, s.De, s.ea_atm, s.chi, s.ro_air_Lv, s.br_coefs, s.liqflag, ...
            s.seb_solver, s.tol, s.maxiter, s.alpha, s.cpl_maxiter, ...
            s.cpl_Ts_tol, s.cpl_seb_tol, s.cpl_alpha, s.cpl_aitken, ...
            s.cpl_jumpmax, false, s.ro_sfc, s.snow_depth, s.opts);

         testCase.assertTrue(ok_seb && ok_ieb && ok);
         testCase.assertTrue(all(isfinite([Ts; T; f_ice; f_liq; k_eff])));

         while testCase.keepMeasuring
            for n = 1:batch_size
               [Ts, T, f_ice, f_liq, k_eff, ok_seb, ok_ieb, ok] = ...
                  SKINEBSOLVE(s.Ts, s.T, s.f_ice, s.f_liq, s.dz, ...
                  s.delz, s.fn, s.opts.dt, s.tair, s.swd, ...
                  s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, s.De, ...
                  s.ea_atm, s.chi, s.ro_air_Lv, s.br_coefs, s.liqflag, ...
                  s.seb_solver, s.tol, s.maxiter, s.alpha, ...
                  s.cpl_maxiter, s.cpl_Ts_tol, s.cpl_seb_tol, ...
                  s.cpl_alpha, s.cpl_aitken, s.cpl_jumpmax, false, ...
                  s.ro_sfc, s.snow_depth, s.opts);
            end
            if ~(ok_seb && ok_ieb && ok) || ...
                  ~all(isfinite([Ts; T; f_ice; f_liq; k_eff]))
               error('SKINEBSOLVE benchmark failed to converge')
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

         [T, f_ice, f_liq, k_eff, ok] = icemodel.column.solve_column_enthalpy( ...
            s.Ts, s.T, s.f_ice, s.f_liq, s.Fc, s.Fp, s.Sc, s.Sp, s.dz, ...
            s.delz, s.fn, s.opts.dt, 1, s.tol, s.maxiter, s.alpha, ...
            s.use_aitken, s.jumpmax, false);

         testCase.assertTrue(ok);
         testCase.assertTrue(all(isfinite([T; f_ice; f_liq; k_eff])));

         while testCase.keepMeasuring
            for n = 1:batch_size
               [T, f_ice, f_liq, k_eff, ok] = ...
                  icemodel.column.solve_column_enthalpy(s.Ts, s.T, ...
                  s.f_ice, s.f_liq, s.Fc, s.Fp, s.Sc, s.Sp, s.dz, ...
                  s.delz, s.fn, s.opts.dt, 1, s.tol, s.maxiter, s.alpha, ...
                  s.use_aitken, s.jumpmax, false);
            end
            if ~ok || ~all(isfinite([T; f_ice; f_liq; k_eff]))
               error('solve_column_enthalpy benchmark failed to converge')
            end
         end
      end

      function testIceEbSolve(testCase)
         % Benchmark the coupled ice-column + SEB solve.
         s = testCase.ice;

         % The coupled Robin icemodel solve carries both column and SEB
         % work, so scale the batch enough to keep variance low without
         % making the suite sluggish.
         batch_size = 128;

         [Ts, T, f_ice, f_liq, k_eff, ok_ieb] = ...
            icemodel.couplers.solve_surface_column_robin(s.Ts, s.T, s.f_ice, ...
            s.f_liq, s.dz, s.delz, s.fn, s.Sc, s.Sp, s.opts.dt, s.tair, ...
            s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, s.De, ...
            s.ea_atm, s.br_coefs, s.ro_air_Lv, s.liqflag, s.chi, 3, s.tol, ...
            s.maxiter, s.alpha, s.use_aitken, s.jumpmax, s.cpl_Ts_tol, ...
            s.cpl_seb_tol, s.cpl_maxiter, s.cpl_alpha, s.cpl_aitken, ...
            s.cpl_jumpmax, s.ro_sfc, s.snow_depth, s.opts);

         testCase.assertTrue(ok_ieb);
         testCase.assertTrue(all(isfinite([Ts; T; f_ice; f_liq; k_eff])));

         while testCase.keepMeasuring
            for n = 1:batch_size
               [Ts, T, f_ice, f_liq, k_eff, ok_ieb] = ...
                  icemodel.couplers.solve_surface_column_robin(s.Ts, s.T, ...
                  s.f_ice, s.f_liq, s.dz, s.delz, s.fn, s.Sc, s.Sp, s.opts.dt, ...
                  s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, ...
                  s.psfc, s.De, s.ea_atm, s.br_coefs, s.ro_air_Lv, s.liqflag, ...
                  s.chi, 3, s.tol, s.maxiter, s.alpha, s.use_aitken, ...
                  s.jumpmax, s.cpl_Ts_tol, s.cpl_seb_tol, s.cpl_maxiter, ...
                  s.cpl_alpha, s.cpl_aitken, s.cpl_jumpmax, ...
                  s.ro_sfc, s.snow_depth, s.opts);
            end
            if ~ok_ieb || ~all(isfinite([Ts; T; f_ice; f_liq; k_eff]))
               error('surface_column_robin benchmark failed to converge')
            end
         end
      end
   end
end
