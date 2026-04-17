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
      br_coefs
      dz
      chi
      Ts0
      ro_atm
      cv_atm
      nu_air
      H_h
      H_e
      hv_atm
      ro_sfc
      snow_depth
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
         testCase.ea_atm = icemodel.vapor.saturation_vapor_pressure( ...
            testCase.Ta, testCase.liqflag) .* testCase.rh / 100;

         [Rd, cp_air, Lv, epsilon] = icemodel.physicalConstant( ...
            'Rd', 'cp_air', 'Lv', 'epsilon');
         testCase.ro_atm = testCase.Pa / (Rd * testCase.Ta);
         testCase.cv_atm = testCase.ro_atm * cp_air;
         testCase.nu_air = 1.5e-5;
         testCase.hv_atm = testCase.ro_atm * Lv;

         z0_bulk = 1e-3;
         z_tair = 3;
         [De_h, testCase.br_coefs] = ...
            icemodel.surface.turbulence.bulk_richardson.exchange_coefficients( ...
            testCase.wspd, z0_bulk, z_tair);
         testCase.H_h = testCase.cv_atm * De_h;
         De_e = De_h * epsilon / testCase.Pa;
         testCase.H_e = testCase.hv_atm * De_e;
         testCase.ro_sfc = icemodel.physicalConstant('ro_ice');
         testCase.snow_depth = 0;

         testCase.dz = 0.04 * ones(500, 1);
         z = cumsum(testCase.dz);
         testCase.T = testCase.Ta * ones(500, 1);
         testCase.T = testCase.T .* exp(-0.1 * z ./ z(end));
      end
   end

   methods (Test)
      function testSolver(testCase, seb_solver)
         % Benchmark one surface-energy-balance root-finder variant on a fixed
         % surface state so only solver-path overhead changes between cases.
         batch_size = sebBenchmarkBatchSize(seb_solver);

         % Benchmark only the solve_surface_energy_balance call for each
         % available root finder.
         opts_sv = struct('seb_solver', seb_solver, 'debug', false, ...
            'turbulent_flux_scheme', 'bulk_richardson');

         [Ts, ok] = icemodel.surface.solve_surface_energy_balance( ...
            testCase.Ts0, testCase.Ta, testCase.Qsi, ...
            testCase.Qli, testCase.albedo, testCase.wspd, ...
            testCase.ppt, testCase.tppt, testCase.Pa, ...
            testCase.ea_atm, testCase.ro_atm, testCase.cv_atm, ...
            testCase.nu_air, testCase.H_h, testCase.H_e, ...
            testCase.hv_atm, testCase.br_coefs, ...
            testCase.liqflag, testCase.chi, testCase.T, ...
            testCase.k_eff, testCase.dz, testCase.ro_sfc, ...
            testCase.snow_depth, opts_sv);

         testCase.assertTrue(ok, ...
            'solve_surface_energy_balance setup did not converge before timing')
         testCase.assertThat(Ts, matlab.unittest.constraints.IsFinite)

         while testCase.keepMeasuring
            for n = 1:batch_size
               [Ts, ok] = icemodel.surface.solve_surface_energy_balance( ...
                  testCase.Ts0, testCase.Ta, testCase.Qsi, ...
                  testCase.Qli, testCase.albedo, testCase.wspd, ...
                  testCase.ppt, testCase.tppt, testCase.Pa, ...
                  testCase.ea_atm, testCase.ro_atm, ...
                  testCase.cv_atm, testCase.nu_air, ...
                  testCase.H_h, testCase.H_e, testCase.hv_atm, ...
                  testCase.br_coefs, testCase.liqflag, ...
                  testCase.chi, testCase.T, testCase.k_eff, ...
                  testCase.dz, testCase.ro_sfc, ...
                  testCase.snow_depth, opts_sv);
            end

            if ~ok || ~isfinite(Ts)
               error(['solve_surface_energy_balance benchmark failed for ' ...
                  'seb_solver=%d'], ...
                  seb_solver)
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
