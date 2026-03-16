function tests = test_spectral_kernels
%TEST_SPECTRAL_KERNELS Verify spectral and radiative-transfer kernels.
   tests = functiontests(localfunctions);
end

function setup(testCase)

   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=24, dt_seconds=900);
   testCase.TestData.workspace = workspace;
   testCase.TestData.state = icemodel.test.fixtures.makeSyntheticColumnState( ...
      workspace, 'icemodel', solver=3, include_spectral=true, ...
      testname='spectral_kernel');
end

function teardown(testCase)

   icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      testCase.TestData.workspace);
end

function test_solarrad_distinguishes_day_and_night(testCase)

   Qday = SOLARRAD(180, 45, 0.2, 12, 180, 5, 0.7);
   Qnight = SOLARRAD(180, 45, 0.2, 0, 180, 5, 0.7);

   testCase.verifyGreaterThan(Qday, 0);
   testCase.verifyEqual(Qnight, 0, 'AbsTol', 1e-12);
end

function test_getdwavl_and_getsolar_return_positive_weights(testCase)

   wavelength = [0.4 0.6 1.0 1.4];
   dwavl = GETDWAVL(wavelength, numel(wavelength));
   solar_in = [0.3 10; 0.8 20; 1.5 30];
   [solar, Q0] = GETSOLAR(solar_in, numel(wavelength), wavelength, dwavl);

   testCase.verifyGreaterThan(min(dwavl), 0);
   testCase.verifyEqual(numel(solar), numel(wavelength));
   testCase.verifyGreaterThan(Q0, 0);
end

function test_specinit_and_getscattercoefs_return_expected_shapes(testCase)

   s = testCase.TestData.state;
   [~, mie, solar, kabs, kice] = SPECINIT(s.opts);
   [g, qext, ss_coalb, wavelength] = GETSCATTERCOEFS(s.opts, mie);

   testCase.verifyEqual(size(g), [s.opts.nradii, s.opts.nwavl]);
   testCase.verifyEqual(size(qext), [s.opts.nradii, s.opts.nwavl]);
   testCase.verifyEqual(size(ss_coalb), [s.opts.nradii, s.opts.nwavl]);
   testCase.verifyEqual(size(wavelength), [s.opts.nradii, s.opts.nwavl]);
   testCase.verifyFalse(isempty(solar));
   if s.opts.kabs_user
      testCase.verifyFalse(isempty(kabs));
      testCase.verifyFalse(isempty(kice));
   end
end

function test_getaandr_and_getupdown_return_finite_fluxes(testCase)

   bulkcoefs = [1.0; 1.2; 1.4; 1.4; 1.4];
   [a, r] = GETAANDR(bulkcoefs, 0.5);
   x = [5; 4; 3; 2];
   z_walls = [0; 0.1; 0.2; 0.3; 0.4];
   [up, down] = GETUPDOWN(a, r, x, 100, z_walls, 3);

   testCase.verifyEqual(numel(up), 5);
   testCase.verifyEqual(numel(down), 5);
   testCase.verifyTrue(all(isfinite(up)));
   testCase.verifyTrue(all(isfinite(down)));
end

function test_solvetwostream_returns_finite_net_flux_profile(testCase)

   bulkcoefs = [1.0; 1.2; 1.4; 1.4; 1.4; 1.4];
   [a, r] = GETAANDR(bulkcoefs, 0.5);
   z_walls = [0; 0.1; 0.2; 0.3; 0.4];
   xynet = SOLVETWOSTREAM(a, r, bulkcoefs, 100, 0.5, z_walls);

   testCase.verifyEqual(numel(xynet), 4);
   testCase.verifyTrue(all(isfinite(xynet)));
   testCase.verifyLessThanOrEqual(xynet(1), -50 + 1e-9);
end

function test_extcoefsinit_and_updateextcoefs_return_finite_terms(testCase)

   s = testCase.TestData.state;
   ro_sno = s.ro_ice * s.f_ice + s.ro_liq * s.f_liq + ...
      s.ro_air * (1 - s.f_ice - s.f_liq);

   [Sc_new, chi] = UPDATEEXTCOEFS(s.swd, s.albedo, s.Q0, s.dz_spect, ...
      s.spect_N, s.spect_S, s.solardwavl, s.Sc, s.dz, ro_sno);

   testCase.verifyEqual(size(Sc_new), size(s.Sc));
   testCase.verifyTrue(all(isfinite(Sc_new)));
   testCase.verifyGreaterThan(chi, 0);
   testCase.verifyLessThanOrEqual(chi, 1 + 1e-9);
end
