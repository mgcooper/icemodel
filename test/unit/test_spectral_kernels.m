function tests = test_spectral_kernels
   %TEST_SPECTRAL_KERNELS Verify spectral and radiative-transfer kernels.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Build one spectral-enabled column state so the spectral helpers all see the
   % same controlled geometry and forcing setup.

   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=96, dt_seconds=900);

   testCase.TestData.workspace = workspace;

   testCase.TestData.state = icemodel.test.fixtures.makeSyntheticColumnState( ...
      workspace, 'icemodel', solver=3, metstep=49, include_spectral=true, ...
      testname='spectral_kernel');
end

function teardown(testCase)
   % Remove the shared spectral workspace after the file-level tests end.

   icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      testCase.TestData.workspace);
end

function test_terrain_adjusted_shortwave_distinguishes_day_and_night(testCase)
   % terrain_adjusted_shortwave_radiation should produce zero at night while
   % remaining positive during daytime on the same site geometry.

   Qday = icemodel.surface.terrain_adjusted_shortwave_radiation( ...
      180, 45, 0.2, 12, 180, 5, 0.7);
   Qnight = icemodel.surface.terrain_adjusted_shortwave_radiation( ...
      180, 45, 0.2, 0, 180, 5, 0.7);

   testCase.verifyGreaterThan(Qday, 0);
   testCase.verifyEqual(Qnight, 0, 'AbsTol', 1e-12);
end

function test_incoming_shortwave_daily_average_matches_hourly_mean(testCase)
   % The daily shortwave fallback should equal the mean of the hourly
   % terrain-adjusted samples it wraps.

   J_day_start = 180;
   step = 1;
   dt = 86400;
   xlat = 45;
   cloud_frac = 0.2;
   slope_az = 180;
   terrain_slope = 5;
   ihrs_day = 24;
   transmiss = 0.7;
   Qhourly = zeros(ihrs_day, 1);

   for hour_sample = 1:ihrs_day
      Qhourly(hour_sample) = ...
         icemodel.surface.terrain_adjusted_shortwave_radiation( ...
         J_day_start, xlat, cloud_frac, hour_sample, slope_az, ...
         terrain_slope, transmiss);
   end

   Qdaily = icemodel.surface.incoming_shortwave_radiation( ...
      J_day_start, step, dt, xlat, cloud_frac, slope_az, terrain_slope, ...
      ihrs_day, transmiss, 12);

   testCase.verifyEqual(Qdaily, mean(Qhourly), 'RelTol', 1e-12);
end

function test_getscattercoefs_and_getsolar_return_positive_weights(testCase)
   % The spectral wavel weights and integrated prototype spectrum should remain
   % positive on a simple monotonic wavel grid.

   s = testCase.TestData.state;
   solar_data = load(fullfile(s.opts.pathinput, 'spectral', 'solar.mat'));

   [qext, g, coalbedo, wavel, dwavel, radii] ...
      = icemodel.radiation.get_scattering_coefficients(s.opts);

   [I0, solar] ...
      = icemodel.radiation.get_solar_spectrum(solar_data.solar, wavel, dwavel);

   testCase.verifyEqual(size(qext), [s.opts.nradii, s.opts.nwavel]);
   testCase.verifyEqual(size(g), [s.opts.nradii, s.opts.nwavel]);
   testCase.verifyEqual(size(coalbedo), [s.opts.nradii, s.opts.nwavel]);
   testCase.verifyEqual(numel(wavel), s.opts.nwavel);
   testCase.verifyGreaterThan(min(dwavel), 0);
   testCase.verifyEqual(numel(radii), s.opts.nradii);
   testCase.verifyEqual(numel(solar), s.opts.nwavel);
   testCase.verifyGreaterThan(I0, 0);
end

function test_getupdown_returns_finite_fluxes(testCase)
   % smooth_twostream_fluxes should reconstruct finite up/down fluxes on a
   % compact hand-built coefficient profile.

   bulkcoefs = [1.0; 1.2; 1.4; 1.4; 1.4];
   albedo = 0.5;
   a = ((1.0 - albedo) / (1.0 + albedo)) * bulkcoefs;
   r = (2.0 * albedo / (1.0 - albedo ^ 2)) * bulkcoefs;
   x = [5; 4; 3; 2];
   z_edges = [0; 0.1; 0.2; 0.3; 0.4];
   [up, down] = icemodel.radiation.smooth_twostream_fluxes(a, r, x, 100, z_edges, 3);

   testCase.verifyEqual(numel(up), 5);
   testCase.verifyEqual(numel(down), 5);
   testCase.verifyTrue(all(isfinite(up)));
   testCase.verifyTrue(all(isfinite(down)));
end

function test_solvetwostream_returns_finite_net_flux_profile(testCase)
   % icemodel.radiation.solvetwostream should return a finite net-flux profile
   % with the top boundary reflecting the imposed incoming radiation.

   bulkcoefs = [1.0; 1.2; 1.4; 1.4; 1.4; 1.4];
   z_edges = [0; 0.1; 0.2; 0.3; 0.4];
   xynet = icemodel.radiation.solvetwostream(100, 0.5, bulkcoefs, z_edges);

   testCase.verifyEqual(numel(xynet), numel(z_edges));
   testCase.verifyTrue(all(isfinite(xynet)));
   testCase.verifyLessThanOrEqual(xynet(1), -50 + 1e-9);
end

function test_extcoefsinit_and_spectralsourceterm_return_finite_terms(testCase)
   % icemodel.column.shortwave_source_term should keep the spectral source term
   % and chi finite on the shared synthetic spectral state.

   s = testCase.TestData.state;

   ro_sno = s.ro_ice * s.f_ice + s.ro_liq * s.f_liq + ...
      s.ro_air * (1 - s.f_ice - s.f_liq);

   % Use an empty lookup table to follow the exact k_bulk path
   k_lookup_empty = struct([]);

   [Sc_new, chi] = icemodel.column.shortwave_source_term( ...
      s.swd, s.albedo, s.I0, s.dz_spect, s.tau_N, s.tau_S, ...
      s.solar_dwavel, s.dz, ro_sno, s.z_nodes, s.z_nodes_spect, ...
      s.z_edges_spect, k_lookup_empty);

   testCase.verifyEqual(size(Sc_new), size(s.Sc));
   testCase.verifyTrue(all(isfinite(Sc_new)));
   testCase.verifyGreaterThan(chi, 0);
   testCase.verifyLessThanOrEqual(chi, 1 + 1e-9);
end

function test_updateextcoefs_rebuilds_current_k_ext_and_tau(testCase)
   % update_extinction_coefficients should reproduce the initialized
   % k_ext/tau state for the configured optical grain-radius index.

   s = testCase.TestData.state;

   [tau_N, tau_S, ~, k_ext] ...
      = icemodel.radiation.update_extinction_coefficients(s.qext, s.g, s.coalbedo, ...
      s.kabs, s.kice, s.wavel, s.radii, s.opts.i_grainradius, ...
      s.z_edges_spect, s.dz_spect, s.solar_dwavel, s.ro_ice, false);

   testCase.verifyEqual(k_ext, s.k_ext, 'AbsTol', 1e-12);
   testCase.verifyEqual(tau_N, s.tau_N, 'AbsTol', 1e-12);
   testCase.verifyEqual(tau_S, s.tau_S, 'AbsTol', 1e-12);
end

function test_spectextcoefs_interpolates_fractional_radius_index(testCase)
   % spectral_extinction_coefficients should linearly interpolate the optical
   % tables when a future grain-size model maps to a fractional table index.

   s = testCase.TestData.state;

   iradius = 10.5;

   k_ext = icemodel.radiation.spectral_extinction_coefficients( ...
      s.qext, s.g, s.coalbedo, s.radii, iradius);

   i0 = floor(iradius);
   i1 = ceil(iradius);
   w1 = iradius - i0;
   w0 = 1.0 - w1;

   r_snow = (w0 * s.radii(i0) + w1 * s.radii(i1)) / 1000.0;
   qext = w0 * s.qext(i0, :) + w1 * s.qext(i1, :);
   g = w0 * s.g(i0, :) + w1 * s.g(i1, :);
   coalbedo = w0 * s.coalbedo(i0, :) + w1 * s.coalbedo(i1, :);
   sigma_e = (3.0 / 4.0) * qext / r_snow;

   k_ext_expected = sigma_e .* sqrt(coalbedo - coalbedo .* g ...
      + coalbedo .^ 2 .* g);

   testCase.verifyTrue(all(isfinite(k_ext)));
   testCase.verifyEqual(k_ext, k_ext_expected, 'AbsTol', 1e-12);
end

function test_bulkextcoefs_matches_inline_transform(testCase)
   % icemodel.radiation.bulk_extinction_coefficients should reproduce the exact
   % spectral extinction transform embedded in the historical inlined spectral
   % source-term implementation.

   s = testCase.TestData.state;

   ro_sno = s.ro_ice * s.f_ice + s.ro_liq * s.f_liq + ...
      s.ro_air * (1 - s.f_ice - s.f_liq);

   ro_sno_spect = max(interp1(s.z_nodes, ro_sno, s.z_nodes_spect, ...
      'nearest', 'extrap'), 300);

   legacy = -log((sum(s.solar_dwavel .* exp( ...
      s.tau_S .* ro_sno_spect), 2)) ./ ...
      (sum(s.solar_dwavel .* exp( ...
      s.tau_N .* ro_sno_spect), 2))) / s.dz_spect(1);
   legacy = [legacy; legacy(end); legacy(end)];

   helper = icemodel.radiation.bulk_extinction_coefficients( ...
      s.dz_spect, ro_sno_spect, s.tau_N, s.tau_S, s.solar_dwavel);

   testCase.verifyEqual(helper, legacy, 'AbsTol', 1e-12);
end

function test_spectralsourceterm_matches_inlined_path(testCase)
   % The exact functions path should match the historical inlined path so the
   % perf study isolates code organization and lookup effects, not behavior.

   s = testCase.TestData.state;
   ro_sno = s.ro_ice * s.f_ice + s.ro_liq * s.f_liq + ...
      s.ro_air * (1 - s.f_ice - s.f_liq);

   % Use an empty lookup table to follow the exact k_bulk path
   k_lookup_empty = struct([]);

   [Sc_inlined, chi_inlined] = SPECTRALSOURCETERM_INLINE(s.swd, s.albedo, ...
      s.I0, s.dz_spect, s.tau_N, s.tau_S, s.solar_dwavel, s.dz, ro_sno, ...
      s.z_nodes, s.z_nodes_spect);

   [Sc_functions, chi_functions] ...
      = icemodel.column.shortwave_source_term(s.swd, s.albedo, ...
      s.I0, s.dz_spect, s.tau_N, s.tau_S, s.solar_dwavel, s.dz, ...
      ro_sno, s.z_nodes, s.z_nodes_spect, s.z_edges_spect, k_lookup_empty);

   testCase.verifyEqual(Sc_functions, Sc_inlined, 'RelTol', 1e-12, ...
      'AbsTol', 1e-12);
   testCase.verifyEqual(chi_functions, chi_inlined, 'AbsTol', 1e-12);
end

function test_bulkextcoefslookup_stays_close_to_exact_transform(testCase)
   % icemodel.radiation.bulk_extinction_coefficients_lookup is approximate, but
   % on the integer-density lookup grid it should stay close to the exact
   % bulk-extinction profile.

   s = testCase.TestData.state;

   ro_sno = s.ro_ice * s.f_ice + s.ro_liq * s.f_liq + ...
      s.ro_air * (1 - s.f_ice - s.f_liq);

   ro_sno_spect = max(interp1(s.z_nodes, ro_sno, s.z_nodes_spect, ...
      'nearest', 'extrap'), 300);

   lookup = icemodel.makeBulkExtCoefsLookup(s.dz_spect, ...
      s.tau_N, s.tau_S, s.solar_dwavel);

   exact = icemodel.radiation.bulk_extinction_coefficients( ...
      s.dz_spect, ro_sno_spect, s.tau_N, s.tau_S, s.solar_dwavel);

   approx = icemodel.radiation.bulk_extinction_coefficients_lookup( ...
      ro_sno_spect, lookup);

   rel_err = max(abs(approx - exact) ./ max(abs(exact), 1e-12));

   testCase.verifyTrue(all(isfinite(approx)));
   testCase.verifyLessThan(rel_err, 0.02);
end
