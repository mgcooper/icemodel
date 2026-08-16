function tests = test_coupled_vapor_transport
   %TEST_COUPLED_VAPOR_TRANSPORT Verify the opt-in coupled vapor mass path.
   %
   % The coupled mode has two steps, and they conserve different things.
   %
   % icemodel.column.couple_vapor_transport moves vapor between cells by
   % Fick's law, which fixes the MASS. Both its boundaries are closed, so it
   % redistributes and creates nothing.
   % icemodel.column.apply_vapor_transport puts that mass into each
   % cell's phase, applying exactly the mass the cell received.
   %
   % The surface exchange is the other step and stays on its own path,
   % because the surface energy balance fixes the ENERGY and lets the mass
   % follow. Mixing the two is what loses mass: a cell told to spend an
   % energy demand can spend part of it at one latent heat and the rest at
   % another, and the mass it applies is then not the mass that arrived.
   %
   % See also: icemodel.column.couple_vapor_transport,
   %  icemodel.column.apply_vapor_transport
   tests = functiontests(localfunctions);
end

function test_isothermal_column_redistributes_nothing(testCase)
   % No interior vapor gradient means no interior flux. The surface exchange
   % is not part of this step, so an isothermal column moves nothing at all.
   % This is what lets the coupled mode reproduce the surface-only result:
   % the surface path is left to act alone, exactly as it does today.

   [dz, delz, fn, T, f_liq] = coupledFixture(8);
   [ro_vap, De] = icemodel.column.accepted_vapor_quantities(T, f_liq);

   [d_vap, dm_vap] = icemodel.column.couple_vapor_transport( ...
      ro_vap, De, dz, delz, fn, 900);

   testCase.verifyEqual(d_vap, zeros(8, 1), 'AbsTol', 1e-30);
   testCase.verifyEqual(dm_vap, zeros(8, 1), 'AbsTol', 1e-30);
end

function test_redistribution_conserves_column_mass(testCase)
   % Both boundaries are closed, so the transport moves mass between cells
   % and creates none. The column total must come back to zero.

   [dz, delz, fn, ~, f_liq] = coupledFixture(8);
   T = gradientColumn(8);
   [ro_vap, De] = icemodel.column.accepted_vapor_quantities(T, f_liq);

   [d_vap, dm_vap] = icemodel.column.couple_vapor_transport( ...
      ro_vap, De, dz, delz, fn, 900);

   % The cells must actually exchange, or the sums below are trivially zero.
   % Compare the net against the size of the exchange: the net is a
   % cancelling sum, so its absolute size means nothing on its own.
   testCase.verifyTrue(any(abs(d_vap) > 0));
   exchange = sum(abs(d_vap) .* dz);
   testCase.verifyGreaterThan(exchange, 0);
   testCase.verifyLessThan(abs(sum(d_vap .* dz)) / exchange, 1e-14);
   testCase.verifyLessThan( ...
      abs(sum(dm_vap .* dz)) / sum(abs(dm_vap) .* dz), 1e-14);
end

function test_applied_mass_equals_redistributed_mass_when_dry(testCase)
   % A dry column takes the exchange into ice. The mass the cells gained must
   % equal the mass the faces moved.

   [dz, delz, fn, ~, f_liq, f_ice, f_res_por] = coupledFixture(8);
   T = gradientColumn(8);
   [ro_vap, De] = icemodel.column.accepted_vapor_quantities(T, f_liq);

   d_vap = icemodel.column.couple_vapor_transport( ...
      ro_vap, De, dz, delz, fn, 900);
   [f_ice_new, f_liq_new, d_sbl_err] = ...
      icemodel.column.apply_vapor_transport( ...
      f_ice, f_liq, d_vap, 0.1, f_res_por);

   testCase.verifyEqual(d_sbl_err, zeros(8, 1), 'AbsTol', 0);
   testCase.verifyEqual(f_liq_new, f_liq, 'AbsTol', 0);
   verifyMassMatches(testCase, f_ice, f_liq, f_ice_new, f_liq_new, ...
      d_vap, dz);
end

function test_applied_mass_equals_redistributed_mass_when_mixed(testCase)
   % A column that mixes phases takes some exchange into liquid and some into
   % ice. The total mass must still match what the faces moved.

   [dz, delz, fn, ~, ~, f_ice, f_res_por] = coupledFixture(8);
   T = gradientColumn(8);
   f_liq = [0; 0; 0.05; 0.05; 0; 0; 0.05; 0];
   [ro_vap, De] = icemodel.column.accepted_vapor_quantities(T, f_liq);

   d_vap = icemodel.column.couple_vapor_transport( ...
      ro_vap, De, dz, delz, fn, 900);
   [f_ice_new, f_liq_new, d_sbl_err] = ...
      icemodel.column.apply_vapor_transport( ...
      f_ice, f_liq, d_vap, 0.1, f_res_por);

   testCase.verifyEqual(d_sbl_err, zeros(8, 1), 'AbsTol', 0);
   testCase.verifyTrue(any(f_liq_new ~= f_liq));
   testCase.verifyTrue(any(f_ice_new ~= f_ice));
   verifyMassMatches(testCase, f_ice, f_liq, f_ice_new, f_liq_new, ...
      d_vap, dz);
end

function test_a_thin_liquid_film_records_what_it_cannot_give(testCase)
   % A wet cell whose liquid margin is smaller than the mass leaving it once
   % lost about twelve percent of that mass: the energy-demand path let it
   % evaporate what it had and spend the rest at the other latent heat, with
   % nothing recording the difference. Applying mass cannot do that.
   %
   % A limited cell does break the column mass balance, because its
   % neighbours still receive what it could not supply. The DesignSpec
   % accepts that and requires the exception to be recorded rather than
   % stopping the run, so what this test pins is the record: the cell gives
   % up exactly what it has, and the shortfall appears in d_sbl_err.

   [dz, delz, fn, ~, ~, f_ice, f_res_por] = coupledFixture(5);

   % A temperature peak at cell 2 makes it the vapor-density maximum, so it
   % loses mass to both neighbours. A monotonic profile would move vapor
   % through that cell instead of out of it, and the film would never drain.
   Tf = icemodel.physicalConstant('Tf');
   T = (Tf - 8) * ones(5, 1);
   T(2) = Tf - 2;

   % Put that cell barely above its residual floor, so the loss exhausts the
   % film at once.
   f_res = icemodel.column.residual_water_fraction( ...
      f_ice, zeros(5, 1), f_res_por);
   f_liq = zeros(5, 1);
   f_liq(2) = f_res(2) + 1e-9;
   testCase.assertTrue(icemodel.column.vapor_exchange_is_wet( ...
      f_ice(2), f_liq(2), f_res_por));

   [ro_vap, De] = icemodel.column.accepted_vapor_quantities(T, f_liq);
   d_vap = icemodel.column.couple_vapor_transport( ...
      ro_vap, De, dz, delz, fn, 900);
   testCase.assertLessThan(d_vap(2), 0);

   [f_ice_new, f_liq_new, d_sbl_err] = ...
      icemodel.column.apply_vapor_transport( ...
      f_ice, f_liq, d_vap, 0.1, f_res_por);

   % The film gave up exactly what it had above the residual floor, and no
   % more. Nothing dipped into the ice of that cell behind the film.
   testCase.verifyEqual(f_liq_new(2), f_res(2), 'AbsTol', 1e-18);
   testCase.verifyEqual(f_ice_new(2), f_ice(2), 'AbsTol', 0);

   % The shortfall is recorded, so nothing goes missing without a trace.
   testCase.verifyLessThan(d_sbl_err(2), 0);
   testCase.verifyEqual(d_sbl_err([1, 3, 4, 5]), zeros(4, 1), 'AbsTol', 0);

   % Check the size of that record, not only its sign. What the column
   % applied plus what it recorded must equal what transport moved, or the
   % conversion from the shortfall mass to d_sbl_err is scaled wrong and a
   % sign check would still pass.
   verifyMassAndRecordMatch(testCase, f_ice, f_liq, f_ice_new, f_liq_new, ...
      d_sbl_err, d_vap, dz, f_res_por);

   % Every unlimited cell still applied exactly the mass it received.
   verifyCellMassMatches(testCase, f_ice, f_liq, f_ice_new, f_liq_new, ...
      d_vap, [1, 3, 4, 5]);
end

function test_max_liquid_change_converts_air_to_liquid_equivalent(testCase)
   % One function owns the condensation limit. The pore volume converts to a
   % liquid-water-equivalent fraction as if it were ice, so a cell that later
   % refreezes the water it took cannot exceed a full control volume.

   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   f_ice = 0.6;
   f_liq = 0.05;

   returned = icemodel.column.max_liquid_fraction_change(f_ice, f_liq);
   expected = ro_ice / ro_liq * (1.0 - f_ice) - f_liq;
   testCase.verifyEqual(returned, expected, 'RelTol', 1e-15);

   % Filling the capacity and refreezing it lands exactly at a full cell.
   f_ice_refrozen = f_ice + (f_liq + returned) * ro_liq / ro_ice;
   testCase.verifyEqual(f_ice_refrozen, 1.0, 'AbsTol', 1e-14);

   % A saturated cell has no room, and an over-full one reports negative so
   % the callers' max(capacity, 0) floor is the thing that clamps.
   testCase.verifyEqual( ...
      icemodel.column.max_liquid_fraction_change(1.0, 0.0), 0.0, 'AbsTol', 0);
   testCase.verifyLessThan( ...
      icemodel.column.max_liquid_fraction_change(0.9, 0.5), 0);

   % It maps over a column, which is how the interior applier uses it.
   f_ice_col = [0.3; 0.6; 0.9];
   f_liq_col = [0.0; 0.05; 0.01];
   testCase.verifyEqual( ...
      icemodel.column.max_liquid_fraction_change(f_ice_col, f_liq_col), ...
      ro_ice / ro_liq * (1.0 - f_ice_col) - f_liq_col, 'RelTol', 1e-15);
end

function test_a_bound_interior_clamp_leaves_a_closure_residual(testCase)
   % Known defect, bead icemodel-55x. icemodel.column.apply_vapor_transport
   % records a clamped cell's shortfall in d_sbl_err, which the ledger spends
   % in its unapplied term. icemodel.column.accumulate_redistribution_budget
   % records only the storage change that survived the clamp. With no surface
   % exchange the identity reduces to
   %
   %   0 = storage + unapplied - redistribution
   %
   % and the shortfall never cancels, so the identity carries a residual the
   % size of the shortfall.
   %
   % This test pins the defect so it cannot be lost. Invert it when the fix
   % lands: the residual must then be zero to roundoff.

   [dz, delz, fn, ~, ~, f_ice, f_res_por] = coupledFixture(5);
   [Ls, Lv, ro_liq, ro_ice] = icemodel.physicalConstant( ...
      'Ls', 'Lv', 'ro_liq', 'ro_ice');
   Tf = icemodel.physicalConstant('Tf');

   % A thin film at cell 2 makes the evaporation clamp bind.
   T = (Tf - 8) * ones(5, 1);
   T(2) = Tf - 2;
   f_res = icemodel.column.residual_water_fraction( ...
      f_ice, zeros(5, 1), f_res_por);
   f_liq = zeros(5, 1);
   f_liq(2) = f_res(2) + 1e-9;

   ledger = icemodel.column.initialize_budget_state();
   [f_ice_new, f_liq_new, d_sbl_err, ledger_new] = ...
      icemodel.column.couple_vapor_step( ...
      T, f_ice, f_liq, zeros(5, 1), dz, delz, fn, 900, 0.1, ...
      f_res_por, ledger, true);

   % The clamp bound, or the case below is vacuous.
   testCase.assertLessThan(d_sbl_err(2), 0);

   [solid_0, liquid_0] = icemodel.column.integrate_column_budget( ...
      T, f_ice, f_liq, dz);
   [solid_1, liquid_1] = icemodel.column.integrate_column_budget( ...
      T, f_ice_new, f_liq_new, dz);
   storage = ro_liq * (Ls * (solid_1 - solid_0) ...
      + Lv * (liquid_1 - liquid_0));
   unapplied = sum(ro_ice * Ls * d_sbl_err(:) .* dz(:));
   residual = storage + unapplied ...
      - ledger_new.mass_budget_vapor_redistribution_j_m2;

   % The residual is the shortfall, not zero.
   testCase.verifyEqual(residual, unapplied, 'RelTol', 1e-12);
   testCase.verifyGreaterThan(abs(residual), 1);
end

function test_the_orchestrator_matches_its_parts(testCase)
   % icemodel.column.couple_vapor_step keeps the driver free of vapor
   % intermediates. It must produce exactly what the same primitives produce
   % when called in order, or the orchestration itself changed the physics.

   [dz, delz, fn, ~, ~, f_ice, f_res_por] = coupledFixture(6);
   T = gradientColumn(6);
   f_liq = [0; 0; 0.05; 0.05; 0; 0];
   dt = 900;
   f_ice_min = 0.1;
   ledger = icemodel.column.initialize_budget_state();

   % The orchestrated path.
   d_sbl_err_in = zeros(6, 1);
   [f_ice_a, f_liq_a, d_sbl_err_a, ledger_a] = ...
      icemodel.column.couple_vapor_step( ...
      T, f_ice, f_liq, d_sbl_err_in, dz, delz, fn, dt, f_ice_min, ...
      f_res_por, ledger, true);

   % The same primitives, called directly in the driver's order.
   [ro_vap, De] = icemodel.column.accepted_vapor_quantities(T, f_liq);
   d_vap = icemodel.column.couple_vapor_transport( ...
      ro_vap, De, dz, delz, fn, dt);
   [solid_r, liquid_r] = icemodel.column.integrate_column_budget( ...
      T, f_ice, f_liq, dz);
   [f_ice_b, f_liq_b, d_sbl_err_b] = ...
      icemodel.column.apply_vapor_transport( ...
      f_ice, f_liq, d_vap, f_ice_min, f_res_por);
   ledger_b = icemodel.column.accumulate_redistribution_budget( ...
      ledger, solid_r, liquid_r, T, f_ice_b, f_liq_b, dz);

   testCase.verifyEqual(f_ice_a, f_ice_b, 'AbsTol', 0);
   testCase.verifyEqual(f_liq_a, f_liq_b, 'AbsTol', 0);
   testCase.verifyEqual(d_sbl_err_a, d_sbl_err_b, 'AbsTol', 0);
   testCase.verifyEqual( ...
      ledger_a.mass_budget_vapor_redistribution_j_m2, ...
      ledger_b.mass_budget_vapor_redistribution_j_m2, 'AbsTol', 0);

   % The column actually exchanged, or the comparison above is trivial.
   testCase.verifyTrue(any(f_ice_a ~= f_ice));
end

function test_the_orchestrator_skips_the_ledger_when_it_is_off(testCase)
   % A production run without the diagnostic profile passes
   % use_mass_budget false. Nothing in the full-model suite reaches that
   % branch, because every coupled run there is a diagnostic run. The state
   % must come out identical and the ledger must stay untouched.

   [dz, delz, fn, ~, ~, f_ice, f_res_por] = coupledFixture(6);
   T = gradientColumn(6);
   f_liq = [0; 0; 0.05; 0.05; 0; 0];
   dt = 900;
   f_ice_min = 0.1;
   ledger = icemodel.column.initialize_budget_state();
   d_sbl_err_in = zeros(6, 1);

   [f_ice_on, f_liq_on, err_on, ledger_on] = ...
      icemodel.column.couple_vapor_step( ...
      T, f_ice, f_liq, d_sbl_err_in, dz, delz, fn, dt, f_ice_min, ...
      f_res_por, ledger, true);
   [f_ice_off, f_liq_off, err_off, ledger_off] = ...
      icemodel.column.couple_vapor_step( ...
      T, f_ice, f_liq, d_sbl_err_in, dz, delz, fn, dt, f_ice_min, ...
      f_res_por, ledger, false);

   % The physics does not depend on whether the ledger is being built.
   testCase.verifyEqual(f_ice_off, f_ice_on, 'AbsTol', 0);
   testCase.verifyEqual(f_liq_off, f_liq_on, 'AbsTol', 0);
   testCase.verifyEqual(err_off, err_on, 'AbsTol', 0);

   % The ledger comes back as it went in, and the on-run moved it, so the
   % check is not passing because nothing happened.
   testCase.verifyEqual(ledger_off, ledger, 'AbsTol', 0);
   testCase.verifyNotEqual( ...
      ledger_on.mass_budget_vapor_redistribution_j_m2, ...
      ledger.mass_budget_vapor_redistribution_j_m2);
end

function test_the_orchestrator_accumulates_into_the_incoming_record(testCase)
   % d_sbl_err threads in and out so the increment accumulates like every
   % other d_* increment. A non-zero record arriving must survive.

   [dz, delz, fn, ~, ~, f_ice, f_res_por] = coupledFixture(5);
   Tf = icemodel.physicalConstant('Tf');
   T = (Tf - 8) * ones(5, 1);
   T(2) = Tf - 2;
   f_res = icemodel.column.residual_water_fraction( ...
      f_ice, zeros(5, 1), f_res_por);
   f_liq = zeros(5, 1);
   f_liq(2) = f_res(2) + 1e-9;
   ledger = icemodel.column.initialize_budget_state();

   [~, ~, from_zero] = icemodel.column.couple_vapor_step( ...
      T, f_ice, f_liq, zeros(5, 1), dz, delz, fn, 900, 0.1, ...
      f_res_por, ledger, true);
   testCase.assertTrue(any(from_zero ~= 0));

   % The same call with a record already in it must return the sum.
   prior = [1e-7; 0; -2e-7; 0; 0];
   [~, ~, from_prior] = icemodel.column.couple_vapor_step( ...
      T, f_ice, f_liq, prior, dz, delz, fn, 900, 0.1, ...
      f_res_por, ledger, true);

   testCase.verifyEqual(from_prior, prior + from_zero, 'AbsTol', 0);
end

function verifyCellMassMatches(testCase, f_ice, f_liq, f_ice_new, ...
      f_liq_new, d_vap, cells)
   %VERIFYCELLMASSMATCHES Compare per-cell applied mass against its share.

   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   for j = cells
      applied = ro_ice * (f_ice_new(j) - f_ice(j)) ...
         + ro_liq * (f_liq_new(j) - f_liq(j));
      moved = ro_liq * d_vap(j);
      testCase.verifyEqual(applied, moved, ...
         'RelTol', 1e-12, 'AbsTol', 1e-18);
   end
end

function test_mass_conserves_in_the_phase_decision_band(testCase)
   % Two ways of asking "is this cell wet" exist in the model, on different
   % thresholds, and they disagree in a band that always exists. A cell in
   % that band once gained thirteen percent or lost twelve percent of its
   % mass. Both sides ask one function.

   [dz, delz, fn] = coupledFixture(3);
   T = gradientColumn(3);
   f_res_por = 0.07;

   % Dense ice: the residual floor sits below the fixed threshold.
   verifyBandConserves(testCase, 0.90 * ones(3, 1), [0; 0.0135; 0], ...
      f_res_por, T, dz, delz, fn, true);

   % Low-density snow: the floor sits above it, so the disagreement flips.
   verifyBandConserves(testCase, 0.30 * ones(3, 1), [0; 0.03; 0], ...
      f_res_por, T, dz, delz, fn, false);
end

function verifyBandConserves(testCase, f_ice, f_liq, f_res_por, T, dz, ...
      delz, fn, expect_wet)
   %VERIFYBANDCONSERVES Check one column's applied mass against its transport.

   testCase.verifyEqual(icemodel.column.vapor_exchange_is_wet( ...
      f_ice(2), f_liq(2), f_res_por), expect_wet);

   [ro_vap, De] = icemodel.column.accepted_vapor_quantities(T, f_liq);
   d_vap = icemodel.column.couple_vapor_transport( ...
      ro_vap, De, dz, delz, fn, 900);
   [f_ice_new, f_liq_new, d_sbl_err] = ...
      icemodel.column.apply_vapor_transport( ...
      f_ice, f_liq, d_vap, 0.1, f_res_por);

   testCase.verifyEqual(d_sbl_err, zeros(numel(f_ice), 1), 'AbsTol', 0);
   verifyMassMatches(testCase, f_ice, f_liq, f_ice_new, f_liq_new, ...
      d_vap, dz);
end

function verifyMassMatches(testCase, f_ice, f_liq, f_ice_new, f_liq_new, ...
      d_vap, dz)
   %VERIFYMASSMATCHES Compare applied mass against redistributed mass.

   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   applied = ro_ice * sum((f_ice_new - f_ice) .* dz) ...
      + ro_liq * sum((f_liq_new - f_liq) .* dz);
   moved = ro_liq * sum(d_vap .* dz);

   % Both are near zero over a closed column, so compare against the size of
   % the exchange rather than against the net.
   scale = ro_liq * sum(abs(d_vap) .* dz);
   testCase.verifyGreaterThan(scale, 0);
   testCase.verifyLessThan(abs(applied - moved), 1e-9 * scale);
end

function test_wet_condensation_caps_at_the_pore_capacity(testCase)
   % A wet cell can only hold so much water. Transport that delivers more
   % than the pore space takes must apply the capacity and record the rest,
   % or the column gains water the pores cannot hold.
   %
   % The applier is called directly here. Driving this branch through
   % couple_vapor_transport would need a temperature field tuned to deliver a
   % specific mass, which tests the fixture rather than the branch.

   f_ice = 0.99;
   f_liq = 0.005;
   f_res_por = 0.02;
   testCase.assertTrue(icemodel.column.vapor_exchange_is_wet( ...
      f_ice, f_liq, f_res_por));

   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   capacity = ro_ice / ro_liq * (1.0 - f_ice) - f_liq;
   testCase.assertGreaterThan(capacity, 0);

   % Deliver twice what the cell can hold.
   d_vap = 2 * capacity;
   [f_ice_new, f_liq_new, d_sbl_err] = ...
      icemodel.column.apply_vapor_transport( ...
      f_ice, f_liq, d_vap, 0.1, f_res_por);

   testCase.verifyEqual(f_liq_new, f_liq + capacity, 'RelTol', 1e-14);
   testCase.verifyEqual(f_ice_new, f_ice, 'AbsTol', 0);

   % The unapplied half is recorded, on the energy basis the ledger spends.
   expected = icemodel.column.potential_sublimation(capacity);
   testCase.verifyEqual(d_sbl_err, expected, 'RelTol', 1e-14);
end

function test_dry_deposition_caps_at_the_air_space(testCase)
   % A dry cell with almost no air space cannot take the vapor arriving at
   % it. The applier must fill the space it has and record the remainder.

   f_ice = 0.999;
   f_liq = 0.0;
   f_res_por = 0.02;
   testCase.assertFalse(icemodel.column.vapor_exchange_is_wet( ...
      f_ice, f_liq, f_res_por));

   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   f_air = 1.0 - f_ice - f_liq;

   % Deliver twice the air space, as the ice fraction the applier converts to.
   d_vap = 2 * f_air * ro_ice / ro_liq;
   [f_ice_new, f_liq_new, d_sbl_err] = ...
      icemodel.column.apply_vapor_transport( ...
      f_ice, f_liq, d_vap, 0.1, f_res_por);

   testCase.verifyEqual(f_ice_new, f_ice + f_air, 'RelTol', 1e-14);
   testCase.verifyEqual(f_liq_new, f_liq, 'AbsTol', 0);

   % A dry cell records the ice fraction of the mass it could not take, which
   % is the leftover air space it had no room for.
   testCase.verifyEqual(d_sbl_err, f_air, 'RelTol', 1e-12);
end

function test_dry_sublimation_floors_at_the_retained_ice(testCase)
   % A cell at the ice floor has nothing left to give. It must stop at the
   % floor rather than go through it, and record what it could not supply.

   f_ice_min = 0.1;
   f_ice = f_ice_min + 1e-6;
   f_liq = 0.0;
   f_res_por = 0.02;
   testCase.assertFalse(icemodel.column.vapor_exchange_is_wet( ...
      f_ice, f_liq, f_res_por));

   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   available = f_ice - f_ice_min;

   % Ask for far more than the cell holds above the floor.
   d_vap = -100 * available * ro_ice / ro_liq;
   [f_ice_new, f_liq_new, d_sbl_err] = ...
      icemodel.column.apply_vapor_transport( ...
      f_ice, f_liq, d_vap, f_ice_min, f_res_por);

   testCase.verifyEqual(f_ice_new, f_ice_min, 'AbsTol', 1e-18);
   testCase.verifyEqual(f_liq_new, f_liq, 'AbsTol', 0);

   % The shortfall is negative and carries the mass the cell never supplied.
   d_ice = d_vap * ro_liq / ro_ice;
   testCase.verifyEqual(d_sbl_err, d_ice + available, 'RelTol', 1e-12);
   testCase.verifyLessThan(d_sbl_err, 0);
end

function verifyMassAndRecordMatch(testCase, f_ice, f_liq, f_ice_new, ...
      f_liq_new, d_sbl_err, d_vap, dz, f_res_por)
   %VERIFYMASSANDRECORDMATCH Applied plus recorded must equal moved.
   %
   % d_sbl_err is an energy-basis channel: the ledger spends it as
   % ro_ice * Ls * d_sbl_err. A dry cell records the ice fraction of the mass
   % it could not take, and Ls converts that same fraction to the right
   % energy. A wet cell records the ice fraction of the same latent ENERGY,
   % which is the unapplied liquid scaled by Lv / Ls, so it is not the ice
   % fraction of that mass. Recovering a mass from the channel therefore has
   % to undo whichever conversion the cell used.
   %
   % Asking icemodel.column.vapor_exchange_is_wet rather than testing a
   % threshold here keeps this inversion tied to the decision the applier
   % made. A second criterion would put a cell in the band the applier's
   % single owner exists to eliminate.

   [ro_ice, ro_liq, Ls, Lv] = icemodel.physicalConstant( ...
      'ro_ice', 'ro_liq', 'Ls', 'Lv');

   applied = ro_ice * sum((f_ice_new - f_ice) .* dz) ...
      + ro_liq * sum((f_liq_new - f_liq) .* dz);
   moved = ro_liq * sum(d_vap .* dz);

   % Invert the conversion per cell, on the state the applier saw.
   wet = icemodel.column.vapor_exchange_is_wet(f_ice, f_liq, f_res_por);
   shortfall = ro_ice * d_sbl_err;
   shortfall(wet) = shortfall(wet) * Ls / Lv;
   recorded = sum(shortfall .* dz);

   scale = ro_liq * sum(abs(d_vap) .* dz);
   testCase.verifyGreaterThan(scale, 0);
   testCase.verifyLessThan(abs(applied + recorded - moved), 1e-9 * scale);
end

function T = gradientColumn(JJ)
   %GRADIENTCOLUMN Return a temperature profile that drives interior vapor.

   Tf = icemodel.physicalConstant('Tf');
   T = (Tf - 8) + linspace(0, 6, JJ)';
end

function [dz, delz, fn, T, f_liq, f_ice, f_res_por] = coupledFixture(JJ)
   %COUPLEDFIXTURE Return one dry, isothermal column on the production mesh.

   Tf = icemodel.physicalConstant('Tf');
   [dz, delz, ~, ~, fn] = ...
      icemodel.column.control_volume_mesh(JJ * 0.04, 0.04);
   dz = dz(1:JJ);
   delz = delz(1:JJ + 1);
   fn = fn(1:JJ + 1);
   T = (Tf - 5) * ones(JJ, 1);
   f_liq = zeros(JJ, 1);
   f_ice = 0.85 * ones(JJ, 1);
   f_res_por = 0.02;
end
