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
      applyTransport(f_ice, f_liq, d_vap, 0.1, f_res_por);

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
      applyTransport(f_ice, f_liq, d_vap, 0.1, f_res_por);

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
      applyTransport(f_ice, f_liq, d_vap, 0.1, f_res_por);

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

function test_a_bound_interior_clamp_keeps_its_own_accounting(testCase)
   % The transport keeps its own accounting when a per-cell limit binds:
   % the redistribution channels equal the realized per-phase storage
   % moves exactly, the shortfall lands in the redistribution's own
   % unapplied channel, and the surface unapplied channel never sees an
   % interior term. That separation is what keeps the surface closure
   % identity free of interior residuals whenever a clamp binds.

   [dz, delz, fn, ~, ~, f_ice, f_res_por] = coupledFixture(5);
   [Ls, ro_ice] = icemodel.physicalConstant('Ls', 'ro_ice');
   Tf = icemodel.physicalConstant('Tf');

   % A thin film at cell 2 makes the evaporation clamp bind.
   T = (Tf - 8) * ones(5, 1);
   T(2) = Tf - 2;
   f_res = icemodel.column.residual_water_fraction( ...
      f_ice, zeros(5, 1), f_res_por);
   f_liq = zeros(5, 1);
   f_liq(2) = f_res(2) + 1e-9;

   % The shortfall the applier records, reproduced from the primitives so
   % the orchestrated ledger has an independent expectation to meet.
   [ro_vap, De] = icemodel.column.accepted_vapor_quantities(T, f_liq);
   d_vap = icemodel.column.couple_vapor_transport( ...
      ro_vap, De, dz, delz, fn, 900);
   [f_ice_new, f_liq_new, d_sbl_err] = ...
      applyTransport(f_ice, f_liq, d_vap, 0.1, f_res_por);
   testCase.assertLessThan(d_sbl_err(2), 0);

   ledger = icemodel.column.initialize_budget_state();
   [~, ~, ~, ~, ~, ledger_new] = icemodel.column.couple_vapor_step( ...
      T, f_ice, f_liq, f_ice, f_liq, zeros(6, 1), 0.0, 0.0, dz, delz, ...
      fn, 900, 0.1, f_res_por, ledger, true);

   % The redistribution channels are the realized per-phase storage moves,
   % to roundoff.
   [solid_0, liquid_0] = icemodel.column.integrate_column_budget( ...
      T, f_ice, f_liq, dz);
   [solid_1, liquid_1] = icemodel.column.integrate_column_budget( ...
      T, f_ice_new, f_liq_new, dz);
   testCase.verifyEqual( ...
      ledger_new.mass_budget_vapor_redistribution_solid_mwe, ...
      solid_1 - solid_0, 'RelTol', 1e-12);
   testCase.verifyEqual( ...
      ledger_new.mass_budget_vapor_redistribution_liquid_mwe, ...
      liquid_1 - liquid_0, 'RelTol', 1e-12);

   % The shortfall lands in the redistribution's own unapplied channel, on
   % the same energy basis the surface channel uses.
   expected = sum(ro_ice * Ls * d_sbl_err(:) .* dz(:));
   testCase.verifyEqual( ...
      ledger_new.mass_budget_vapor_redistribution_unapplied_j_m2, ...
      expected, 'RelTol', 1e-12);
   testCase.verifyGreaterThan(abs(expected), 1);

   % The surface unapplied channel never sees the interior shortfall.
   testCase.verifyEqual( ...
      ledger_new.mass_budget_unapplied_vapor_j_m2, 0, 'AbsTol', 0);
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
   d_vap_faces_in = zeros(7, 1);
   [f_ice_a, f_liq_a, d_vap_faces_a, ~, ~, ledger_a] = ...
      icemodel.column.couple_vapor_step( ...
      T, f_ice, f_liq, f_ice, f_liq, d_vap_faces_in, 0.0, 0.0, dz, ...
      delz, fn, dt, f_ice_min, f_res_por, ledger, true);

   % The same primitives, called directly in the driver's order.
   ro_liq = icemodel.physicalConstant('ro_liq');
   [ro_vap, De] = icemodel.column.accepted_vapor_quantities(T, f_liq);
   [d_vap, ~, U_vap_faces] = icemodel.column.couple_vapor_transport( ...
      ro_vap, De, dz, delz, fn, dt);
   [solid_r, liquid_r] = icemodel.column.integrate_column_budget( ...
      T, f_ice, f_liq, dz);
   [f_ice_b, f_liq_b, d_sbl_err_b] = ...
      applyTransport(f_ice, f_liq, d_vap, f_ice_min, f_res_por);
   ledger_b = icemodel.column.accumulate_redistribution_budget( ...
      ledger, solid_r, liquid_r, T, f_ice_b, f_liq_b, dz, d_sbl_err_b);

   testCase.verifyEqual(f_ice_a, f_ice_b, 'AbsTol', 0);
   testCase.verifyEqual(f_liq_a, f_liq_b, 'AbsTol', 0);
   testCase.verifyEqual(d_vap_faces_a, abs(U_vap_faces) * dt / ro_liq, ...
      'AbsTol', 0);
   testCase.verifyEqual( ...
      ledger_a.mass_budget_vapor_redistribution_solid_mwe, ...
      ledger_b.mass_budget_vapor_redistribution_solid_mwe, 'AbsTol', 0);
   testCase.verifyEqual( ...
      ledger_a.mass_budget_vapor_redistribution_liquid_mwe, ...
      ledger_b.mass_budget_vapor_redistribution_liquid_mwe, 'AbsTol', 0);
   testCase.verifyEqual( ...
      ledger_a.mass_budget_vapor_redistribution_unapplied_j_m2, ...
      ledger_b.mass_budget_vapor_redistribution_unapplied_j_m2, ...
      'AbsTol', 0);

   % The column actually exchanged, or the comparison above is trivial.
   testCase.verifyTrue(any(f_ice_a ~= f_ice));
   testCase.verifyTrue(any(d_vap_faces_a > 0));
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
   d_vap_faces_in = zeros(7, 1);

   [f_ice_on, f_liq_on, faces_on, ~, ~, ledger_on] = ...
      icemodel.column.couple_vapor_step( ...
      T, f_ice, f_liq, f_ice, f_liq, d_vap_faces_in, 0.0, 0.0, dz, ...
      delz, fn, dt, f_ice_min, f_res_por, ledger, true);
   [f_ice_off, f_liq_off, faces_off, ~, ~, ledger_off] = ...
      icemodel.column.couple_vapor_step( ...
      T, f_ice, f_liq, f_ice, f_liq, d_vap_faces_in, 0.0, 0.0, dz, ...
      delz, fn, dt, f_ice_min, f_res_por, ledger, false);

   % The physics does not depend on whether the ledger is being built.
   testCase.verifyEqual(f_ice_off, f_ice_on, 'AbsTol', 0);
   testCase.verifyEqual(f_liq_off, f_liq_on, 'AbsTol', 0);
   testCase.verifyEqual(faces_off, faces_on, 'AbsTol', 0);

   % The ledger comes back as it went in, and the on-run moved it, so the
   % check is not passing because nothing happened.
   testCase.verifyEqual(ledger_off, ledger, 'AbsTol', 0);
   testCase.verifyNotEqual( ...
      ledger_on.mass_budget_vapor_redistribution_solid_mwe, ...
      ledger.mass_budget_vapor_redistribution_solid_mwe);
end

function test_the_orchestrator_accumulates_into_the_incoming_record(testCase)
   % d_vap_faces threads in and out so the gross face exchange accumulates
   % across substeps like every other d_* increment. A record arriving from
   % earlier substeps must survive, and every added increment is a
   % magnitude, so the accumulation never decreases.

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
      T, f_ice, f_liq, f_ice, f_liq, zeros(6, 1), 0.0, 0.0, dz, delz, ...
      fn, 900, 0.1, f_res_por, ledger, true);
   testCase.assertTrue(any(from_zero ~= 0));
   testCase.verifyTrue(all(from_zero >= 0));

   % The same call with a record already in it must return the sum.
   prior = [1e-7; 0; 2e-7; 0; 0; 3e-8];
   [~, ~, from_prior] = icemodel.column.couple_vapor_step( ...
      T, f_ice, f_liq, f_ice, f_liq, prior, 0.0, 0.0, dz, delz, ...
      fn, 900, 0.1, f_res_por, ledger, true);

   testCase.verifyEqual(from_prior, prior + from_zero, 'AbsTol', 0);
end

function test_transport_stays_out_of_the_surface_vapor_channels(testCase)
   % The surface vapor channels score surface exchange alone, because the
   % ablation comparator reads them as surface loss. The transport runs
   % after the vapor budget closes, so with no surface exchange at all the
   % surface vapor channels must stay exactly zero while the transport
   % moves real mass and its own channels record it. This exercises the
   % driver's accepted-substep order end to end.

   [dz, delz, fn, ~, ~, f_ice, f_res_por] = coupledFixture(6);
   Tf = icemodel.physicalConstant('Tf');
   T = (Tf - 8) + linspace(0, 6, 6)';
   f_liq = [0; 0; 0.05; 0.05; 0; 0];
   ledger = icemodel.column.initialize_budget_state();

   % The driver's accepted-substep tail, in its order: surface budgets
   % first, transport after. No surface exchange: d_pevp is zero.
   [solid_p, liquid_p] = icemodel.column.integrate_column_budget( ...
      T, f_ice, f_liq, dz);
   [T2, f_ice2, f_liq2, d_liq, ~, d_rof, d_sbl_err] = ...
      icemodel.column.budget_surface_mass_balance( ...
      T, f_ice, f_liq, f_liq, 0.0, 0.0, 0.0, 0.0, f_res_por, 0.1);
   ledger = icemodel.column.accumulate_vapor_budget(ledger, ...
      solid_p, liquid_p, T2, f_ice2, f_liq2, dz, 0.0, d_rof, d_sbl_err);
   [f_ice3, ~, ~, ~, ~, ledger] = icemodel.column.couple_vapor_step( ...
      T2, f_ice2, f_liq2, f_ice2, f_liq2, zeros(7, 1), 0.0, 0.0, dz, ...
      delz, fn, 900, 0.1, f_res_por, ledger, true);

   % The transport moved mass, and none of it reached the surface channels.
   testCase.assertTrue(any(f_ice3 ~= f_ice2));
   testCase.verifyEqual(ledger.mass_budget_vapor_solid_mwe, 0, 'AbsTol', 0);
   testCase.verifyEqual(ledger.mass_budget_vapor_liquid_mwe, 0, 'AbsTol', 0);
   testCase.verifyNotEqual( ...
      ledger.mass_budget_vapor_redistribution_solid_mwe, 0);

   % The melt and freeze reading stays clean: no phase change ran, so the
   % accumulated d_liq is zero in every cell, and the transport that moved
   % f_liq after the surface budgets closed can never enter it.
   testCase.verifyEqual(d_liq, zeros(6, 1), 'AbsTol', 0);
end

function test_the_realized_surface_exchange_excludes_rejected_demand(testCase)
   % Grain growth consumes the realized surface exchange, so the
   % realized-exchange output must carry only the mass that actually
   % crossed the surface. Demand a clamp rejects never crossed, and
   % growing grains from it would grow them from mass the column never
   % took.

   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   f_ice_min = 0.1;
   f_res_por = 0.02;

   % A dry top cell just above the ice floor, asked for far more than it
   % holds. The clamp binds and most of the demand is rejected.
   JJ = 3;
   T = (icemodel.physicalConstant('Tf') - 5) * ones(JJ, 1);
   f_ice = [f_ice_min + 1e-6; 0.85; 0.85];
   f_liq = zeros(JJ, 1);
   d_pevp = -1e-3;

   [~, f_ice_new, ~, ~, ~, ~, d_sbl_err, d_applied] = ...
      icemodel.column.budget_surface_mass_balance( ...
      T, f_ice, f_liq, f_liq, d_pevp, 0.0, 0.0, 0.0, ...
      f_res_por, f_ice_min);

   % The clamp bound: the cell stopped at the floor and recorded the rest.
   testCase.assertEqual(f_ice_new(1), f_ice_min, 'AbsTol', 1e-15);
   testCase.assertLessThan(d_sbl_err, 0);

   % The realized exchange is the tiny sliver above the floor, on the
   % liquid-equivalent basis, not the demand.
   expected = -(f_ice(1) - f_ice_min) * ro_ice / ro_liq;
   testCase.verifyEqual(d_applied, expected, 'RelTol', 1e-12);
   testCase.verifyLessThan(abs(d_applied), ...
      abs(icemodel.surface.potential_surface_vapor_exchange( ...
      d_pevp, f_ice(1), f_liq(1), f_res_por)));
end

function test_routing_reads_the_solve_state_not_the_current_state(testCase)
   % The phase a cell exchanges with comes from the solve state, because
   % that is the state the solve's face energy carried, while the amounts
   % clamp against the current fractions. A cell the surface exchange
   % dried between the two states must still route as wet: it then has no
   % mobile liquid to give, so its ice stays untouched and the whole
   % demand lands in the shortfall on the wet-branch basis. Routing from
   % the current state would sublimate ice instead.

   f_res_por = 0.02;
   f_ice = 0.85 * ones(3, 1);
   f_liq = zeros(3, 1);
   f_liq_solve = [0; 0.05; 0];
   d_vap = [0; -1e-5; 0];

   [wet, f_res] = icemodel.column.vapor_exchange_is_wet( ...
      f_ice, f_liq_solve, f_res_por);
   testCase.assertTrue(wet(2));

   [f_ice_solve_route, ~, err_solve_route] = ...
      icemodel.column.apply_vapor_transport( ...
      f_ice, f_liq, d_vap, 0.1, wet, f_res);

   % Wet routing with no current liquid: the ice is untouched and the
   % demand is recorded on the wet-branch basis.
   testCase.verifyEqual(f_ice_solve_route, f_ice, 'AbsTol', 0);
   testCase.verifyEqual(err_solve_route(2), ...
      icemodel.column.potential_sublimation(d_vap(2)), 'RelTol', 1e-14);

   % The distinction is observable: current-state routing classifies the
   % cell dry and sublimates its ice.
   [f_ice_now_route] = applyTransport(f_ice, f_liq, d_vap, 0.1, f_res_por);
   testCase.verifyLessThan(f_ice_now_route(2), f_ice(2));

   % The same distinction through the orchestrator, because
   % couple_vapor_step owns the predicate evaluation: with a dry current
   % state and a wet solve state, wet routing keeps every exchange in the
   % liquid phase, so the ice is untouched while the transport still
   % moves liquid. A regression that evaluated the predicate from the
   % current fractions would route the same transport to the ice phase,
   % as the contrast run shows.
   [dz, delz, fn, ~, ~, f_ice_col, f_res_por_c] = coupledFixture(6);
   T = gradientColumn(6);
   f_liq_cur = zeros(6, 1);
   f_liq_solve = 0.05 * ones(6, 1);
   ledger = icemodel.column.initialize_budget_state();

   [f_ice_wet, f_liq_wet] = icemodel.column.couple_vapor_step( ...
      T, f_ice_col, f_liq_cur, f_ice_col, f_liq_solve, zeros(7, 1), ...
      0.0, 0.0, dz, delz, fn, 900, 0.1, f_res_por_c, ledger, true);
   testCase.verifyEqual(f_ice_wet, f_ice_col, 'AbsTol', 0);
   testCase.verifyGreaterThan(max(abs(f_liq_wet - f_liq_cur)), 0);

   [f_ice_dry] = icemodel.column.couple_vapor_step( ...
      T, f_ice_col, f_liq_cur, f_ice_col, f_liq_cur, zeros(7, 1), ...
      0.0, 0.0, dz, delz, fn, 900, 0.1, f_res_por_c, ledger, true);
   testCase.verifyNotEqual(f_ice_dry, f_ice_col);
end

function test_the_context_gains_exactly_the_transport_increments(testCase)
   % The vapor storage-change context threads through the coupled step so
   % the remesh budget's endpoint gross sees the transport. The context
   % must gain exactly the per-phase increments the redistribution
   % channels record, no more and no less.

   [dz, delz, fn, ~, ~, f_ice, f_res_por] = coupledFixture(6);
   T = gradientColumn(6);
   f_liq = [0; 0; 0.05; 0.05; 0; 0];
   ledger = icemodel.column.initialize_budget_state();

   context_in_solid = 1.25e-4;
   context_in_liquid = -3.5e-5;
   [~, ~, ~, vs_out, vl_out, ledger_out] = ...
      icemodel.column.couple_vapor_step( ...
      T, f_ice, f_liq, f_ice, f_liq, zeros(7, 1), context_in_solid, ...
      context_in_liquid, dz, delz, fn, 900, 0.1, f_res_por, ledger, true);

   testCase.assertNotEqual( ...
      ledger_out.mass_budget_vapor_redistribution_solid_mwe, 0);
   testCase.verifyEqual(vs_out - context_in_solid, ...
      ledger_out.mass_budget_vapor_redistribution_solid_mwe, 'AbsTol', 0);
   testCase.verifyEqual(vl_out - context_in_liquid, ...
      ledger_out.mass_budget_vapor_redistribution_liquid_mwe, 'AbsTol', 0);
end

function test_the_realized_exchange_includes_condensation_overflow(testCase)
   % Overflow condensate crossed the surface before it ran off, so the
   % realized exchange must count it: with the demand above capacity,
   % everything crossed and nothing was rejected, so the realized
   % exchange equals the full demand rather than the retained capacity.

   f_res_por = 0.02;
   f_ice = 0.99;
   f_liq = 0.005;
   testCase.assertTrue(icemodel.column.vapor_exchange_is_wet( ...
      f_ice, f_liq, f_res_por));

   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   capacity = ro_ice / ro_liq * (1.0 - f_ice) - f_liq;
   testCase.assertGreaterThan(capacity, 0);
   d_pevp = 2 * capacity;

   [~, f_liq_new, d_rof, d_sbl_err, d_applied] = ...
      icemodel.surface.apply_surface_vapor_exchange( ...
      f_ice, f_liq, 0.0, d_pevp, 0.1, f_res_por);

   testCase.assertGreaterThan(d_rof, 0);
   testCase.verifyEqual(d_sbl_err, zeros(1, 1), 'AbsTol', 0);
   testCase.verifyEqual(f_liq_new, f_liq + capacity, 'RelTol', 1e-14);
   testCase.verifyEqual(d_applied, d_pevp, 'RelTol', 1e-14);
end

function test_the_gross_shortfall_survives_opposite_signs(testCase)
   % One substep can reject deposition in one cell and starve sublimation
   % in another. The signed shortfall channel nets the two; the gross
   % channel must keep both magnitudes, and both must accumulate across
   % substeps rather than reset, or opposite-sign shortfalls disappear
   % from the record.

   dz = 0.04 * ones(3, 1);
   f_res_por = 0.02;
   f_ice_min = 0.1;
   [ro_ice, ro_liq, Ls, Tf] = icemodel.physicalConstant( ...
      'ro_ice', 'ro_liq', 'Ls', 'Tf');
   T = (Tf - 5) * ones(3, 1);

   % Cell 1 is nearly full and rejects deposition. Cell 2 sits at the ice
   % floor and cannot supply sublimation. Cell 3 idles.
   f_ice = [0.999; f_ice_min + 1e-6; 0.85];
   f_liq = zeros(3, 1);
   f_air = 1.0 - f_ice(1);
   d_vap = [2 * f_air * ro_ice / ro_liq; -1e-3; 0];

   [~, ~, d_sbl_err] = applyTransport( ...
      f_ice, f_liq, d_vap, f_ice_min, f_res_por);
   testCase.assertGreaterThan(d_sbl_err(1), 0);
   testCase.assertLessThan(d_sbl_err(2), 0);

   % Passing the pre-transport state as both baseline and after-state
   % isolates the shortfall channels: the per-phase increments stay zero.
   [solid_r, liquid_r] = icemodel.column.integrate_column_budget( ...
      T, f_ice, f_liq, dz);
   ledger = icemodel.column.initialize_budget_state();
   ledger = icemodel.column.accumulate_redistribution_budget( ...
      ledger, solid_r, liquid_r, T, f_ice, f_liq, dz, d_sbl_err);

   weighted = ro_ice * Ls * d_sbl_err .* dz;
   testCase.verifyEqual( ...
      ledger.mass_budget_vapor_redistribution_unapplied_j_m2, ...
      sum(weighted), 'RelTol', 1e-14);
   testCase.verifyEqual( ...
      ledger.mass_budget_vapor_redistribution_unapplied_gross_j_m2, ...
      sum(abs(weighted)), 'RelTol', 1e-14);
   testCase.verifyGreaterThan( ...
      ledger.mass_budget_vapor_redistribution_unapplied_gross_j_m2, ...
      abs(ledger.mass_budget_vapor_redistribution_unapplied_j_m2));

   % A second substep accumulates rather than resets.
   ledger = icemodel.column.accumulate_redistribution_budget( ...
      ledger, solid_r, liquid_r, T, f_ice, f_liq, dz, d_sbl_err);
   testCase.verifyEqual( ...
      ledger.mass_budget_vapor_redistribution_unapplied_gross_j_m2, ...
      2 * sum(abs(weighted)), 'RelTol', 1e-14);
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
      applyTransport(f_ice, f_liq, d_vap, 0.1, f_res_por);

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
      applyTransport(f_ice, f_liq, d_vap, 0.1, f_res_por);

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
      applyTransport(f_ice, f_liq, d_vap, 0.1, f_res_por);

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
      applyTransport(f_ice, f_liq, d_vap, f_ice_min, f_res_por);

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

function [f_ice, f_liq, d_sbl_err] = applyTransport( ...
      f_ice, f_liq, d_vap, f_ice_min, f_res_por)
   %APPLYTRANSPORT Apply transport with the predicate at the applied state.
   %
   % The production applier takes the phase decision as an input, so its
   % caller controls the evaluation state. These fixtures apply no surface
   % exchange between the decision and the application, so the applied
   % state is the decision state, which is what the coupled step supplies
   % in production.

   [wet, f_res] = icemodel.column.vapor_exchange_is_wet( ...
      f_ice, f_liq, f_res_por);
   [f_ice, f_liq, d_sbl_err] = icemodel.column.apply_vapor_transport( ...
      f_ice, f_liq, d_vap, f_ice_min, wet, f_res);
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
