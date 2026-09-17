function tests = test_coupled_vapor_transport
   %TEST_COUPLED_VAPOR_TRANSPORT Verify vapor mass transport.
   %
   % Interior transport and surface exchange conserve different quantities.
   %
   % The accepted face flux moves vapor between cells by Fick's law, which
   % fixes the MASS. Because both transport boundary faces are closed,
   % transport redistributes mass and creates none.
   % icemodel.column.apply_vapor_transfer removes mass from the donor cell and
   % deposits the same phase in the receiver cell.
   %
   % Surface exchange remains separate from interior transport. The surface
   % energy balance fixes the ENERGY and lets the mass follow. Combining
   % surface exchange with interior transport loses mass.
   % A cell can spend part of an energy demand at one latent heat and the rest
   % at another. The applied mass then differs from the mass that arrived.
   %
   % See also: icemodel.column.vapor_transport_terms,
   %  icemodel.column.apply_vapor_transfer
   tests = functiontests(localfunctions);
end

function test_isothermal_column_redistributes_nothing(testCase)
   % No interior vapor gradient means no interior flux. The surface exchange
   % is not part of this step, so an isothermal column moves nothing at all.

   [dz, delz, fn, T_ice, f_liq, f_ice, f_res_por] = coupledFixture(8);
   [d_vap, dm_vap] = acceptedTransport( ...
      T_ice, f_ice, f_liq, dz, delz, fn, f_res_por, 900);

   testCase.verifyEqual(d_vap, zeros(8, 1), 'AbsTol', 1e-30);
   testCase.verifyEqual(dm_vap, zeros(8, 1), 'AbsTol', 1e-30);
end

function test_redistribution_conserves_column_mass(testCase)
   % Both boundaries are closed, so the transport moves mass between cells
   % and creates none. The column total must come back to zero.

   [dz, delz, fn, ~, f_liq, f_ice, f_res_por] = coupledFixture(8);
   T_ice = gradientColumn(8);
   [d_vap, dm_vap] = acceptedTransport( ...
      T_ice, f_ice, f_liq, dz, delz, fn, f_res_por, 900);

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
   T_ice = gradientColumn(8);
   d_vap = acceptedTransport( ...
      T_ice, f_ice, f_liq, dz, delz, fn, f_res_por, 900);
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
   T_ice = gradientColumn(8);
   f_liq = [0; 0; 0.05; 0.05; 0; 0; 0.05; 0];
   d_vap = acceptedTransport( ...
      T_ice, f_ice, f_liq, dz, delz, fn, f_res_por, 900);
   [f_ice_new, f_liq_new, d_sbl_err] = ...
      applyTransport(f_ice, f_liq, d_vap, 0.1, f_res_por);

   testCase.verifyEqual(d_sbl_err, zeros(8, 1), 'AbsTol', 0);
   testCase.verifyTrue(any(f_liq_new ~= f_liq));
   testCase.verifyTrue(any(f_ice_new ~= f_ice));
   verifyMassMatches(testCase, f_ice, f_liq, f_ice_new, f_liq_new, ...
      d_vap, dz);
end

function test_a_thin_liquid_film_records_what_it_cannot_give(testCase)
   % A wet cell cannot supply more liquid mass than it holds above its
   % residual floor. Apply the available mass and record the shortfall.
   %
   % Transport still supplies the neighbouring cells after the wet cell
   % reaches its residual floor. D_SBL_ERR returns the donor shortfall, which
   % explains the applied column storage change.

   [dz, delz, fn, ~, ~, f_ice, f_res_por] = coupledFixture(5);

   % A temperature peak at cell 2 makes it the vapor-density maximum, so it
   % loses mass to both neighbours. A monotonic profile would move vapor
   % through that cell instead of out of it, and the film would never drain.
   Tf = icemodel.physicalConstant('Tf');
   T_ice = (Tf - 8) * ones(5, 1);
   T_ice(2) = Tf - 2;

   % Put that cell barely above its residual floor, so the loss exhausts the
   % film at once.
   f_res = icemodel.column.residual_water_fraction( ...
      f_ice, zeros(5, 1), f_res_por);
   f_liq = zeros(5, 1);
   f_liq(2) = f_res(2) + 1e-9;
   testCase.assertTrue(icemodel.column.vapor_exchange_is_wet( ...
      f_ice(2), f_liq(2), f_res_por));

   d_vap = acceptedTransport( ...
      T_ice, f_ice, f_liq, dz, delz, fn, f_res_por, 900);
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
   % max_liquid_fraction_change is the condensation limit for every
   % caller. The pore volume converts to a liquid-water-equivalent
   % fraction as if it were ice, so a cell that later refreezes the water
   % it took cannot exceed a full control volume.

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
   % When a cell reaches a storage limit, the transport budget records only
   % the applied phase changes. APPLY_VAPOR_TRANSFER returns the unapplied
   % increments, and the surface-exchange budget fields remain zero.

   [dz, delz, fn, ~, ~, f_ice, f_res_por] = coupledFixture(5);
   [Ls, ro_ice] = icemodel.physicalConstant('Ls', 'ro_ice');
   Tf = icemodel.physicalConstant('Tf');

   % A thin film at cell 2 makes the evaporation clamp bind.
   T_ice = (Tf - 8) * ones(5, 1);
   T_ice(2) = Tf - 2;
   f_res = icemodel.column.residual_water_fraction( ...
      f_ice, zeros(5, 1), f_res_por);
   f_liq = zeros(5, 1);
   f_liq(2) = f_res(2) + 1e-9;

   % The shortfall the applier records, reproduced from the primitives so
   % the orchestrated budget has an independent expectation to meet.
   [~, ~, U_vap_faces, L_vap_faces] = acceptedTransport( ...
      T_ice, f_ice, f_liq, dz, delz, fn, f_res_por, 900);
   [f_ice_new, f_liq_new, d_sbl_err] = applyDonorTransport( ...
      f_ice, f_liq, U_vap_faces, L_vap_faces, 900, dz, 0.1, f_res_por);
   testCase.assertLessThan(d_sbl_err(2), 0);

   budget = icemodel.column.initialize_budget_state(T_ice, f_ice, f_liq, dz);
   [~, ~, ~, ~, budget_new] = icemodel.column.couple_vapor_step( ...
      f_ice, f_liq, U_vap_faces, L_vap_faces, zeros(5, 1), zeros(5, 1), ...
      dz, 900, 0.1, f_res_por, budget);

   % The transport channels are the realized per-phase storage moves, to
   % roundoff.
   [solid_1, liquid_1] = icemodel.column.integrate_column_budget( ...
      T_ice, f_ice_new, f_liq_new, dz);
   testCase.verifyEqual( ...
      budget_new.mass_budget_vapor_transport_solid_mwe, ...
      solid_1 - budget.mass_budget_solid_start_mwe, 'RelTol', 1e-12);
   testCase.verifyEqual( ...
      budget_new.mass_budget_vapor_transport_liquid_mwe, ...
      liquid_1 - budget.mass_budget_liquid_start_mwe, 'RelTol', 1e-12);

   % The independent calculation expresses the shortfall as ro_ice * Ls *
   % d_sbl_err. Require a nontrivial value so this test exercises the limit.
   expected_shortfall_energy = sum(ro_ice * Ls * d_sbl_err(:) .* dz(:));
   testCase.verifyGreaterThan(abs(expected_shortfall_energy), 1);

   % COUPLE_VAPOR_STEP does not add an interior shortfall to the surface
   % exchange budget fields.
   testCase.verifyEqual( ...
      budget_new.mass_budget_vapor_solid_mwe, 0, 'AbsTol', 0);
   testCase.verifyEqual( ...
      budget_new.mass_budget_vapor_liquid_mwe, 0, 'AbsTol', 0);
end

function test_the_orchestrator_matches_its_parts(testCase)
   % Compare COUPLE_VAPOR_STEP with its transfer and budget functions called
   % in the same order.

   [dz, delz, fn, ~, ~, f_ice, f_res_por] = coupledFixture(6);
   T_ice = gradientColumn(6);
   f_liq = [0; 0; 0.05; 0.05; 0; 0];
   dt = 900;
   f_ice_min = 0.1;
   budget = icemodel.column.initialize_budget_state(T_ice, f_ice, f_liq, dz);

   % Run the combined function.
   [~, ~, U_vap_faces, L_vap_faces] = acceptedTransport( ...
      T_ice, f_ice, f_liq, dz, delz, fn, f_res_por, dt);
   [f_ice_a, f_liq_a, ~, ~, budget_a] = icemodel.column.couple_vapor_step( ...
      f_ice, f_liq, U_vap_faces, L_vap_faces, zeros(6, 1), zeros(6, 1), ...
      dz, dt, f_ice_min, f_res_por, budget);

   % The same primitives, called directly in the driver's order.
   [f_ice_b, f_liq_b, ~, d_liq_applied, d_ice_applied] = applyDonorTransport( ...
      f_ice, f_liq, U_vap_faces, L_vap_faces, dt, dz, f_ice_min, f_res_por);
   budget_b = icemodel.column.accumulate_vapor_transport( ...
      budget, d_liq_applied, d_ice_applied, dz);

   testCase.verifyEqual(f_ice_a, f_ice_b, 'AbsTol', 0);
   testCase.verifyEqual(f_liq_a, f_liq_b, 'AbsTol', 0);
   testCase.verifyEqual( ...
      budget_a.mass_budget_vapor_transport_solid_mwe, ...
      budget_b.mass_budget_vapor_transport_solid_mwe, 'AbsTol', 0);
   testCase.verifyEqual( ...
      budget_a.mass_budget_vapor_transport_liquid_mwe, ...
      budget_b.mass_budget_vapor_transport_liquid_mwe, 'AbsTol', 0);

   % The column actually exchanged, or the comparison above is trivial.
   testCase.verifyTrue(any(f_ice_a ~= f_ice));
end

function test_the_budget_accumulates_across_repeated_transport_calls(testCase)
   % ACCUMULATE_VAPOR_TRANSPORT adds each accepted substep to the forcing-step
   % budget. Two calls must retain the sum of both transport increments.

   [dz, delz, fn, ~, ~, f_ice, f_res_por] = coupledFixture(6);
   T_ice = gradientColumn(6);
   f_liq = [0; 0; 0.05; 0.05; 0; 0];
   budget = icemodel.column.initialize_budget_state(T_ice, f_ice, f_liq, dz);

   [~, ~, U_vap_1, L_vap_1] = acceptedTransport( ...
      T_ice, f_ice, f_liq, dz, delz, fn, f_res_por, 900);
   [f_ice_1, f_liq_1, d_vap_liq_1, d_vap_ice_1, budget_1] = ...
      icemodel.column.couple_vapor_step( ...
      f_ice, f_liq, U_vap_1, L_vap_1, zeros(6, 1), zeros(6, 1), dz, ...
      900, 0.1, f_res_por, budget);
   first_solid = budget_1.mass_budget_vapor_transport_solid_mwe;
   first_liquid = budget_1.mass_budget_vapor_transport_liquid_mwe;
   testCase.assertNotEqual(first_solid, 0);

   % A second accepted substep transports again from the state the first
   % substep left behind, chaining the accumulator outputs forward the way
   % the driver chains substeps.
   [~, ~, U_vap_2, L_vap_2] = acceptedTransport( ...
      T_ice, f_ice_1, f_liq_1, dz, delz, fn, f_res_por, 900);
   [~, ~, ~, ~, budget_2] = icemodel.column.couple_vapor_step( ...
      f_ice_1, f_liq_1, U_vap_2, L_vap_2, d_vap_liq_1, d_vap_ice_1, dz, ...
      900, 0.1, f_res_por, budget_1);

   % A fresh single-call budget isolates the second transport's own
   % contribution, so the chained total can be checked against the sum of
   % two independent calls.
   budget_solo_2 = icemodel.column.initialize_budget_state( ...
      T_ice, f_ice_1, f_liq_1, dz);
   [~, ~, ~, ~, budget_solo_2] = icemodel.column.couple_vapor_step( ...
      f_ice_1, f_liq_1, U_vap_2, L_vap_2, zeros(6, 1), zeros(6, 1), dz, ...
      900, 0.1, f_res_por, budget_solo_2);
   second_solid = budget_solo_2.mass_budget_vapor_transport_solid_mwe;
   second_liquid = budget_solo_2.mass_budget_vapor_transport_liquid_mwe;
   testCase.assertNotEqual(second_solid, 0);

   testCase.verifyEqual( ...
      budget_2.mass_budget_vapor_transport_solid_mwe, ...
      first_solid + second_solid, 'AbsTol', 0);
   testCase.verifyEqual( ...
      budget_2.mass_budget_vapor_transport_liquid_mwe, ...
      first_liquid + second_liquid, 'AbsTol', 0);
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
   T_ice = (Tf - 8) + linspace(0, 6, 6)';
   f_liq = [0; 0; 0.05; 0.05; 0; 0];
   budget = icemodel.column.initialize_budget_state(T_ice, f_ice, f_liq, dz);

   % Run the accepted-substep order: surface budget first, then transport.
   % D_PEVP is zero, so the surface call returns zero vapor increments.
   [T_ice2, f_ice2, f_liq2, d_liq, ~, ~, d_vap_liq, d_vap_ice, ~, budget] = ...
      icemodel.column.budget_surface_mass_balance( ...
      T_ice, f_ice, f_liq, f_liq, 0.0, zeros(6, 1), zeros(6, 1), 0.0, ...
      zeros(6, 1), zeros(6, 1), f_res_por, 0.1, budget, dz);
   [~, ~, U_vap_faces, L_vap_faces] = acceptedTransport( ...
      T_ice2, f_ice2, f_liq2, dz, delz, fn, f_res_por, 900);
   [f_ice3, ~, ~, ~, budget] = icemodel.column.couple_vapor_step( ...
      f_ice2, f_liq2, U_vap_faces, L_vap_faces, d_vap_liq, d_vap_ice, ...
      dz, 900, 0.1, f_res_por, budget);

   % Transport moved mass without changing the surface budget fields.
   testCase.assertTrue(any(f_ice3 ~= f_ice2));
   testCase.verifyEqual(budget.mass_budget_vapor_solid_mwe, 0, 'AbsTol', 0);
   testCase.verifyEqual(budget.mass_budget_vapor_liquid_mwe, 0, 'AbsTol', 0);
   testCase.verifyNotEqual( ...
      budget.mass_budget_vapor_transport_solid_mwe, 0);

   % The melt and freeze reading stays clean: no phase change ran, so the
   % accumulated d_liq is zero in every cell, and the transport that moved
   % f_liq after the surface budgets closed can never enter it.
   testCase.verifyEqual(d_liq, zeros(6, 1), 'AbsTol', 0);
end

function test_the_realized_surface_exchange_excludes_rejected_demand(testCase)
   % Grain growth uses the realized surface exchange. The output must exclude
   % demand that the column cannot supply after every cell reaches its floor.

   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   f_ice_min = 0.1;
   f_res_por = 0.02;

   % Every cell just above the ice floor, asked for far more than the whole
   % column holds. The cascade draws every cell down to the floor in turn
   % and still cannot satisfy the demand.
   JJ = 3;
   T_ice = (icemodel.physicalConstant('Tf') - 5) * ones(JJ, 1);
   f_ice = (f_ice_min + 1e-6) * ones(JJ, 1);
   f_liq = zeros(JJ, 1);
   d_pevp = -1e-3;
   dz = 0.04 * ones(JJ, 1);
   budget = icemodel.column.initialize_budget_state(T_ice, f_ice, f_liq, dz);

   [~, f_ice_new, ~, ~, ~, ~, ~, ~, d_applied, ~] = ...
      icemodel.column.budget_surface_mass_balance( ...
      T_ice, f_ice, f_liq, f_liq, d_pevp, zeros(JJ, 1), zeros(JJ, 1), 0.0, ...
      zeros(JJ, 1), zeros(JJ, 1), f_res_por, f_ice_min, budget, dz);

   % Every cell stops at the floor. The realized exchange shows how much mass
   % the column supplied; no budget field stores the unmet remainder.
   testCase.assertEqual(f_ice_new, f_ice_min * ones(JJ, 1), 'AbsTol', 1e-15);

   % The realized exchange is what the whole column could give above the
   % floor, on the liquid-equivalent basis, not the demand.
   expected = -sum(f_ice - f_ice_min) * ro_ice / ro_liq;
   testCase.verifyEqual(d_applied, expected, 'RelTol', 1e-12);
   testCase.verifyLessThan(abs(d_applied), abs(d_pevp));
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

   [f_ice_solve_route, ~, liq_unapplied, ice_unapplied] = ...
      icemodel.column.apply_vapor_transfer( ...
      f_ice, f_liq, d_vap .* wet, d_vap .* ~wet, 0.1, f_res);
   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   err_solve_route = icemodel.column.potential_sublimation(liq_unapplied) ...
      + ice_unapplied * ro_liq / ro_ice;

   % Wet routing with no current liquid: the ice is untouched and the
   % demand is recorded on the wet-branch basis.
   testCase.verifyEqual(f_ice_solve_route, f_ice, 'AbsTol', 0);
   testCase.verifyEqual(err_solve_route(2), ...
      icemodel.column.potential_sublimation(d_vap(2)), 'RelTol', 1e-14);

   % The distinction is observable: current-state routing classifies the
   % cell dry and sublimates its ice.
   [f_ice_now_route] = applyTransport(f_ice, f_liq, d_vap, 0.1, f_res_por);
   testCase.verifyLessThan(f_ice_now_route(2), f_ice(2));

   % COUPLE_VAPOR_STEP selects the donor phase from L_VAP_FACES. A wet solve
   % state therefore keeps this transfer liquid even when the current cell is
   % dry. Recomputing the phase from the current fractions would remove ice.
   [dz, delz, fn, ~, ~, f_ice_col, f_res_por_c] = coupledFixture(6);
   T_ice = gradientColumn(6);
   f_liq_cur = zeros(6, 1);
   f_liq_solve = 0.05 * ones(6, 1);
   budget = icemodel.column.initialize_budget_state( ...
      T_ice, f_ice_col, f_liq_cur, dz);
   [~, ~, U_vap_wet, L_vap_wet] = acceptedTransport( ...
      T_ice, f_ice_col, f_liq_solve, dz, delz, fn, f_res_por_c, 900);
   [~, ~, U_vap_dry, L_vap_dry] = acceptedTransport( ...
      T_ice, f_ice_col, f_liq_cur, dz, delz, fn, f_res_por_c, 900);

   [f_ice_wet, f_liq_wet] = icemodel.column.couple_vapor_step( ...
      f_ice_col, f_liq_cur, U_vap_wet, L_vap_wet, zeros(6, 1), ...
      zeros(6, 1), dz, 900, 0.1, f_res_por_c, budget);
   testCase.verifyEqual(f_ice_wet, f_ice_col, 'AbsTol', 0);
   testCase.verifyGreaterThan(max(abs(f_liq_wet - f_liq_cur)), 0);

   [f_ice_dry] = icemodel.column.couple_vapor_step( ...
      f_ice_col, f_liq_cur, U_vap_dry, L_vap_dry, zeros(6, 1), ...
      zeros(6, 1), dz, 900, 0.1, f_res_por_c, budget);
   testCase.verifyNotEqual(f_ice_dry, f_ice_col);
end

function test_cross_phase_faces_preserve_the_donor_phase(testCase)
   % A dry donor sends vapor into a wet receiver. Both the removal and the
   % deposition must use the dry donor phase because the face-energy
   % calculation uses `Ls`: couple_vapor_step trusts L_vap over either
   % node's own current wetness.

   [~, ~, ~, ~, ~, f_ice, f_res_por] = coupledFixture(2);
   dz = [0.03; 0.05];
   T_ice = gradientColumn(2);
   f_liq = [0; 0.05];
   U_vap_faces = [0; 1e-5; 0];
   Ls = icemodel.physicalConstant('Ls');
   L_vap_faces = Ls * ones(3, 1);
   dt = 100;
   budget = icemodel.column.initialize_budget_state(T_ice, f_ice, f_liq, dz);

   [f_ice_new, f_liq_new] = icemodel.column.couple_vapor_step( ...
      f_ice, f_liq, U_vap_faces, L_vap_faces, zeros(2, 1), zeros(2, 1), ...
      dz, dt, 0.1, f_res_por, budget);

   ro_ice = icemodel.physicalConstant('ro_ice');
   expected_ice_change = 1e-5 * dt / dz(1) / ro_ice;
   expected_receiver_change = 1e-5 * dt / dz(2) / ro_ice;
   testCase.verifyEqual(f_ice(1) - f_ice_new(1), expected_ice_change, ...
      'RelTol', 1e-10);
   testCase.verifyEqual(f_ice_new(2) - f_ice(2), expected_receiver_change, ...
      'RelTol', 1e-10);
   testCase.verifyEqual(f_liq_new, f_liq, 'AbsTol', 0);
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

   [~, f_liq_new, d_rof, d_vap_liq, d_vap_ice, d_applied] = ...
      icemodel.surface.apply_surface_vapor_exchange( ...
      f_ice, f_liq, 0.0, d_pevp, 1, 0.1, f_res_por);

   testCase.assertGreaterThan(d_rof, 0);
   % The applied liquid increment is exactly the water capacity; the rest of
   % the demand ran off as overflow, and no ice was touched.
   testCase.verifyEqual(d_vap_liq, capacity, 'RelTol', 1e-14);
   testCase.verifyEqual(d_vap_ice, zeros(1, 1), 'AbsTol', 0);
   testCase.verifyEqual(f_liq_new, f_liq + capacity, 'RelTol', 1e-14);
   testCase.verifyEqual(d_applied, d_pevp, 'RelTol', 1e-14);
end

function test_the_shortfall_survives_opposite_signs(testCase)
   % One substep can reject deposition in one cell and starve sublimation in
   % another. APPLY_VAPOR_TRANSFER returns each cell's unapplied increment.
   % D_SBL_ERR keeps the opposite signs separate; their scalar sum cancels.

   f_res_por = 0.02;
   f_ice_min = 0.1;
   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');

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

   % The per-cell record keeps both signs. A single summed value would let
   % a positive and a negative shortfall of similar size cancel out.
   testCase.verifyGreaterThan( ...
      abs(d_sbl_err(1)) + abs(d_sbl_err(2)), abs(sum(d_sbl_err)));
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
   T_ice = gradientColumn(3);
   f_res_por = 0.07;

   % Dense ice: the residual floor sits below the fixed threshold.
   verifyBandConserves(testCase, 0.90 * ones(3, 1), [0; 0.0135; 0], ...
      f_res_por, T_ice, dz, delz, fn, true);

   % Low-density snow: the floor sits above it, so the disagreement flips.
   verifyBandConserves(testCase, 0.30 * ones(3, 1), [0; 0.03; 0], ...
      f_res_por, T_ice, dz, delz, fn, false);
end

function verifyBandConserves(testCase, f_ice, f_liq, f_res_por, T_ice, dz, ...
      delz, fn, expect_wet)
   %VERIFYBANDCONSERVES Check one column's applied mass against its transport.

   testCase.verifyEqual(icemodel.column.vapor_exchange_is_wet( ...
      f_ice(2), f_liq(2), f_res_por), expect_wet);

   d_vap = acceptedTransport( ...
      T_ice, f_ice, f_liq, dz, delz, fn, f_res_por, 900);
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
   % Call the applier directly. Testing this branch through the accepted face
   % flux needs a temperature field tuned to deliver a specific mass. That
   % setup tests the fixture rather than the branch.

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

function test_dry_deposition_reserves_residual_liquid_expansion(testCase)
   % Residual liquid prevents a dry cell from using all air space for vapor
   % deposition. Reserve capacity for retained liquid to expand during
   % refreezing so the next phase projection does not delete deposited mass.

   [ro_ice, ro_liq, Tf] = ...
      icemodel.physicalConstant('ro_ice', 'ro_liq', 'Tf');
   f_ice = 0.6;
   f_liq = 0.05;
   f_res_por = 0.2;
   [wet, f_res] = icemodel.column.vapor_exchange_is_wet( ...
      f_ice, f_liq, f_res_por);
   testCase.assertFalse(wet);

   capacity_lwe = icemodel.column.max_liquid_fraction_change(f_ice, f_liq);
   raw_air_lwe = (1.0 - f_ice - f_liq) * ro_ice / ro_liq;
   testCase.assertGreaterThan(raw_air_lwe, capacity_lwe);
   demand_lwe = 2 * raw_air_lwe;

   [f_ice_new, f_liq_new, ~, d_ice_unapplied] = ...
      icemodel.column.apply_vapor_transfer( ...
      f_ice, f_liq, 0.0, demand_lwe, 0.1, f_res);
   testCase.verifyEqual(f_liq_new, f_liq, 'AbsTol', 0);
   testCase.verifyEqual(d_ice_unapplied, demand_lwe - capacity_lwe, ...
      'RelTol', 1e-14);

   f_wat = icemodel.column.water_fraction(f_ice_new, f_liq_new);
   testCase.verifyEqual(f_wat, ro_ice / ro_liq, 'AbsTol', 2 * eps);
   testCase.verifyTrue(icemodel.column.assert_max_water(f_ice_new, f_liq_new));
   [~, f_ice_projected, f_liq_projected] = ...
      icemodel.column.liquid_fraction_function( ...
      Tf - 5, f_ice_new, f_liq_new);
   testCase.verifyEqual( ...
      icemodel.column.water_fraction(f_ice_projected, f_liq_projected), ...
      f_wat, 'AbsTol', 2 * eps);
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
   %VERIFYMASSANDRECORDMATCH Applied plus rejected must equal requested.
   %
   % D_SBL_ERR uses ice-fraction units. A dry cell records the unapplied ice
   % fraction. A wet cell scales its
   % unapplied liquid fraction by Lv/Ls before converting it to ice fraction.
   % Recover the mass with the conversion used for that phase.
   %
   % Use VAPOR_EXCHANGE_IS_WET so this inverse calculation uses the same phase
   % criterion as APPLY_VAPOR_TRANSFER.

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
   % These fixtures apply no surface exchange between the phase decision and
   % the mutation, so the applied state is also the decision state.

   [wet, f_res] = icemodel.column.vapor_exchange_is_wet( ...
      f_ice, f_liq, f_res_por);
   [f_ice, f_liq, liq_unapplied, ice_unapplied] = ...
      icemodel.column.apply_vapor_transfer(f_ice, f_liq, ...
      d_vap .* wet, d_vap .* ~wet, f_ice_min, f_res);

   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   d_sbl_err = icemodel.column.potential_sublimation(liq_unapplied) ...
      + ice_unapplied * ro_liq / ro_ice;
end

function [d_vap_nodes, dm_vap_nodes, U_vap_faces, L_vap_faces] = ...
      acceptedTransport(T_ice, f_ice, f_liq, dz, delz, fn, f_res_por, dt)
   %ACCEPTEDTRANSPORT Compute the accepted face flux used by transport fixtures.

   ro_liq = icemodel.physicalConstant('ro_liq');
   [ro_vap, dro_vapdT] = ...
      icemodel.vapor.saturation_vapor_density(T_ice, f_liq);
   [~, De] = icemodel.vapor.vapor_thermal_conductivity( ...
      T_ice, f_liq, dro_vapdT);
   k_eff = icemodel.column.bulk_thermal_conductivity( ...
      T_ice, f_ice, f_liq, zeros(size(T_ice)));
   [~, ~, ~, U_vap_faces, L_vap_faces] = ...
      icemodel.column.vapor_transport_terms( ...
      T_ice, f_ice, f_liq, k_eff, ro_vap, dro_vapdT, De, delz, fn, ...
      f_res_por);

   dm_vap_nodes = ...
      (U_vap_faces(1:end-1) - U_vap_faces(2:end)) ./ dz;
   d_vap_nodes = dm_vap_nodes * dt / ro_liq;
end

function T_ice = gradientColumn(JJ)
   %GRADIENTCOLUMN Return a temperature profile that drives interior vapor.

   Tf = icemodel.physicalConstant('Tf');
   T_ice = (Tf - 8) + linspace(0, 6, JJ)';
end

function [f_ice, f_liq, d_sbl_err, d_liq_applied, d_ice_applied] = ...
      applyDonorTransport( ...
      f_ice, f_liq, U_vap_faces, L_vap_faces, dt, dz, f_ice_min, f_res_por)
   %APPLYDONORTRANSPORT Reproduce couple_vapor_step's donor routing directly.
   %
   % Independent of icemodel.column.couple_vapor_step: it applies the same
   % L_VAP-selected donor phase through icemodel.column.apply_vapor_transfer.
   % L_VAP_FACES supplies the donor phase. TEST_VAPOR_FACE_DISCRETIZATION
   % verifies that solver output independently.
   %
   % D_LIQ_APPLIED and D_ICE_APPLIED are the demand minus the unapplied
   % remainder, the same accumulator inputs
   % icemodel.column.accumulate_vapor_transport expects.

   [Lv, ro_ice, ro_liq] = ...
      icemodel.physicalConstant('Lv', 'ro_ice', 'ro_liq');
   JJ = numel(f_ice);
   [~, f_res] = icemodel.column.vapor_exchange_is_wet( ...
      f_ice, f_liq, f_res_por);
   donor_wet = L_vap_faces(2:JJ) == Lv;
   d_vap_face = U_vap_faces(2:JJ) * dt / ro_liq;
   d_liq = zeros(JJ, 1);
   d_ice = zeros(JJ, 1);
   d_liq(1:JJ-1) = d_liq(1:JJ-1) ...
      - d_vap_face ./ dz(1:JJ-1) .* donor_wet;
   d_liq(2:JJ) = d_liq(2:JJ) ...
      + d_vap_face ./ dz(2:JJ) .* donor_wet;
   d_ice(1:JJ-1) = d_ice(1:JJ-1) ...
      - d_vap_face ./ dz(1:JJ-1) .* ~donor_wet;
   d_ice(2:JJ) = d_ice(2:JJ) ...
      + d_vap_face ./ dz(2:JJ) .* ~donor_wet;
   [f_ice, f_liq, d_liq_unapplied, d_ice_unapplied] = ...
      icemodel.column.apply_vapor_transfer( ...
      f_ice, f_liq, d_liq, d_ice, f_ice_min, f_res);
   d_sbl_err = icemodel.column.potential_sublimation(d_liq_unapplied) ...
      + ro_liq / ro_ice * d_ice_unapplied;
   d_liq_applied = d_liq - d_liq_unapplied;

   % The transport budget stores ice increments on the ice-fraction basis.
   d_ice_applied = (d_ice - d_ice_unapplied) * ro_liq / ro_ice;
end

function [dz, delz, fn, T_ice, f_liq, f_ice, f_res_por] = coupledFixture(JJ)
   %COUPLEDFIXTURE Return one dry, isothermal control-volume column.

   Tf = icemodel.physicalConstant('Tf');
   [dz, delz, ~, ~, fn] = ...
      icemodel.column.control_volume_mesh(JJ * 0.04, 0.04);
   dz = dz(1:JJ);
   delz = delz(1:JJ + 1);
   fn = fn(1:JJ + 1);
   T_ice = (Tf - 5) * ones(JJ, 1);
   f_liq = zeros(JJ, 1);
   f_ice = 0.85 * ones(JJ, 1);
   f_res_por = 0.02;
end
