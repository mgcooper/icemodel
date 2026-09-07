function tests = test_mass_budget_bookkeeping
   %TEST_MASS_BUDGET_BOOKKEEPING Verify the physical and remesh ledgers close.
   tests = functiontests(localfunctions);
end

function test_budget_state_uses_documented_references(testCase)
   % Storage must use metres water equivalent and the production dry-reference
   % bulk enthalpy integrated over the supplied mesh.

   [Tf, ro_ice, ro_liq] = ...
      icemodel.physicalConstant('Tf', 'ro_ice', 'ro_liq');
   T = [Tf - 2; Tf - 1];
   f_ice = [0.8; 0.7];
   f_liq = [0.01; 0.02];
   dz = [0.04; 0.06];

   [solid_mwe, liquid_mwe, enthalpy_j_m2] = ...
      icemodel.column.integrate_column_budget(T, f_ice, f_liq, dz);

   % Compare against the explicit accepted storage definitions.
   testCase.verifyEqual(solid_mwe, ...
      ro_ice / ro_liq * sum(f_ice .* dz), 'AbsTol', 1e-15);
   testCase.verifyEqual(liquid_mwe, sum(f_liq .* dz), 'AbsTol', 1e-15);
   f_wat = icemodel.column.water_fraction(f_ice, f_liq);
   H = icemodel.column.bulk_enthalpy(T, f_ice, f_liq, f_wat);
   testCase.verifyEqual(enthalpy_j_m2, sum(H .* dz), 'AbsTol', 1e-9);
end

function test_budgetoutputs_are_one_partitioned_contract(testCase)
   % Every diagnostic channel must belong to exactly one retime class.

   first_fields = icemodel.namelists.budgetoutputs('first');
   last_fields = icemodel.namelists.budgetoutputs('last');
   sum_fields = icemodel.namelists.budgetoutputs('sum');
   all_fields = icemodel.namelists.budgetoutputs();

   % The all-list order is the stable diagnostic-profile append order.
   % budgetoutputs lists 21 channels: signed nets and endpoints only,
   % with no gross channels, no enthalpy channels, and a solid-only
   % remesh decomposition.
   testCase.verifyEqual(all_fields, [first_fields, last_fields, sum_fields]);
   testCase.verifyEqual(numel(unique(all_fields)), numel(all_fields));
   testCase.verifyEqual(numel(all_fields), 21);
   testCase.verifyError( ...
      @() icemodel.namelists.budgetoutputs('median'), ...
      'icemodel:namelists:budgetoutputs:kind');
end

function test_budget_ledger_is_fixed_codegen_schema(testCase)
   % The kernel ledger must use one literal scalar-double layout whose order
   % remains synchronized with icemodel.namelists.budgetoutputs.

   Tf = icemodel.physicalConstant('Tf');
   T = Tf - 2;
   f_ice = 0.8;
   f_liq = 0.01;
   dz = 0.04;
   [solid_start, liquid_start] = ...
      icemodel.column.integrate_column_budget(T, f_ice, f_liq, dz);

   ledger = icemodel.column.initialize_budget_state(T, f_ice, f_liq, dz);
   budget_fields = icemodel.namelists.budgetoutputs('all');
   ledger_fields = transpose(fieldnames(ledger));

   % The struct holds exactly the 21 mass_budget_ channels, in
   % budgetoutputs order, with no other field. Every producer receives its
   % increments directly, so the struct carries no scratch baseline.
   testCase.verifyEqual(ledger_fields, budget_fields);
   testCase.verifyEqual(numel(ledger_fields), 21);

   channel_values = struct2cell(ledger);
   testCase.verifyTrue(all(cellfun( ...
      @(value) isa(value, 'double') && isequal(size(value), [1, 1]), ...
      channel_values)));

   % Every channel resets to zero except the entry-state storage endpoints,
   % which initialize_budget_state records from the supplied state instead of
   % zeroing.
   expected_values = zeros(1, numel(budget_fields));
   expected_values(1) = solid_start;
   expected_values(2) = liquid_start;
   testCase.verifyEqual([channel_values{:}], expected_values, 'AbsTol', 1e-15);

   % Keep runtime field-name construction out of the #codegen kernel. The
   % ordered field-name registry remains available only to MATLAB consumers.
   % Coder rejects dynamic field names and cell2struct, so a source check is
   % the only way to catch their reintroduction from a MATLAB-only suite.
   core_source = fileread(which('icemodel'));
   testCase.verifyFalse(contains(core_source, 'cell2struct'));
   testCase.verifyFalse(contains(core_source, 'surface_state.('));
   testCase.verifyFalse(contains(core_source, 'budgetoutputs('));
end

function test_diagnostic_model_payload_matches_canonical_budget_registry(testCase)
   % The runtime diagnostic payload must contain exactly the
   % budgetoutputs channel list.

   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=4, dt_seconds=3600);
   cleanup = onCleanup(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'icemodel', 2016, output_profile='diagnostic', solver=1);

   % Exercise initialization, timestep accumulation, payload assembly, and
   % raw output storage instead of checking configuration declarations alone.
   [ice1, ~] = icemodel.test.helpers.runSmbModel(opts);
   budget_fields = icemodel.namelists.budgetoutputs('all');
   payload_fields = transpose(fieldnames(ice1));
   payload_fields = payload_fields(startsWith(payload_fields, 'mass_budget_'));

   % Exact ordered equality detects both a missing registry field and any
   % independently added payload field that could drift from the registry.
   testCase.verifyEqual(payload_fields, budget_fields);
   for n = 1:numel(budget_fields)
      testCase.verifyClass(ice1.(budget_fields{n}), 'double');
      testCase.verifySize(ice1.(budget_fields{n}), [4, 1]);
      testCase.verifyTrue(all(isfinite(ice1.(budget_fields{n}))));
   end
   clear cleanup
end

function test_use_ro_glc_changes_initialization_not_diagnostic_basis(testCase)
   % Equal-density initialization must remain distinct while every diagnostic
   % ledger continues to use the solver's physical intrinsic-density basis.

   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=1, dt_seconds=3600);
   cleanup = onCleanup(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));
   base_opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'icemodel', 2016, output_profile='diagnostic', solver=1);
   [ro_ice, ro_liq, Ls, Lv] = ...
      icemodel.physicalConstant('ro_ice', 'ro_liq', 'Ls', 'Lv');
   ro_glc = (ro_ice + ro_liq) / 2;
   [dz, ~, ~, ~, ~] = icemodel.column.control_volume_mesh( ...
      base_opts.z0_thermal, base_opts.dz_thermal);
   solid_start = zeros(1, 2);

   % Exercise both initialization choices through the production model while
   % deriving diagnostic storage from physical solver densities in both cases.
   for k = 1:2
      use_ro_glc = logical(k - 1);
      opts = icemodel.resetopts(base_opts, 'use_ro_glc', use_ro_glc);
      [ice1, ~] = icemodel.test.helpers.runSmbModel(opts);
      initialization_density = ro_ice;
      if use_ro_glc
         initialization_density = ro_glc;
      end
      expected_initial_fraction = ...
         opts.ro_ice_init / initialization_density;
      expected_solid_mwe = ro_ice / ro_liq ...
         * expected_initial_fraction * sum(dz);
      solid_start(k) = ice1.mass_budget_solid_start_mwe(1);
      testCase.verifyEqual( ...
         solid_start(k), expected_solid_mwe, ...
         'AbsTol', 1e-10);

      % Physical phase change and mixed-latent-heat vapor accounting must close
      % independently of the initialization-only option.
      testCase.verifyEqual( ...
         ice1.mass_budget_phase_solid_mwe ...
         + ice1.mass_budget_phase_liquid_mwe, ...
         zeros(size(ice1.mass_budget_phase_solid_mwe)), 'AbsTol', 1e-12);
      vapor_accounted = ro_liq * ( ...
         Ls * ice1.mass_budget_vapor_solid_mwe ...
         + Lv * ice1.mass_budget_vapor_liquid_mwe ...
         + Lv * ice1.mass_budget_condensation_overflow_mwe);
      % Allow only subtraction roundoff from column-integrated checkpoints.
      testCase.verifyEqual( ...
         ice1.mass_budget_vapor_potential_j_m2, vapor_accounted, ...
         'AbsTol', 2e-5);
   end
   testCase.verifyNotEqual(solid_start(1), solid_start(2));
   testCase.verifyEqual( ...
      solid_start(2) / solid_start(1), ro_ice / ro_glc, 'AbsTol', 1e-12);
   clear cleanup
end

function test_no_merge_returns_zero_event_and_exchange_ledger(testCase)
   % An eligibility-free call must leave state, d_lyr, and the remesh budget
   % channels unchanged.

   [T, f_ice, f_liq, Sc, Sp, d_lyr] = mergeFixture([0.5; 0.6; 0.7]);
   dz = 0.04;
   budget = icemodel.column.initialize_budget_state(T, f_ice, f_liq, dz);
   [T_new, f_ice_new, f_liq_new, Sc_new, Sp_new, d_lyr_new, budget] = ...
      icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, dz, 0.0, d_lyr, 0.1, budget);

   % No inferred event or storage exchange is allowed on the early return.
   testCase.verifyEqual(T_new, T);
   testCase.verifyEqual(f_ice_new, f_ice);
   testCase.verifyEqual(f_liq_new, f_liq);
   testCase.verifyEqual(Sc_new, Sc);
   testCase.verifyEqual(Sp_new, Sp);
   testCase.verifyEqual(d_lyr_new, d_lyr);
   testCase.verifyEqual(budget.mass_budget_remesh_solid_mwe, 0, 'AbsTol', 0);
   testCase.verifyEqual(budget.mass_budget_remesh_liquid_mwe, 0, 'AbsTol', 0);
   testCase.verifyEqual( ...
      budget.mass_budget_cloned_bottom_solid_mwe, 0, 'AbsTol', 0);
   testCase.verifyEqual( ...
      budget.mass_budget_merge_export_solid_mwe, 0, 'AbsTol', 0);
   testCase.verifyEqual(budget.mass_budget_top_deletion_count, 0, 'AbsTol', 0);
   testCase.verifyEqual( ...
      budget.mass_budget_top_deletion_height_m, 0, 'AbsTol', 0);
   testCase.verifyEqual( ...
      budget.mass_budget_top_export_solid_mwe, 0, 'AbsTol', 0);
   testCase.verifyEqual( ...
      budget.mass_budget_top_export_liquid_mwe, 0, 'AbsTol', 0);
   testCase.verifyEqual( ...
      budget.mass_budget_interior_merge_count, 0, 'AbsTol', 0);
end

function test_top_export_is_the_mass_a_surface_removal_removes(testCase)
   % A merge retains the mean of the two combined cells, so the mass leaving
   % the column is their combined mass minus what the merged cell keeps. That
   % quantity, not the quantized cell height, is the surface-loss comparator.

   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   f_ice = [0.05; 0.6; 0.7];
   [T, ~, f_liq, Sc, Sp, d_lyr] = mergeFixture(f_ice);
   dz = 0.04;
   budget = icemodel.column.initialize_budget_state(T, f_ice, f_liq, dz);
   [~, f_ice_new, f_liq_new, ~, ~, ~, budget] = ...
      icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, dz, 0.0, d_lyr, 0.1, budget);

   % Hand-compute the removed mass from the pair and the surviving merged cell.
   pair_solid = ro_ice / ro_liq * (f_ice(1) + f_ice(2)) * dz;
   retained_solid = ro_ice / ro_liq * f_ice_new(1) * dz;
   pair_liquid = (f_liq(1) + f_liq(2)) * dz;
   retained_liquid = f_liq_new(1) * dz;
   testCase.verifyEqual(budget.mass_budget_top_export_solid_mwe, ...
      pair_solid - retained_solid, 'AbsTol', 1e-15);
   testCase.verifyEqual(budget.mass_budget_top_export_liquid_mwe, ...
      pair_liquid - retained_liquid, 'AbsTol', 1e-15);

   % A top removal exports mass and also translates the grid, but the two are
   % different quantities and must not be interconverted.
   testCase.verifyEqual(budget.mass_budget_top_deletion_count, 1);
   testCase.verifyGreaterThan(budget.mass_budget_top_export_solid_mwe, 0);
   testCase.verifyNotEqual(budget.mass_budget_top_export_solid_mwe, ...
      budget.mass_budget_top_deletion_height_m);
end

function test_runoff_credits_condensation_overflow(testCase)
   % Condensation the top cell could not store never entered the reservoir, so
   % it runs off directly. Dropping it would take that water out of the budget
   % entirely.

   opts = struct('dz_thermal', 0.04, 'tlag', 2);
   n_steps = 10;
   ice2 = struct('df_liq', zeros(3, n_steps), 'df_evp', zeros(3, n_steps));
   ice1 = struct('df_rof', zeros(n_steps, 1));

   base = icemodel.column.diagnose_column_runoff(ice1, ice2, opts);
   ice1.df_rof(5) = 0.25;
   returned = icemodel.column.diagnose_column_runoff(ice1, ice2, opts);

   % df_rof is a top-cell liquid fraction, so it scales by the cell thickness.
   expected = 0.25 * opts.dz_thermal;
   testCase.verifyEqual(base.runoff(end), 0, 'AbsTol', 1e-15);
   testCase.verifyEqual(returned.runoff(end) - base.runoff(end), expected, ...
      'AbsTol', 1e-15);
end

function test_runoff_subtracts_evaporated_pore_water(testCase)
   % Evaporation removes liquid that runoff would otherwise have carried, so
   % counting it as runoff would overstate the water that actually drained.

   opts = struct('dz_thermal', 0.04, 'tlag', 2);
   n_steps = 10;
   ice2 = struct('df_liq', zeros(3, n_steps), 'df_evp', zeros(3, n_steps));
   ice1 = struct('df_rof', zeros(n_steps, 1));
   ice2.df_liq(1, 3) = 0.5;

   melt_only = icemodel.column.diagnose_column_runoff(ice1, ice2, opts);
   ice2.df_evp(1, 6) = -0.2;
   returned = icemodel.column.diagnose_column_runoff(ice1, ice2, opts);

   expected = -0.2 * opts.dz_thermal;
   testCase.verifyEqual(melt_only.runoff(end), 0.5 * opts.dz_thermal, ...
      'AbsTol', 1e-15);
   testCase.verifyEqual(returned.runoff(end) - melt_only.runoff(end), ...
      expected, 'AbsTol', 1e-15);

   % Condensation is the opposite sign: it adds liquid that later drains.
   ice2.df_evp(1, 6) = 0.2;
   condensed = icemodel.column.diagnose_column_runoff(ice1, ice2, opts);
   testCase.verifyEqual(condensed.runoff(end) - melt_only.runoff(end), ...
      -expected, 'AbsTol', 1e-15);
end

function test_runoff_is_independent_of_output_profile(testCase)
   % Runoff is a physical diagnostic, so every icemodel profile must carry the
   % channels it consumes. Gating df_rof behind the diagnostic profile would
   % have made the same run report different runoff depending on its outputs.

   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=8, dt_seconds=3600);
   cleanup = onCleanup(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));

   for output_profile = ["minimal", "standard", "diagnostic"]
      opts = icemodel.test.helpers.buildSyntheticOpts(workspace, ...
         'icemodel', 2016, output_profile=char(output_profile), solver=1);
      testCase.verifyTrue(ismember('df_rof', opts.vars1), ...
         sprintf('%s profile must carry df_rof', output_profile));
   end
   clear cleanup
end

function test_d_lyr_carries_total_merge_export_not_liquid_only(testCase)
   % d_lyr must accumulate the full water-equivalent mass a merge removes.
   % df_lyr carries solid plus liquid, which is why the derived dlayer series
   % cannot be reconciled against melt and runoff.

   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   f_ice = [0.05; 0.6; 0.7];
   [T, ~, f_liq, Sc, Sp, d_lyr] = mergeFixture(f_ice);
   dz = 0.04;
   budget = icemodel.column.initialize_budget_state(T, f_ice, f_liq, dz);
   [~, f_ice_new, f_liq_new, ~, ~, d_lyr_new, budget] = ...
      icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, dz, 0.0, d_lyr, 0.1, budget);

   % Hand-compute the removed solid and liquid mass from the pair and the
   % surviving merged cell, independent of the budget channels, so the
   % comparator does not read back its own oracle.
   pair_solid = ro_ice / ro_liq * (f_ice(1) + f_ice(2)) * dz;
   retained_solid = ro_ice / ro_liq * f_ice_new(1) * dz;
   pair_liquid = (f_liq(1) + f_liq(2)) * dz;
   retained_liquid = f_liq_new(1) * dz;
   expected_export = (pair_solid - retained_solid) + (pair_liquid - retained_liquid);

   % Scaling by the cell thickness must reproduce the independently computed
   % total export, which combines solid and liquid mass removed.
   exported_mwe = sum(d_lyr_new) * dz;
   testCase.verifyEqual(exported_mwe, expected_export, 'AbsTol', 1e-15);

   % The liquid difference alone is strictly smaller than the total export,
   % so a liquid-only accounting cannot satisfy the identity above.
   liquid_only = pair_liquid - retained_liquid;
   testCase.verifyLessThan(liquid_only, exported_mwe);

   % The top-removal solid export is a subset of the all-merge solid export.
   % The liquid side has no matching nesting check: the schema carries only
   % the net remesh_liquid_mwe and this fixture's top_export_liquid_mwe, so
   % there is no separate all-merge liquid channel to nest
   % against.
   testCase.verifyLessThanOrEqual(budget.mass_budget_top_export_solid_mwe, ...
      budget.mass_budget_merge_export_solid_mwe + 1e-15);
end

function test_interior_merge_exports_no_surface_mass(testCase)
   % Interior remeshing moves mass without lowering the surface, so it must
   % leave the surface-loss comparator untouched.

   [T, f_ice, f_liq, Sc, Sp, d_lyr] = ...
      mergeFixture([0.6; 0.05; 0.7; 0.8]);
   dz = 0.04;
   budget = icemodel.column.initialize_budget_state(T, f_ice, f_liq, dz);
   [~, ~, ~, ~, ~, ~, budget] = icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, dz, 0.0, d_lyr, 0.1, budget);

   testCase.verifyEqual(budget.mass_budget_interior_merge_count, 1);
   testCase.verifyEqual(budget.mass_budget_top_deletion_count, 0);
   testCase.verifyEqual(budget.mass_budget_top_export_solid_mwe, 0, 'AbsTol', 0);
   testCase.verifyEqual( ...
      budget.mass_budget_top_export_liquid_mwe, 0, 'AbsTol', 0);

   % The interior event still contributes to the total merge export, so the
   % two channels are genuinely distinct rather than duplicates.
   testCase.verifyGreaterThan(budget.mass_budget_merge_export_solid_mwe, 0);
end

function test_top_merge_counts_actual_grid_translation(testCase)
   % Removing the actual top cell contributes exactly one uniform-grid dz.

   [T, f_ice, f_liq, Sc, Sp, d_lyr] = mergeFixture([0.05; 0.6; 0.7]);
   dz = 0.04;
   budget = icemodel.column.initialize_budget_state(T, f_ice, f_liq, dz);
   [~, ~, ~, ~, ~, d_lyr_new, budget] = icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, dz, 0.0, d_lyr, 0.1, budget);

   % Geometry and the legacy non-geometric diagnostic remain distinct.
   testCase.verifyEqual(budget.mass_budget_top_deletion_count, 1);
   testCase.verifyEqual( ...
      budget.mass_budget_top_deletion_height_m, 0.04, 'AbsTol', 0);
   testCase.verifyEqual(budget.mass_budget_interior_merge_count, 0);
   testCase.verifyGreaterThan(sum(d_lyr_new), 0);
   verifyRemeshIdentity(testCase, budget);
end

function test_interior_merge_does_not_translate_top_grid(testCase)
   % Removing an interior cell is numerical remeshing, not surface lowering.

   [T, f_ice, f_liq, Sc, Sp, d_lyr] = ...
      mergeFixture([0.6; 0.05; 0.7; 0.8]);
   dz = 0.04;
   budget = icemodel.column.initialize_budget_state(T, f_ice, f_liq, dz);
   [~, ~, ~, ~, ~, ~, budget] = icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, dz, 0.0, d_lyr, 0.1, budget);

   % The actual-event ledger must not infer height from the eligibility mask.
   testCase.verifyEqual(budget.mass_budget_top_deletion_count, 0);
   testCase.verifyEqual(budget.mass_budget_top_deletion_height_m, 0, 'AbsTol', 0);
   testCase.verifyEqual(budget.mass_budget_interior_merge_count, 1);
   verifyRemeshIdentity(testCase, budget);
end

function test_multiple_top_merges_follow_index_drift_once_per_original_flag(testCase)
   % Adjacent flagged top cells shift to index one as each prior cell is removed.

   [T, f_ice, f_liq, Sc, Sp, d_lyr] = ...
      mergeFixture([0.05; 0.04; 0.6; 0.7]);
   dz = 0.04;
   budget = icemodel.column.initialize_budget_state(T, f_ice, f_liq, dz);
   [~, ~, ~, ~, ~, ~, budget] = icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, dz, 0.0, d_lyr, 0.1, budget);

   % Appended bottom clones are not reprocessed as new eligible cells.
   testCase.verifyEqual(budget.mass_budget_top_deletion_count, 2);
   testCase.verifyEqual( ...
      budget.mass_budget_top_deletion_height_m, 0.08, 'AbsTol', 0);
   testCase.verifyEqual(budget.mass_budget_interior_merge_count, 0);
   verifyRemeshIdentity(testCase, budget);
end

function test_full_diagnostic_grid_staircase_contracts_at_two_resolutions(testCase)
   % The accepted resolution evidence must traverse production accumulation and
   % postprocessing, not stop at the remeshing kernel that creates each event.
   target_height = 0.26;
   dz = [0.10, 0.04];
   grid_height = zeros(size(dz));
   discrepancy = zeros(size(dz));
   for k = 1:numel(dz)
      [ledger, n_top] = runGridResolutionFixture(dz(k), target_height);

      % Postprocessed interval ledgers preserve actual top and interior event
      % provenance and the exact uniform-grid height conversion.
      top_count = sum(ledger.mass_budget_top_deletion_count);
      interior_count = sum(ledger.mass_budget_interior_merge_count);
      grid_height(k) = sum(ledger.mass_budget_top_deletion_height_m);
      testCase.verifyEqual(top_count, n_top);
      testCase.verifyEqual(interior_count, 1);
      testCase.verifyEqual(ledger.mass_budget_top_deletion_height_m, ...
         dz(k) * ledger.mass_budget_top_deletion_count, 'AbsTol', 1e-15);
      testCase.verifyEqual(grid_height(k), n_top * dz(k), 'AbsTol', 1e-15);

      % R = B - O must still close after forcing-step accumulation and the
      % budget-aware hourly retime used by downstream comparison and
      % reporting. Only the solid side carries this cloned/exported
      % decomposition; the liquid side carries only the net
      % mass_budget_remesh_liquid_mwe.
      testCase.verifyEqual(ledger.mass_budget_remesh_solid_mwe, ...
         ledger.mass_budget_cloned_bottom_solid_mwe ...
         - ledger.mass_budget_merge_export_solid_mwe, 'AbsTol', 1e-15);

      discrepancy(k) = abs(target_height - grid_height(k));
      testCase.verifyLessThanOrEqual(discrepancy(k), dz(k));
   end

   testCase.verifyLessThan(discrepancy(2), discrepancy(1));
end

function test_depleted_bottom_removal_matches_minimal_transition(testCase)
   % Requesting the budget must never change the state transition, including
   % when the deepest cell is itself the merge-eligible layer.

   [T, f_ice, f_liq, Sc, Sp, d_lyr] = ...
      mergeFixture([0.6; 0.7; 0.05]);
   dz = 0.04;
   budget = icemodel.column.initialize_budget_state(T, f_ice, f_liq, dz);
   [~, f_ice_new, ~, ~, ~, ~, budget] = icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, dz, 0.0, d_lyr, 0.1, budget);

   % The deepest cell is removed, then the surviving bottom state is cloned.
   % The removed cell must not come back: cloning before the deletion copied
   % it into the column, so it stayed while the cell above it kept only half
   % the pair's mass.
   testCase.verifyNotEqual(f_ice_new(end), f_ice(end));
   testCase.verifyEqual(f_ice_new(end), f_ice_new(end - 1), 'AbsTol', 0);
   testCase.verifyGreaterThanOrEqual(min(f_ice_new), 0.1);

   % A deepest-cell removal is not a surface removal, so it records interior
   % merge activity and contributes no quantized grid translation.
   testCase.verifyEqual(budget.mass_budget_top_deletion_count, 0);
   testCase.verifyEqual(budget.mass_budget_interior_merge_count, 1);
   testCase.verifyEqual(budget.mass_budget_top_deletion_height_m, 0.0);
   verifyRemeshIdentity(testCase, budget);
end

function test_bottom_adjacent_merge_matches_minimal_transition(testCase)
   % A bottom-adjacent merge clones an already-depleted reservoir, and the
   % ledger must still close its domain-exchange identity.

   [T, f_ice, f_liq, Sc, Sp, d_lyr] = mergeFixture([0.6; 0.0; 0.15]);
   f_liq(:) = 0.0;
   [j1, j2] = icemodel.column.merge_layer_indices(2, f_ice);
   testCase.verifyEqual([j1, j2], [2, 3]);
   testCase.verifyGreaterThan(f_ice(end), 0.1);

   dz = 0.04;
   budget = icemodel.column.initialize_budget_state(T, f_ice, f_liq, dz);
   [~, f_ice_new, ~, ~, ~, ~, budget] = icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, dz, 0.0, d_lyr, 0.1, budget);

   testCase.verifyLessThanOrEqual(f_ice_new(end), 0.1);
   testCase.verifyEqual(budget.mass_budget_top_deletion_count, 0);
   testCase.verifyEqual(budget.mass_budget_interior_merge_count, 1);
   verifyRemeshIdentity(testCase, budget);
end

function test_merge_prediction_uses_physical_vapor_basis(testCase)
   % Merge prediction and d_pevp must share the physical-density Lv-to-Ls
   % conversion used by the production surface-vapor kernel.

   [Ls, Lv, ro_ice, ro_liq] = ...
      icemodel.physicalConstant('Ls', 'Lv', 'ro_ice', 'ro_liq');
   d_pevp = -0.02;
   potential_ice_change = d_pevp * (Lv * ro_liq) / (Ls * ro_ice);
   f_ice_min = 0.1;
   f_top = f_ice_min - potential_ice_change / 2;
   [T, f_ice, f_liq, Sc, Sp, d_lyr] = ...
      mergeFixture([f_top; 0.6; 0.7]);

   % The current top is retained, but the physical vapor prediction crosses the
   % floor and therefore marks only that top cell for merging.
   testCase.verifyGreaterThan(f_ice(1), f_ice_min);
   dz = 0.04;
   budget = icemodel.column.initialize_budget_state(T, f_ice, f_liq, dz);
   [~, ~, ~, ~, ~, ~, budget] = icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, dz, d_pevp, d_lyr, f_ice_min, budget);
   testCase.verifyEqual(budget.mass_budget_top_deletion_count, 1);
   testCase.verifyEqual(budget.mass_budget_interior_merge_count, 0);
end

function test_potential_sublimation_converts_on_latent_heat(testCase)
   % One helper owns the liquid-to-ice conversion. Compare it against the
   % conversion written out from the constants, not against itself.

   [Ls, Lv, ro_ice, ro_liq] = ...
      icemodel.physicalConstant('Ls', 'Lv', 'ro_ice', 'ro_liq');

   % Compare bit for bit, and write the conversion with the same association
   % the callers use. Multiplying by a cached factor instead of dividing
   % after the multiply differs in the last bit. That would move every
   % default-mode result for no physical reason.
   %
   % Sublimation and deposition must both pass, because the callers apply the
   % helper to a signed tendency and the sign carries the direction.
   for d_pevp = [-0.02, 0, 0.03]
      returned = icemodel.column.potential_sublimation(d_pevp);
      expected = d_pevp * (Lv * ro_liq) / (Ls * ro_ice);
      testCase.verifyEqual(returned, expected, 'AbsTol', 0);
   end

   % A column of tendencies must convert elementwise, because the coupled
   % vapor path applies the same conversion to every cell.
   d_pevp = [-0.02; 0.01];
   testCase.verifyEqual( ...
      icemodel.column.potential_sublimation(d_pevp), ...
      d_pevp * (Lv * ro_liq) / (Ls * ro_ice), 'AbsTol', 0);
end

function test_merge_prediction_reaches_the_top_layer_only(testCase)
   % Only the top layer receives d_pevp, so the merge look-ahead must not flag
   % near-threshold interior layers. Broadcasting the scalar across the column
   % deletes interior layers because of a mass change they never receive.

   [Ls, Lv, ro_ice, ro_liq] = ...
      icemodel.physicalConstant('Ls', 'Lv', 'ro_ice', 'ro_liq');
   d_pevp = -0.02;
   potential_ice_change = d_pevp * (Lv * ro_liq) / (Ls * ro_ice);
   f_ice_min = 0.1;

   % Place all three layers the same small distance above the floor, at half
   % the predicted ice loss. Every layer therefore crosses the floor under the
   % broadcast rule, and only the top layer crosses under the top-only rule.
   f_near = f_ice_min - potential_ice_change / 2;
   [T, f_ice, f_liq, Sc, Sp, d_lyr] = ...
      mergeFixture([f_near; f_near; f_near]);
   testCase.verifyGreaterThan(min(f_ice), f_ice_min);

   dz = 0.04;
   budget = icemodel.column.initialize_budget_state(T, f_ice, f_liq, dz);
   [~, ~, ~, ~, ~, ~, budget] = icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, dz, d_pevp, d_lyr, f_ice_min, budget);
   testCase.verifyEqual(budget.mass_budget_top_deletion_count, 1);
   testCase.verifyEqual(budget.mass_budget_interior_merge_count, 0);

   % A zero tendency must leave every layer unflagged, which shows the counts
   % above came from the prediction and not from the floor test.
   budget_no_vapor = icemodel.column.initialize_budget_state(T, f_ice, f_liq, dz);
   [~, ~, ~, ~, ~, ~, budget_no_vapor] = icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, dz, 0.0, d_lyr, f_ice_min, budget_no_vapor);
   testCase.verifyEqual(budget_no_vapor.mass_budget_top_deletion_count, 0);
   testCase.verifyEqual(budget_no_vapor.mass_budget_interior_merge_count, 0);
end

function test_pending_flags_match_legacy_multi_event_result(testCase)
   % A non-top pending mask must retain the historical appended bottom flag so
   % state and d_lyr evolution stay exactly unchanged. Requesting the budget
   % must never alter the transition it observes.

   [T, f_ice, f_liq, Sc, Sp, d_lyr] = ...
      mergeFixture([0.6; 0.05; 0.04]);
   dz = 0.04;
   budget = icemodel.column.initialize_budget_state(T, f_ice, f_liq, dz);
   [T_new, f_ice_new, f_liq_new, Sc_new, Sp_new, d_lyr_new, budget] = ...
      icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, dz, 0.0, d_lyr, 0.1, budget);
   [T_legacy, f_ice_legacy, f_liq_legacy, Sc_legacy, Sp_legacy, ...
      d_lyr_legacy, legacy_mask] = legacyMergeThinLayers( ...
      T, f_ice, f_liq, Sc, Sp, dz, 0.0, d_lyr, 0.1);

   % Compare state against the legacy transition oracle, including its mask.
   testCase.verifyEqual(legacy_mask, [false; true; true]);
   testCase.verifyEqual(T_new, T_legacy);
   testCase.verifyEqual(f_ice_new, f_ice_legacy);
   testCase.verifyEqual(f_liq_new, f_liq_legacy);
   testCase.verifyEqual(Sc_new, Sc_legacy);
   testCase.verifyEqual(Sp_new, Sp_legacy);
   testCase.verifyEqual(d_lyr_new, d_lyr_legacy);

   % Neither flagged cell is the original top cell, so the two-event merge
   % records only interior activity and no grid translation.
   testCase.verifyEqual(budget.mass_budget_top_deletion_count, 0);
   testCase.verifyEqual(budget.mass_budget_interior_merge_count, 2);
end

function test_legacy_parity_holds_with_a_nonzero_vapor_tendency(testCase)
   % Every other legacy-parity call passes d_pevp = 0, which leaves the
   % oracle's vapor term unexercised: a production rule change would not move
   % the comparison. This case drives the term on both sides, so the oracle
   % fails if the two eligibility rules diverge.

   [Ls, Lv, ro_ice, ro_liq] = ...
      icemodel.physicalConstant('Ls', 'Lv', 'ro_ice', 'ro_liq');
   d_pevp = -0.02;
   potential_ice_change = d_pevp * (Lv * ro_liq) / (Ls * ro_ice);
   f_ice_min = 0.1;

   % The top layer crosses the floor through the vapor term alone. The
   % second layer is already below the floor. The transition therefore does a
   % real multi-event merge, not a comparison of two empty results.
   f_top = f_ice_min - potential_ice_change / 2;
   [T, f_ice, f_liq, Sc, Sp, d_lyr] = mergeFixture([f_top; 0.05; 0.7]);
   dz = 0.04;
   budget = icemodel.column.initialize_budget_state(T, f_ice, f_liq, dz);
   [T_new, f_ice_new, f_liq_new, Sc_new, Sp_new, d_lyr_new, ~] = ...
      icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, dz, d_pevp, d_lyr, f_ice_min, budget);
   [T_legacy, f_ice_legacy, f_liq_legacy, Sc_legacy, Sp_legacy, ...
      d_lyr_legacy, legacy_mask] = legacyMergeThinLayers( ...
      T, f_ice, f_liq, Sc, Sp, dz, d_pevp, d_lyr, f_ice_min);

   % The vapor term must be what flags the top layer, or the case would not
   % reach the branch it exists to cover.
   testCase.verifyEqual(legacy_mask, [true; true; false]);
   testCase.verifyEqual(T_new, T_legacy);
   testCase.verifyEqual(f_ice_new, f_ice_legacy);
   testCase.verifyEqual(f_liq_new, f_liq_legacy);
   testCase.verifyEqual(Sc_new, Sc_legacy);
   testCase.verifyEqual(Sp_new, Sp_legacy);
   testCase.verifyEqual(d_lyr_new, d_lyr_legacy);
end

function test_vapor_identity_partitions_wet_evaporation_and_sublimation(testCase)
   % Evaporation that exhausts mobile liquid must use Lv for liquid and Ls for
   % the remaining realized solid sublimation.

   d_rof = verifyVaporIdentity(testCase, 0.5, 0.05, -0.05, 0.1, 0.02);
   testCase.verifyEqual(d_rof, 0, 'AbsTol', 0);
end

function test_vapor_identity_includes_condensation_overflow(testCase)
   % Wet condensation beyond cell capacity must be exposed without routing it.

   d_rof = verifyVaporIdentity(testCase, 0.99, 0.005, 0.01, 0.1, 0.02);
   testCase.verifyGreaterThan(d_rof, 0);
end

function test_vapor_identity_cascades_to_the_next_cell(testCase)
   % A top cell too thin to supply the full sublimation demand draws the
   % remainder from the cell below it, cell by cell, until the demand is
   % satisfied. No remainder is left unaccounted, so the vapor-energy
   % identity closes exactly across both cells.

   [Tf, Ls, Lv, ro_liq] = ...
      icemodel.physicalConstant('Tf', 'Ls', 'Lv', 'ro_liq');
   f_ice_min = 0.1;
   f_res_por = 0.02;
   dz = 0.04 * ones(2, 1);
   T = (Tf - 2) * ones(2, 1);
   f_ice = [f_ice_min + 1e-4; 0.8];
   f_liq = zeros(2, 1);

   % Demand far larger than what cell 1 can give above the retained floor.
   d_pevp = -0.05;

   budget = icemodel.column.initialize_budget_state(T, f_ice, f_liq, dz);
   [~, f_ice_new, ~, ~, ~, d_rof, ~, ~, ~, budget] = ...
      icemodel.column.budget_surface_mass_balance( ...
      T, f_ice, f_liq, f_liq, d_pevp, zeros(2, 1), zeros(2, 1), 0, ...
      zeros(2, 1), zeros(2, 1), f_res_por, f_ice_min, budget, dz);

   % Cell 1 stops at the floor; cell 2 supplies the rest of the demand.
   testCase.verifyEqual(f_ice_new(1), f_ice_min, 'AbsTol', 1e-12);
   testCase.verifyLessThan(f_ice_new(2), f_ice(2));
   testCase.verifyEqual(d_rof, 0, 'AbsTol', 0);

   % The demanded energy and the energy the cascade actually moved match
   % exactly: no cell is left holding an unaccounted remainder.
   potential = budget.mass_budget_vapor_potential_j_m2;
   accepted = ro_liq * ( ...
      Ls * budget.mass_budget_vapor_solid_mwe ...
      + Lv * budget.mass_budget_vapor_liquid_mwe ...
      + Lv * budget.mass_budget_condensation_overflow_mwe);
   testCase.verifyEqual(potential, accepted, 'AbsTol', 1e-7);
end

function test_vapor_identity_accepts_dry_sublimation(testCase)
   % Dry sublimation must close on the physical ro_liq/Lv to ro_ice/Ls basis.

   d_rof = verifyVaporIdentity(testCase, 0.5, 0.0, -0.02, 0.1, 0.02);
   testCase.verifyEqual(d_rof, 0, 'AbsTol', 0);
end

function test_vapor_identity_accepts_dry_deposition(testCase)
   % Dry deposition must add solid mass on the same physical latent-heat basis.

   d_rof = verifyVaporIdentity(testCase, 0.5, 0.0, 0.02, 0.1, 0.02);
   testCase.verifyEqual(d_rof, 0, 'AbsTol', 0);
end

function test_hourly_retime_uses_budget_aggregation_classes(testCase)
   % State endpoints select first/last while signed increments sum.

   time = datetime(2020, 1, 1, 'TimeZone', 'UTC') ...
      + minutes(15) * transpose(0:7);
   signed_phase = repmat([1; -1], 4, 1);
   positive_vapor = abs(signed_phase);
   TT = timetable(time, transpose(1:8), transpose(10:17), ...
      transpose(20:27), signed_phase, positive_vapor, ...
      'VariableNames', {'ordinary_mean', ...
      'mass_budget_solid_start_mwe', 'mass_budget_solid_end_mwe', ...
      'mass_budget_phase_solid_mwe', ...
      'mass_budget_vapor_solid_mwe'});

   TT = icemodel.retimeHourlyFixedStep(TT);

   % The default mean remains intact for every unrelated output channel.
   testCase.verifyEqual(TT.ordinary_mean, [2.5; 6.5], 'AbsTol', 0);
   testCase.verifyEqual(TT.mass_budget_solid_start_mwe, [10; 14]);
   testCase.verifyEqual(TT.mass_budget_solid_end_mwe, [23; 27]);
   % A signed sum-class channel sums to zero on this alternating pattern.
   testCase.verifyEqual(TT.mass_budget_phase_solid_mwe, [0; 0]);
   % A second, independently named sum-class channel is also summed, which
   % shows the aggregation rule applies across the whole 'sum' class rather
   % than to one hardcoded field.
   testCase.verifyEqual(TT.mass_budget_vapor_solid_mwe, [4; 4]);
end

function test_hourly_retime_preserves_variable_classes_and_double_precision(testCase)
   % A single-valued flag must not down-cast double scientific diagnostics.

   time = datetime(2020, 1, 1, 'TimeZone', 'UTC') ...
      + minutes(15) * transpose(0:3);
   precise = 1 + 1e-10 * transpose(1:4);
   TT = timetable(time, precise, single(ones(4, 1)), ...
      'VariableNames', {'ordinary_double', 'single_flag'});

   TT = icemodel.retimeHourlyFixedStep(TT);

   % Preserve each class independently and retain information below single eps.
   testCase.verifyClass(TT.ordinary_double, 'double');
   testCase.verifyClass(TT.single_flag, 'single');
   testCase.verifyEqual(TT.ordinary_double, mean(precise), 'AbsTol', 1e-15);
   testCase.verifyGreaterThan(TT.ordinary_double, 1);
end

function test_hourly_retime_accepts_unaligned_short_partial_bin(testCase)
   % Valid 15-minute windows need not start on the hour or contain four rows.

   time = datetime(2020, 1, 1, 0, 15, 0, 'TimeZone', 'UTC') ...
      + minutes(15) * transpose(0:2);
   TT = timetable(time, [1; 2; 3], [10; 11; 12], [20; 21; 22], ...
      [1; -1; 2], [1; 1; 2], ...
      'VariableNames', {'ordinary_mean', ...
      'mass_budget_solid_start_mwe', 'mass_budget_solid_end_mwe', ...
      'mass_budget_phase_solid_mwe', ...
      'mass_budget_vapor_solid_mwe'});

   TT = icemodel.retimeHourlyFixedStep(TT);

   % Native hourly bins keep endpoint values and sum without cancelling.
   testCase.verifyEqual(TT.Properties.RowTimes, ...
      datetime(2020, 1, 1, 'TimeZone', 'UTC'));
   testCase.verifyEqual(TT.ordinary_mean, 2);
   testCase.verifyEqual(TT.mass_budget_solid_start_mwe, 10);
   testCase.verifyEqual(TT.mass_budget_solid_end_mwe, 22);
   testCase.verifyEqual(TT.mass_budget_phase_solid_mwe, 2);
   testCase.verifyEqual(TT.mass_budget_vapor_solid_mwe, 4);
end

function test_non_diagnostic_model_paths_omit_budget_channels(testCase)
   % Standard and minimal production runs must omit diagnostic bookkeeping
   % channels and preserve their externally visible profile schemas, even
   % though budget accumulation itself is unconditional: only emission is
   % profile-bound. The profiler proves the accumulators still ran on every
   % accepted substep even though the diagnostic-only emission step never
   % copies their channels into ice1.

   profiles = {'standard', 'minimal'};
   for n = 1:numel(profiles)
      workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
         configure=true, nsteps=4, dt_seconds=3600);
      cleanup = onCleanup(@() ...
         icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));
      opts = icemodel.test.helpers.buildSyntheticOpts( ...
         workspace, 'icemodel', 2016, output_profile=profiles{n}, solver=1);

      % The profiler status is saved and restored because a test must not
      % clear a profiling run its caller started.
      prior = profile('status');
      restore_profiler = onCleanup(@() restoreProfiler(prior));
      profile off
      profile clear
      profile on
      [ice1, ~, run_opts] = icemodel.test.helpers.runSmbModel(opts);
      profile off
      called = string({profile('info').FunctionTable.FunctionName});
      testCase.verifyTrue( ...
         any(contains(called, 'accumulate_phase_budget')), profiles{n});
      testCase.verifyTrue( ...
         any(contains(called, 'budget_surface_mass_balance')), profiles{n});
      clear restore_profiler

      % Exercise each non-diagnostic branch rather than only its option list.
      testCase.verifyEqual(transpose(fieldnames(ice1)), run_opts.vars1, ...
         profiles{n});
      testCase.verifyFalse(any(ismember(fieldnames(ice1), ...
         icemodel.namelists.budgetoutputs())), profiles{n});
      clear cleanup
   end
end

function [ledger, n_top] = runGridResolutionFixture(dz_thermal, target_height)
   %RUNGRIDRESOLUTIONFIXTURE Run one restart-driven diagnostic model interval.

   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=1, dt_seconds=3600);
   cleanup = onCleanup(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'icemodel', 2016, output_profile='diagnostic', solver=0);
   opts = icemodel.resetopts(opts, ...
      'dz_thermal', dz_thermal, 'z0_thermal', 1.2);

   % Seed known low-density top cells plus one separated interior cell. The
   % restart path enters the same solver, accumulation, and remeshing branches as
   % a normal diagnostic run while keeping the synthetic event count exact.
   [mesh_dz, ~, ~, ~, ~] = icemodel.column.control_volume_mesh( ...
      opts.z0_thermal, opts.dz_thermal);
   n_top = floor(target_height / dz_thermal);
   n_cells = numel(mesh_dz);
   Tf = icemodel.physicalConstant('Tf');
   restart = struct();
   restart.T = (Tf - 5) * ones(n_cells, 1);
   restart.f_ice = 0.7 * ones(n_cells, 1);
   restart.f_ice(1:n_top) = 0.05;
   restart.f_ice(n_top + 2) = 0.05;
   restart.f_liq = zeros(n_cells, 1);
   restart.Ts = Tf - 5;
   restart.r_eff = 1e-3 * ones(n_cells, 1);
   restart_file = fullfile(workspace.rootdir, 'grid-resolution-restart.mat');
   save(restart_file, 'restart')

   % Reconfigure after the mesh and restart overrides, then pass the raw output
   % through the same hourly postprocessor consumed by the evaluation runner.
   opts = icemodel.resetopts(opts, ...
      'use_restart', true, 'restartfile', restart_file);
   opts = icemodel.configureRun(opts);
   [ice1, ice2, opts] = icemodel.test.helpers.runSmbModel(opts);
   [ledger, ~] = icemodel.postprocess( ...
      ice1, ice2, opts, opts.output_years);
   clear cleanup
end

function [T, f_ice, f_liq, Sc, Sp, d_lyr] = mergeFixture(f_ice)
   %MERGEFIXTURE Return a compact cold-column remesh fixture.

   Tf = icemodel.physicalConstant('Tf');
   T = (Tf - 2) * ones(size(f_ice));
   f_liq = 0.01 * ones(size(f_ice));
   Sc = zeros(size(f_ice));
   Sp = zeros(size(f_ice));
   d_lyr = zeros(size(f_ice));
end

function verifyRemeshIdentity(testCase, budget)
   %VERIFYREMESHIDENTITY Check R = B - O for the accumulated remesh budget.
   %
   % Only the solid side carries the cloned/exported decomposition:
   % remesh_solid = cloned_bottom_solid - merge_export_solid. The 21-channel
   % schema has no liquid or enthalpy decomposition and no gross
   % accumulator.

   testCase.verifyEqual(budget.mass_budget_remesh_solid_mwe, ...
      budget.mass_budget_cloned_bottom_solid_mwe ...
      - budget.mass_budget_merge_export_solid_mwe, 'AbsTol', 1e-15);
end

function d_rof = verifyVaporIdentity( ...
      testCase, f_ice, f_liq, d_pevp, f_ice_min, f_res_por)
   %VERIFYVAPORIDENTITY Check the demanded energy equals the realized energy.
   %
   % A single cell with enough ice and liquid above its floors satisfies
   % the whole demand, so the vapor-energy identity closes exactly: no
   % remainder is drawn from a cell below and none is left unaccounted.

   [Tf, Ls, Lv, ro_liq] = icemodel.physicalConstant('Tf', 'Ls', 'Lv', 'ro_liq');
   dz = 0.04;
   T = Tf - 2;

   budget = icemodel.column.initialize_budget_state(T, f_ice, f_liq, dz);
   [~, ~, ~, ~, ~, d_rof, ~, ~, ~, budget] = ...
      icemodel.column.budget_surface_mass_balance( ...
      T, f_ice, f_liq, f_liq, d_pevp, 0, 0, 0, 0, 0, ...
      f_res_por, f_ice_min, budget, dz);

   % Potential = realized solid + realized liquid + overflow, all expressed
   % with the physical solver densities used by d_pevp and sublimation.
   potential = budget.mass_budget_vapor_potential_j_m2;
   accepted = ro_liq * ( ...
      Ls * budget.mass_budget_vapor_solid_mwe ...
      + Lv * budget.mass_budget_vapor_liquid_mwe ...
      + Lv * budget.mass_budget_condensation_overflow_mwe);
   testCase.verifyEqual(potential, accepted, 'AbsTol', 1e-7);
end

function [T, f_ice, f_liq, Sc, Sp, d_lyr, merge_mask] = ...
      legacyMergeThinLayers( ...
      T, f_ice, f_liq, Sc, Sp, dz_therm, d_pevp, d_lyr, f_ice_min)
   %LEGACYMERGETHINLAYERS Reproduce the pre-ledger seven-output transition.
   %
   % This oracle isolates the ledger from the state transition, so its
   % eligibility rule must track production. The vapor prediction applies to
   % the top layer only, because that is the only layer d_pevp reaches.

   [Ls, Lv, ro_ice, ro_liq] = ...
      icemodel.physicalConstant('Ls', 'Lv', 'ro_ice', 'ro_liq');
   merge_mask = f_ice <= f_ice_min;
   merge_mask(1) = merge_mask(1) || ...
      (f_ice(1) + d_pevp * (Lv * ro_liq) / (Ls * ro_ice)) <= f_ice_min;
   do_merge = merge_mask;

   % Follow the historical index drift, clone, delete, and pending-flag order.
   ii = 0;
   for j = 1:numel(f_ice)
      ji = j + ii;
      if ~do_merge(ji)
         continue
      end
      [j1, j2] = icemodel.column.merge_layer_indices(ji, f_ice);
      [T(j2), f_ice(j2), f_liq(j2), Sc(j2), Sp(j2), d_lyr] = ...
         icemodel.column.merge_layers( ...
         T, f_ice, f_liq, Sc, Sp, j1, j2, d_lyr, dz_therm);
      T = dropCellCloneBottom(T, j1);
      Sc = dropCellCloneBottom(Sc, j1);
      Sp = dropCellCloneBottom(Sp, j1);
      f_ice = dropCellCloneBottom(f_ice, j1);
      f_liq = dropCellCloneBottom(f_liq, j1);
      do_merge = dropCellCloneBottom(do_merge, j1);
      ii = ii - 1;
   end
end

function values = dropCellCloneBottom(values, j1)
   %DROPCELLCLONEBOTTOM Remove cell j1 and clone the bottom to keep the length.
   %
   % This mirrors what merge_thin_layers does in the solver: the array length
   % is fixed, so removing a merged cell requires inserting a node at the
   % bottom. Assigning through values(:) keeps the size constant.

   values(:) = [values(1:j1 - 1); values(j1 + 1:end); values(end)];
end

function test_condensation_overflow_is_a_step_total_not_a_substep_sum(testCase)
   % d_rof is reset once per forcing step and then accumulated across
   % substeps, unlike d_pevp which is per-substep. The ledger must
   % therefore record the running total: adding it on every accepted
   % substep would inflate the overflow channel and the vapor closure
   % identity by the substep count.

   dz = 0.04;
   n_substeps = 4;
   overflow_fraction = 1e-4;
   T_ice = 273.0;
   f_ice = 0.9;
   f_liq = 0.02;

   % One overflow event on the first substep, then three quiet substeps that
   % still carry the accumulated total forward. Every substep here reports
   % zero exchange, so this test exercises the overflow channel alone, not
   % the vapor_solid/vapor_liquid storage terms.
   budget = icemodel.column.initialize_budget_state(T_ice, f_ice, f_liq, dz);
   d_rof = 0.0;
   for n = 1:n_substeps
      if n == 1
         d_rof = d_rof + overflow_fraction;
      end
      budget = icemodel.column.accumulate_vapor_exchange(budget, ...
         0.0, 0.0, 0.0, d_rof, dz);
   end

   returned = budget.mass_budget_condensation_overflow_mwe;
   expected = overflow_fraction * dz;
   testCase.verifyEqual(returned, expected, AbsTol=1e-15)
end

function restoreProfiler(prior)
   %RESTOREPROFILER Put the profiler back the way the caller had it.
   %
   % A test must not clear a profiling run the caller started. profile('info')
   % needs the data collected during the test, so the data cannot be preserved
   % as well; restoring the on/off state is the most that can be given back.

   profile off
   profile clear
   if strcmp(prior.ProfilerStatus, 'on')
      profile on
   end
end
