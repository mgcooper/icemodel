function tests = test_mass_budget_bookkeeping
   %TEST_MASS_BUDGET_BOOKKEEPING Verify physical and remesh ledger contracts.
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

function test_budget_output_fields_are_one_partitioned_contract(testCase)
   % Every diagnostic channel must belong to exactly one retime class.

   first_fields = icemodel.column.budget_output_fields('first');
   last_fields = icemodel.column.budget_output_fields('last');
   sum_fields = icemodel.column.budget_output_fields('sum');
   all_fields = icemodel.column.budget_output_fields();

   % The all-list order is the stable diagnostic-profile append order.
   testCase.verifyEqual(all_fields, [first_fields, last_fields, sum_fields]);
   testCase.verifyEqual(numel(unique(all_fields)), numel(all_fields));
   testCase.verifyEqual(numel(all_fields), 43);
   testCase.verifyError( ...
      @() icemodel.column.budget_output_fields('median'), ...
      'icemodel:column:budgetOutputFields:kind');
end

function test_budget_ledger_is_fixed_codegen_schema(testCase)
   % The kernel ledger must use one literal scalar-double layout whose order
   % remains synchronized with the MATLAB-side output registry.

   ledger = icemodel.column.initialize_budget_state();
   budget_fields = icemodel.column.budget_output_fields('all');
   ledger_fields = transpose(fieldnames(ledger));
   ledger_values = struct2cell(ledger);

   % Exact order protects profile selection, while the value checks protect
   % the scalar type and zero-reset contract used at each forcing step.
   testCase.verifyEqual(ledger_fields, budget_fields);
   testCase.verifyTrue(all(cellfun( ...
      @(value) isa(value, 'double') && isequal(size(value), [1, 1]), ...
      ledger_values)));
   testCase.verifyEqual([ledger_values{:}], zeros(1, numel(ledger_values)));

   % Keep runtime field-name construction out of the #codegen kernel. The
   % ordered field-name registry remains available only to MATLAB consumers.
   core_source = fileread(which('icemodel'));
   testCase.verifyFalse(contains(core_source, 'cell2struct'));
   testCase.verifyFalse(contains(core_source, 'budget_output_fields'));
   testCase.verifyFalse(contains(core_source, 'surface_state.('));
end

function test_diagnostic_model_payload_matches_canonical_budget_registry(testCase)
   % The runtime diagnostic payload must contain exactly the canonical ledger.

   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=4, dt_seconds=3600);
   cleanup = onCleanup(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'icemodel', 2016, output_profile='diagnostic', solver=1);

   % Exercise initialization, timestep accumulation, payload assembly, and
   % raw output storage instead of checking configuration declarations alone.
   [ice1, ~] = icemodel.test.helpers.runSmbModel(opts);
   budget_fields = icemodel.column.budget_output_fields('all');
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
         + Lv * ice1.mass_budget_condensation_overflow_mwe) ...
         + ice1.mass_budget_unapplied_vapor_j_m2;
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
   % An eligibility-free call must leave state and legacy d_lyr unchanged.

   [T, f_ice, f_liq, Sc, Sp, d_lyr] = mergeFixture([0.5; 0.6; 0.7]);
   [T_new, f_ice_new, f_liq_new, Sc_new, Sp_new, d_lyr_new, mask, diag] = ...
      icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, 0.04, 0.0, d_lyr, 0.1);

   % No inferred event or storage exchange is allowed on the early return.
   testCase.verifyFalse(any(mask));
   testCase.verifyEqual(T_new, T);
   testCase.verifyEqual(f_ice_new, f_ice);
   testCase.verifyEqual(f_liq_new, f_liq);
   testCase.verifyEqual(Sc_new, Sc);
   testCase.verifyEqual(Sp_new, Sp);
   testCase.verifyEqual(d_lyr_new, d_lyr);
   testCase.verifyEqual(cell2mat(struct2cell(diag)), zeros(23, 1));
end

function test_top_export_is_the_mass_a_surface_removal_removes(testCase)
   % A merge retains the mean of the two combined cells, so the mass leaving
   % the column is their combined mass minus what the merged cell keeps. That
   % quantity, not the quantized cell height, is the surface-loss comparator.

   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   f_ice = [0.05; 0.6; 0.7];
   [T, ~, f_liq, Sc, Sp, d_lyr] = mergeFixture(f_ice);
   [~, f_ice_new, f_liq_new, ~, ~, ~, ~, diag] = ...
      icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, 0.04, 0.0, d_lyr, 0.1);

   % Hand-compute the removed mass from the pair and the surviving merged cell.
   pair_solid = ro_ice / ro_liq * (f_ice(1) + f_ice(2)) * 0.04;
   retained_solid = ro_ice / ro_liq * f_ice_new(1) * 0.04;
   pair_liquid = (f_liq(1) + f_liq(2)) * 0.04;
   retained_liquid = f_liq_new(1) * 0.04;
   testCase.verifyEqual(diag.top_export_solid_mwe, ...
      pair_solid - retained_solid, 'AbsTol', 1e-15);
   testCase.verifyEqual(diag.top_export_liquid_mwe, ...
      pair_liquid - retained_liquid, 'AbsTol', 1e-15);

   % A top removal exports mass and also translates the grid, but the two are
   % different quantities and must not be interconverted.
   testCase.verifyEqual(diag.top_deletion_count, 1);
   testCase.verifyGreaterThan(diag.top_export_solid_mwe, 0);
   testCase.verifyNotEqual( ...
      diag.top_export_solid_mwe, diag.top_deletion_height_m);
end

function test_runoff_credits_condensation_overflow(testCase)
   % Condensation the top cell could not store never entered the reservoir, so
   % it runs off directly. It was previously computed and then discarded, which
   % dropped that water out of the budget entirely.

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

   for profile = ["minimal", "standard", "diagnostic"]
      opts = icemodel.test.helpers.buildSyntheticOpts(workspace, ...
         'icemodel', 2016, output_profile=char(profile), solver=1);
      testCase.verifyTrue(ismember('df_rof', opts.vars1), ...
         sprintf('%s profile must carry df_rof', profile));
   end
   clear cleanup
end

function test_d_lyr_carries_total_merge_export_not_liquid_only(testCase)
   % d_lyr must accumulate the full water-equivalent mass a merge removes.
   % It previously recorded only the liquid difference, which is why the
   % derived dlayer series could never be reconciled against melt and runoff.

   [T, f_ice, f_liq, Sc, Sp, d_lyr] = mergeFixture([0.05; 0.6; 0.7]);
   dz = 0.04;
   [~, ~, ~, ~, ~, d_lyr_new, ~, diag] = ...
      icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, dz, 0.0, d_lyr, 0.1);

   % Scaling by the cell thickness must reproduce the ledger's independently
   % computed export, which is derived from column storage rather than from
   % the merged cell's fractions.
   exported_mwe = sum(d_lyr_new) * dz;
   testCase.verifyEqual(exported_mwe, ...
      diag.collapse_export_solid_mwe + diag.collapse_export_liquid_mwe, ...
      'AbsTol', 1e-15);

   % The liquid difference alone is strictly smaller, so the old definition
   % cannot satisfy the identity above.
   liquid_only = max(f_liq(1) + f_liq(2) - f_liq(1), 0) * dz;
   testCase.verifyLessThan(liquid_only, exported_mwe);

   % The three export views must nest rather than duplicate: the top-removal
   % channels are a subset of the all-merge collapse export, which in turn is
   % what d_lyr totals. A future reader must be able to tell these apart.
   testCase.verifyLessThanOrEqual(diag.top_export_solid_mwe, ...
      diag.collapse_export_solid_mwe + 1e-15);
   testCase.verifyLessThanOrEqual(diag.top_export_liquid_mwe, ...
      diag.collapse_export_liquid_mwe + 1e-15);
end

function test_interior_merge_exports_no_surface_mass(testCase)
   % Interior remeshing moves mass without lowering the surface, so it must
   % leave the surface-loss comparator untouched.

   [T, f_ice, f_liq, Sc, Sp, d_lyr] = ...
      mergeFixture([0.6; 0.05; 0.7; 0.8]);
   [~, ~, ~, ~, ~, ~, ~, diag] = icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, 0.04, 0.0, d_lyr, 0.1);

   testCase.verifyEqual(diag.interior_merge_count, 1);
   testCase.verifyEqual(diag.top_deletion_count, 0);
   testCase.verifyEqual(diag.top_export_solid_mwe, 0, 'AbsTol', 0);
   testCase.verifyEqual(diag.top_export_liquid_mwe, 0, 'AbsTol', 0);

   % The interior event still contributes to the total collapse export, so the
   % two channels are genuinely distinct rather than duplicates.
   testCase.verifyGreaterThan(diag.collapse_export_solid_mwe, 0);
end

function test_top_merge_counts_actual_grid_translation(testCase)
   % Removing the actual top cell contributes exactly one uniform-grid dz.

   [T, f_ice, f_liq, Sc, Sp, d_lyr] = mergeFixture([0.05; 0.6; 0.7]);
   [~, ~, ~, ~, ~, d_lyr_new, mask, diag] = ...
      icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, 0.04, 0.0, d_lyr, 0.1);

   % Geometry and the legacy non-geometric diagnostic remain distinct.
   testCase.verifyTrue(mask(1));
   testCase.verifyEqual(diag.top_deletion_count, 1);
   testCase.verifyEqual(diag.top_deletion_height_m, 0.04, 'AbsTol', 0);
   testCase.verifyEqual(diag.interior_merge_count, 0);
   testCase.verifyGreaterThan(sum(d_lyr_new), 0);
   verifyRemeshIdentity(testCase, diag);
end

function test_interior_merge_does_not_translate_top_grid(testCase)
   % Removing an interior cell is numerical remeshing, not surface lowering.

   [T, f_ice, f_liq, Sc, Sp, d_lyr] = ...
      mergeFixture([0.6; 0.05; 0.7; 0.8]);
   [~, ~, ~, ~, ~, ~, mask, diag] = ...
      icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, 0.04, 0.0, d_lyr, 0.1);

   % The actual-event ledger must not infer height from the eligibility mask.
   testCase.verifyTrue(mask(2));
   testCase.verifyEqual(diag.top_deletion_count, 0);
   testCase.verifyEqual(diag.top_deletion_height_m, 0, 'AbsTol', 0);
   testCase.verifyEqual(diag.interior_merge_count, 1);
   verifyRemeshIdentity(testCase, diag);
end

function test_multiple_top_merges_follow_index_drift_once_per_original_flag(testCase)
   % Adjacent flagged top cells shift to index one as each prior cell is removed.

   [T, f_ice, f_liq, Sc, Sp, d_lyr] = ...
      mergeFixture([0.05; 0.04; 0.6; 0.7]);
   [~, ~, ~, ~, ~, ~, mask, diag] = ...
      icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, 0.04, 0.0, d_lyr, 0.1);

   % Appended bottom clones are not reprocessed as new eligible cells.
   testCase.verifyEqual(sum(mask), 2);
   testCase.verifyEqual(diag.top_deletion_count, 2);
   testCase.verifyEqual(diag.top_deletion_height_m, 0.08, 'AbsTol', 0);
   testCase.verifyEqual(diag.interior_merge_count, 0);
   verifyRemeshIdentity(testCase, diag);
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
      % budget-aware hourly retime used by downstream comparison and reporting.
      testCase.verifyEqual(ledger.mass_budget_remesh_solid_mwe, ...
         ledger.mass_budget_cloned_bottom_solid_mwe ...
         - ledger.mass_budget_collapse_export_solid_mwe, 'AbsTol', 1e-15);
      testCase.verifyEqual(ledger.mass_budget_remesh_liquid_mwe, ...
         ledger.mass_budget_cloned_bottom_liquid_mwe ...
         - ledger.mass_budget_collapse_export_liquid_mwe, 'AbsTol', 1e-15);
      testCase.verifyEqual(ledger.mass_budget_remesh_enthalpy_j_m2, ...
         ledger.mass_budget_cloned_bottom_enthalpy_j_m2 ...
         - ledger.mass_budget_collapse_export_enthalpy_j_m2, 'AbsTol', 1e-8);

      discrepancy(k) = abs(target_height - grid_height(k));
      testCase.verifyLessThanOrEqual(discrepancy(k), dz(k));
   end

   testCase.verifyLessThan(discrepancy(2), discrepancy(1));
end

function test_multi_event_remesh_retains_non_cancelling_throughput(testCase)
   % Opposite-signed event exchanges must not disappear in the signed net.

   [T, f_ice, f_liq, Sc, Sp, d_lyr] = ...
      mergeFixture([0.05; 0.04; 0.8; 0.2]);
   [~, ~, ~, ~, ~, ~, ~, diag] = ...
      icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, 0.04, 0.0, d_lyr, 0.1);

   % The two top events exchange solid storage in opposite directions for this
   % column, so event throughput must strictly exceed the absolute signed net.
   testCase.verifyEqual(diag.top_deletion_count, 2);
   testCase.verifyGreaterThan( ...
      diag.solid_throughput_mwe, abs(diag.solid_mwe));
   testCase.verifyGreaterThanOrEqual( ...
      diag.cloned_bottom_solid_throughput_mwe, ...
      abs(diag.cloned_bottom_solid_mwe));
   testCase.verifyGreaterThanOrEqual( ...
      diag.collapse_export_solid_throughput_mwe, ...
      abs(diag.collapse_export_solid_mwe));
   verifyRemeshIdentity(testCase, diag);
end

function test_depleted_bottom_removal_matches_minimal_transition(testCase)
   % Requesting the ledger must never change the state transition, including
   % when the deepest cell is itself the merge-eligible layer.

   [T, f_ice, f_liq, Sc, Sp, d_lyr] = ...
      mergeFixture([0.6; 0.7; 0.05]);

   [T_min, f_ice_min_out, f_liq_min, Sc_min, Sp_min, d_lyr_min, mask] = ...
      icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, 0.04, 0.0, d_lyr, 0.1);
   [T_diag, f_ice_diag, f_liq_diag, Sc_diag, Sp_diag, d_lyr_diag, ~, diag] = ...
      icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, 0.04, 0.0, d_lyr, 0.1);

   % The deepest cell is removed and the surviving bottom state is cloned, so
   % the returned column keeps the legacy clone/delete result on both paths.
   testCase.verifyTrue(mask(end));
   testCase.verifyEqual(f_ice_min_out(end), f_ice(end), 'AbsTol', 0);
   testCase.verifyEqual(T_diag, T_min);
   testCase.verifyEqual(f_ice_diag, f_ice_min_out);
   testCase.verifyEqual(f_liq_diag, f_liq_min);
   testCase.verifyEqual(Sc_diag, Sc_min);
   testCase.verifyEqual(Sp_diag, Sp_min);
   testCase.verifyEqual(d_lyr_diag, d_lyr_min);

   % A deepest-cell removal is not a surface removal, so it records interior
   % merge activity and contributes no quantized grid translation.
   testCase.verifyEqual(diag.top_deletion_count, 0);
   testCase.verifyEqual(diag.interior_merge_count, 1);
   testCase.verifyEqual(diag.top_deletion_height_m, 0.0);
   verifyRemeshIdentity(testCase, diag);
end

function test_bottom_adjacent_merge_matches_minimal_transition(testCase)
   % A bottom-adjacent merge clones an already-depleted reservoir. Both paths
   % must still agree, and the ledger must close its domain-exchange identity.

   [T, f_ice, f_liq, Sc, Sp, d_lyr] = ...
      mergeFixture([0.6; 0.0; 0.15]);
   f_liq(:) = 0.0;
   [j1, j2] = icemodel.column.merge_layer_indices(2, f_ice);
   testCase.verifyEqual([j1, j2], [2, 3]);
   testCase.verifyGreaterThan(f_ice(end), 0.1);

   [T_min, f_ice_min_out, f_liq_min, Sc_min, Sp_min, d_lyr_min, mask] = ...
      icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, 0.04, 0.0, d_lyr, 0.1);
   [T_diag, f_ice_diag, f_liq_diag, Sc_diag, Sp_diag, d_lyr_diag, ~, diag] = ...
      icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, 0.04, 0.0, d_lyr, 0.1);

   testCase.verifyTrue(mask(2));
   testCase.verifyLessThanOrEqual(f_ice_min_out(end), 0.1);
   testCase.verifyEqual(T_diag, T_min);
   testCase.verifyEqual(f_ice_diag, f_ice_min_out);
   testCase.verifyEqual(f_liq_diag, f_liq_min);
   testCase.verifyEqual(Sc_diag, Sc_min);
   testCase.verifyEqual(Sp_diag, Sp_min);
   testCase.verifyEqual(d_lyr_diag, d_lyr_min);
   testCase.verifyEqual(diag.top_deletion_count, 0);
   testCase.verifyEqual(diag.interior_merge_count, 1);
   verifyRemeshIdentity(testCase, diag);
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
   [~, ~, ~, ~, ~, ~, mask] = icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, 0.04, d_pevp, d_lyr, f_ice_min);
   testCase.verifyEqual(mask, [true; false; false]);
end

function test_seven_output_merge_skips_scientific_storage_integration(testCase)
   % Existing solver callers must not execute the opt-in eighth-output ledger.

   [T, f_ice, f_liq, Sc, Sp, d_lyr] = mergeFixture([0.05; 0.6; 0.7]);
   profile clear
   cleanup = onCleanup(@() profile('off'));
   profile on
   [~, ~, ~, ~, ~, ~, mask] = icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, 0.04, 0.0, d_lyr, 0.1);
   profile off
   info = profile('info');

   % A real merge proves the early no-event return did not mask the branch.
   names = string({info.FunctionTable.FunctionName});
   testCase.verifyTrue(any(mask));
   testCase.verifyFalse(any(contains(names, 'budget_state')));
   clear cleanup
   profile clear
end

function test_pending_flags_match_legacy_multi_event_result(testCase)
   % A non-top pending mask must retain the historical appended bottom flag so
   % state and d_lyr evolution stay exactly unchanged, on BOTH output paths.
   % Requesting the ledger must never alter the transition it observes.

   [T, f_ice, f_liq, Sc, Sp, d_lyr] = ...
      mergeFixture([0.6; 0.05; 0.04]);
   [T_new, f_ice_new, f_liq_new, Sc_new, Sp_new, d_lyr_new, mask] = ...
      icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, 0.04, 0.0, d_lyr, 0.1);
   [T_diag, f_ice_diag, f_liq_diag, Sc_diag, Sp_diag, d_lyr_diag, ~, ~] = ...
      icemodel.column.merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, 0.04, 0.0, d_lyr, 0.1);
   [T_legacy, f_ice_legacy, f_liq_legacy, Sc_legacy, Sp_legacy, ...
      d_lyr_legacy, legacy_mask] = legacyMergeThinLayers( ...
      T, f_ice, f_liq, Sc, Sp, 0.04, 0.0, d_lyr, 0.1);

   % Compare every seven-output contract value, including the original mask.
   testCase.verifyEqual(mask, [false; true; true]);
   testCase.verifyEqual(mask, legacy_mask);
   testCase.verifyEqual(T_new, T_legacy);
   testCase.verifyEqual(f_ice_new, f_ice_legacy);
   testCase.verifyEqual(f_liq_new, f_liq_legacy);
   testCase.verifyEqual(Sc_new, Sc_legacy);
   testCase.verifyEqual(Sp_new, Sp_legacy);
   testCase.verifyEqual(d_lyr_new, d_lyr_legacy);

   % The diagnostic path shares one pending-flag rule with the legacy result.
   testCase.verifyEqual(T_diag, T_legacy);
   testCase.verifyEqual(f_ice_diag, f_ice_legacy);
   testCase.verifyEqual(f_liq_diag, f_liq_legacy);
   testCase.verifyEqual(Sc_diag, Sc_legacy);
   testCase.verifyEqual(Sp_diag, Sp_legacy);
   testCase.verifyEqual(d_lyr_diag, d_lyr_legacy);
end

function test_vapor_identity_partitions_wet_evaporation_and_sublimation(testCase)
   % Evaporation that exhausts mobile liquid must use Lv for liquid and Ls for
   % the remaining realized solid sublimation.

   [d_rof, d_sbl_err] = verifyVaporIdentity( ...
      testCase, 0.5, 0.05, -0.05, 0.1, 0.02);
   testCase.verifyEqual(d_rof, 0, 'AbsTol', 0);
   testCase.verifyEqual(d_sbl_err, 0, 'AbsTol', 0);
end

function test_vapor_identity_includes_condensation_overflow(testCase)
   % Wet condensation beyond cell capacity must be exposed without routing it.

   [d_rof, d_sbl_err] = verifyVaporIdentity( ...
      testCase, 0.99, 0.005, 0.01, 0.1, 0.02);
   testCase.verifyGreaterThan(d_rof, 0);
   testCase.verifyEqual(d_sbl_err, 0, 'AbsTol', 0);
end

function test_vapor_identity_retains_negative_unapplied_dry_exhaustion(testCase)
   % Sublimation demand beyond remaining dry ice must retain a negative signed
   % unapplied-energy remainder rather than disappear from the budget.

   [d_rof, d_sbl_err] = verifyVaporIdentity( ...
      testCase, 0.05, 0.0, -0.1, 0.1, 0.02);
   testCase.verifyEqual(d_rof, 0, 'AbsTol', 0);
   testCase.verifyLessThan(d_sbl_err, 0);
end

function test_vapor_identity_accepts_dry_sublimation(testCase)
   % Dry sublimation must close on the physical ro_liq/Lv to ro_ice/Ls basis.

   [d_rof, d_sbl_err] = verifyVaporIdentity( ...
      testCase, 0.5, 0.0, -0.02, 0.1, 0.02);
   testCase.verifyEqual(d_rof, 0, 'AbsTol', 0);
   testCase.verifyEqual(d_sbl_err, 0, 'AbsTol', 0);
end

function test_vapor_identity_accepts_dry_deposition(testCase)
   % Dry deposition must add solid mass on the same physical latent-heat basis.

   [d_rof, d_sbl_err] = verifyVaporIdentity( ...
      testCase, 0.5, 0.0, 0.02, 0.1, 0.02);
   testCase.verifyEqual(d_rof, 0, 'AbsTol', 0);
   testCase.verifyEqual(d_sbl_err, 0, 'AbsTol', 0);
end

function test_hourly_retime_uses_budget_aggregation_classes(testCase)
   % State endpoints select first/last while signed increments sum.

   time = datetime(2020, 1, 1, 'TimeZone', 'UTC') ...
      + minutes(15) * transpose(0:7);
   signed_phase = repmat([1; -1], 4, 1);
   TT = timetable(time, transpose(1:8), transpose(10:17), ...
      transpose(20:27), signed_phase, abs(signed_phase), ...
      'VariableNames', {'ordinary_mean', ...
      'mass_budget_solid_start_mwe', 'mass_budget_solid_end_mwe', ...
      'mass_budget_phase_solid_mwe', ...
      'mass_budget_phase_solid_throughput_mwe'});

   TT = icemodel.retimeHourlyFixedStep(TT);

   % The default mean remains intact for every unrelated output channel.
   testCase.verifyEqual(TT.ordinary_mean, [2.5; 6.5], 'AbsTol', 0);
   testCase.verifyEqual(TT.mass_budget_solid_start_mwe, [10; 14]);
   testCase.verifyEqual(TT.mass_budget_solid_end_mwe, [23; 27]);
   testCase.verifyEqual(TT.mass_budget_phase_solid_mwe, [0; 0]);
   testCase.verifyEqual( ...
      TT.mass_budget_phase_solid_throughput_mwe, [4; 4]);
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
      'mass_budget_phase_solid_throughput_mwe'});

   TT = icemodel.retimeHourlyFixedStep(TT);

   % Native hourly bins retain endpoint and non-cancelling sum semantics.
   testCase.verifyEqual(TT.Properties.RowTimes, ...
      datetime(2020, 1, 1, 'TimeZone', 'UTC'));
   testCase.verifyEqual(TT.ordinary_mean, 2);
   testCase.verifyEqual(TT.mass_budget_solid_start_mwe, 10);
   testCase.verifyEqual(TT.mass_budget_solid_end_mwe, 22);
   testCase.verifyEqual(TT.mass_budget_phase_solid_mwe, 2);
   testCase.verifyEqual(TT.mass_budget_phase_solid_throughput_mwe, 4);
end

function test_non_diagnostic_model_paths_omit_budget_channels(testCase)
   % Standard and minimal production runs must omit diagnostic bookkeeping
   % channels and preserve their externally visible profile schemas.

   profiles = {'standard', 'minimal'};
   for n = 1:numel(profiles)
      workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
         configure=true, nsteps=4, dt_seconds=3600);
      cleanup = onCleanup(@() ...
         icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));
      opts = icemodel.test.helpers.buildSyntheticOpts( ...
         workspace, 'icemodel', 2016, output_profile=profiles{n}, solver=1);

      % Exercise each non-diagnostic branch rather than only its option list.
      [ice1, ~, run_opts] = icemodel.test.helpers.runSmbModel(opts);
      testCase.verifyEqual(transpose(fieldnames(ice1)), run_opts.vars1, ...
         profiles{n});
      testCase.verifyFalse(any(ismember(fieldnames(ice1), ...
         icemodel.column.budget_output_fields())), profiles{n});
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

function verifyRemeshIdentity(testCase, diag)
   %VERIFYREMESHIDENTITY Check R = B - O for every accepted remesh reference.

   testCase.verifyEqual(diag.solid_mwe, ...
      diag.cloned_bottom_solid_mwe - diag.collapse_export_solid_mwe, ...
      'AbsTol', 1e-15);
   testCase.verifyEqual(diag.liquid_mwe, ...
      diag.cloned_bottom_liquid_mwe - diag.collapse_export_liquid_mwe, ...
      'AbsTol', 1e-15);
   testCase.verifyEqual(diag.enthalpy_j_m2, ...
      diag.cloned_bottom_enthalpy_j_m2 - diag.collapse_export_enthalpy_j_m2, ...
      'AbsTol', 1e-8);

   % Event-level absolute accumulators must bound every corresponding net.
   testCase.verifyGreaterThanOrEqual( ...
      diag.solid_throughput_mwe, abs(diag.solid_mwe));
   testCase.verifyGreaterThanOrEqual( ...
      diag.liquid_throughput_mwe, abs(diag.liquid_mwe));
   testCase.verifyGreaterThanOrEqual( ...
      diag.enthalpy_throughput_j_m2, abs(diag.enthalpy_j_m2));
   testCase.verifyGreaterThanOrEqual( ...
      diag.cloned_bottom_solid_throughput_mwe, ...
      abs(diag.cloned_bottom_solid_mwe));
   testCase.verifyGreaterThanOrEqual( ...
      diag.cloned_bottom_liquid_throughput_mwe, ...
      abs(diag.cloned_bottom_liquid_mwe));
   testCase.verifyGreaterThanOrEqual( ...
      diag.cloned_bottom_enthalpy_throughput_j_m2, ...
      abs(diag.cloned_bottom_enthalpy_j_m2));
   testCase.verifyGreaterThanOrEqual( ...
      diag.collapse_export_solid_throughput_mwe, ...
      abs(diag.collapse_export_solid_mwe));
   testCase.verifyGreaterThanOrEqual( ...
      diag.collapse_export_liquid_throughput_mwe, ...
      abs(diag.collapse_export_liquid_mwe));
   testCase.verifyGreaterThanOrEqual( ...
      diag.collapse_export_enthalpy_throughput_j_m2, ...
      abs(diag.collapse_export_enthalpy_j_m2));
end

function [d_rof, d_sbl_err] = verifyVaporIdentity( ...
      testCase, f_ice, f_liq, d_pevp, f_ice_min, f_res_por)
   %VERIFYVAPORIDENTITY Check the phase-aware accepted latent-energy identity.

   [Tf, Ls, Lv, ro_ice, ro_liq] = icemodel.physicalConstant( ...
      'Tf', 'Ls', 'Lv', 'ro_ice', 'ro_liq');
   dz = 0.04;
   [solid_p, liquid_p] = ...
      icemodel.column.integrate_column_budget(Tf - 2, f_ice, f_liq, dz);

   % Apply one accepted vapor increment through the production budget kernel.
   [~, f_ice_v, f_liq_v, ~, ~, d_rof, d_sbl_err] = ...
      icemodel.column.budget_surface_mass_balance( ...
      Tf - 2, f_ice, f_liq, f_liq, d_pevp, 0, 0, 0, ...
      f_res_por, f_ice_min);
   [solid_v, liquid_v] = ...
      icemodel.column.integrate_column_budget(Tf - 2, f_ice_v, f_liq_v, dz);

   % Potential = realized solid + realized liquid + overflow + unapplied, all
   % expressed with the physical solver densities used by d_pevp and sublimation.
   potential = ro_liq * Lv * d_pevp * dz;
   realized_and_overflow = ro_liq * ( ...
      Ls * (solid_v - solid_p) ...
      + Lv * (liquid_v - liquid_p) ...
      + Lv * d_rof * dz);
   unapplied = ro_ice * Ls * d_sbl_err * dz;
   testCase.verifyEqual(potential, realized_and_overflow + unapplied, ...
      'AbsTol', 1e-7);
end

function [T, f_ice, f_liq, Sc, Sp, d_lyr, merge_mask] = ...
      legacyMergeThinLayers( ...
      T, f_ice, f_liq, Sc, Sp, dz_therm, d_pevp, d_lyr, f_ice_min)
   %LEGACYMERGETHINLAYERS Reproduce the pre-ledger seven-output transition.

   [Ls, Lv, ro_ice, ro_liq] = ...
      icemodel.physicalConstant('Ls', 'Lv', 'ro_ice', 'ro_liq');
   merge_mask = f_ice <= f_ice_min | ...
      (f_ice + d_pevp * (Lv * ro_liq) / (Ls * ro_ice)) <= f_ice_min;
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
   % substeps, unlike d_pevp and d_sbl_err which are per-substep. The ledger
   % must therefore record the running total, not add it on every accepted
   % substep. Adding it once per substep inflated the overflow channel and the
   % vapor closure identity by the substep count.

   dz = 0.04;
   n_substeps = 4;
   overflow_fraction = 1e-4;

   % One overflow event on the first substep, then three quiet substeps that
   % still carry the accumulated total forward, which is the real failure mode.
   ledger = icemodel.column.initialize_budget_state();
   d_rof = 0.0;
   for n = 1:n_substeps
      if n == 1
         d_rof = d_rof + overflow_fraction;
      end
      state = struct('T_ice', 273.0, 'f_ice', 0.9, 'f_liq', 0.02);
      ledger = icemodel.column.accumulate_vapor_budget(ledger, ...
         0.9, 0.02, state.T_ice, state.f_ice, state.f_liq, dz, ...
         0.0, d_rof, 0.0);
   end

   returned = ledger.mass_budget_condensation_overflow_mwe;
   expected = overflow_fraction * dz;
   testCase.verifyEqual(returned, expected, AbsTol=1e-15)

   % The throughput magnitude must match for the same reason.
   testCase.verifyEqual( ...
      ledger.mass_budget_condensation_overflow_throughput_mwe, ...
      expected, AbsTol=1e-15)
end
