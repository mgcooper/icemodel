function tests = test_colbeck_case_definition
   tests = functiontests(localfunctions);
end

function test_three_experiments(testCase)
   def = icemodel.verification.colbeck.caseDefinition();
   names = {def.experiments.name};
   verifyEqual(testCase, sort(names), {'exp1', 'exp2', 'exp3'});
   descriptions = {def.experiments.description};
   verifyTrue(testCase, ismember('ripe', descriptions));
   verifyTrue(testCase, ismember('refrozen', descriptions));
   verifyTrue(testCase, ismember('fresh', descriptions));
end

function test_grid_and_time_window(testCase)
   def = icemodel.verification.colbeck.caseDefinition();
   verifyEqual(testCase, def.grid.n_layers, 100);
   verifyEqual(testCase, def.grid.dz_m, 0.01);
   verifyEqual(testCase, def.grid.total_depth_m, ...
      def.grid.n_layers * def.grid.dz_m, 'AbsTol', 1e-12);

   verifyEqual(testCase, def.time.dt_s, 60);
   verifyEqual(testCase, def.time.n_steps, 600);
   verifyEqual(testCase, def.time.t_end_s, ...
      def.time.dt_s * def.time.n_steps);
   verifyEqual(testCase, def.time.rain_window_s, 10800);
   verifyLessThan(testCase, def.time.rain_window_s, def.time.t_end_s);
end

function test_surface_boundary_units(testCase)
   def = icemodel.verification.colbeck.caseDefinition();
   verifyEqual(testCase, def.boundary.q_top_kg_per_m2_per_s, 0.01, ...
      'AbsTol', 1e-12);

   ro_liq = icemodel.physicalConstant('ro_liq');
   expected = 0.01 / ro_liq;
   verifyEqual(testCase, def.boundary.q_top_m_per_s, expected, ...
      'AbsTol', 1e-12);
end

function test_initial_state_physical_bounds(testCase)
   def = icemodel.verification.colbeck.caseDefinition();
   Tf = icemodel.physicalConstant('Tf');
   for i = 1:numel(def.experiments)
      e = def.experiments(i);
      verifyGreaterThanOrEqual(testCase, e.f_ice, 0);
      verifyGreaterThanOrEqual(testCase, e.f_liq, 0);
      verifyLessThanOrEqual(testCase, e.f_ice + e.f_liq, 1);
      if strcmp(e.description, 'ripe')
         verifyEqual(testCase, e.T_K, 273.116, 'AbsTol', 1e-6);
         verifyGreaterThan(testCase, e.f_liq, 0);
      else
         verifyLessThanOrEqual(testCase, e.T_K, Tf);
         verifyEqual(testCase, e.f_liq, 0);
      end
      verifyGreaterThan(testCase, e.permeability_m2, 0);
   end
end

function test_params_present(testCase)
   def = icemodel.verification.colbeck.caseDefinition();
   verifyEqual(testCase, def.params.f_res_pore, 0.07, 'AbsTol', 1e-12);
   verifyFalse(testCase, isfield(def.params, 'm_kernel_icemodel'), ...
      'kernel exponent must come from parameterLookup, not caseDefinition');
   verifyFalse(testCase, isfield(def.params, 'mw_exp_summa'), ...
      'SUMMA Mualem exponent must not be carried in caseDefinition');
end

function test_grainsz_per_experiment(testCase)
   % Clark 2017 Table 1: 2 mm grain for exp1/exp2 (ripe + refrozen),
   % 0.2 mm for exp3 (fresh snow).
   def = icemodel.verification.colbeck.caseDefinition();
   for i = 1:numel(def.experiments)
      e = def.experiments(i);
      verifyTrue(testCase, isfield(e, 'grainsz_m'), ...
         sprintf('%s: grainsz_m field required', e.name));
      verifyGreaterThan(testCase, e.grainsz_m, 0);
   end
   verifyEqual(testCase, def.experiments(1).grainsz_m, 2.0e-3, 'AbsTol', 1e-12);
   verifyEqual(testCase, def.experiments(2).grainsz_m, 2.0e-3, 'AbsTol', 1e-12);
   verifyEqual(testCase, def.experiments(3).grainsz_m, 2.0e-4, 'AbsTol', 1e-12);
end

function test_provenance_strings(testCase)
   def = icemodel.verification.colbeck.caseDefinition();
   verifyTrue(testCase, isfield(def.provenance, 'clark_2021_eqs'));
   verifyTrue(testCase, isfield(def.provenance, 'colbeck_1976'));
   verifyTrue(testCase, contains(def.provenance.summa_initcond, ...
      "summa_zInitialCond_colbeck1976"));
end
