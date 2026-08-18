function tests = test_coupled_vapor_model_run
   %TEST_COUPLED_VAPOR_MODEL_RUN Verify the coupled vapor mode in a real run.
   %
   % The kernel tests in test_coupled_vapor_transport check the transport and
   % the application in isolation. These run the model itself, which is what
   % the DesignSpec's acceptance policy asks for: the coupled mode must not
   % regress the default, and its ledger must close.
   %
   % See also: test_coupled_vapor_transport,
   %  icemodel.column.couple_vapor_transport
   tests = functiontests(localfunctions);
end

function test_the_flag_is_off_by_default(testCase)
   % Nothing may enable the coupled mode. The DesignSpec keeps production
   % wiring as a separate deliverable (icemodel-bhk.5), and the regression
   % baselines depend on that.

   opts = icemodel.setopts("icemodel", "kanm", 2016, "kanm");
   testCase.verifyFalse(opts.use_coupled_vapor);
end

function test_flag_off_reproduces_the_default_run(testCase)
   % Setting the flag false explicitly must give exactly the default result.
   % This is the guard on the whole opt-in design: every branch the coupled
   % mode adds has to be invisible when it is off.

   [base_opts, cleanup] = syntheticRunOpts(testCase);
   [ice1_default, ~] = icemodel.test.helpers.runSmbModel(base_opts);
   [ice1_off, ~] = icemodel.test.helpers.runSmbModel( ...
      icemodel.resetopts(base_opts, 'use_coupled_vapor', false));

   fields = fieldnames(ice1_default);
   for k = 1:numel(fields)
      value = ice1_default.(fields{k});
      if isnumeric(value)
         testCase.verifyEqual(ice1_off.(fields{k}), value, ...
            sprintf('field %s moved with the flag off', fields{k}));
      end
   end
   clear cleanup
end

function test_an_opts_struct_without_the_flag_still_runs(testCase)
   % An opts struct built before this flag existed carries no
   % use_coupled_vapor field, and icemodel.configureRun does not add one.
   % Reading it with getopts would throw before the default path could run,
   % so the driver reads it defensively. A saved struct must still run.

   [base_opts, cleanup] = syntheticRunOpts(testCase);
   testCase.assertTrue(isfield(base_opts, 'use_coupled_vapor'));

   older_opts = rmfield(base_opts, 'use_coupled_vapor');
   testCase.assertFalse(isfield(older_opts, 'use_coupled_vapor'));

   [ice1_older, ~] = icemodel.test.helpers.runSmbModel(older_opts);
   [ice1_default, ~] = icemodel.test.helpers.runSmbModel(base_opts);

   % The struct without the field takes the default path, so every numeric
   % output must match the run that carries the field set false.
   fields = fieldnames(ice1_default);
   for k = 1:numel(fields)
      value = ice1_default.(fields{k});
      if isnumeric(value)
         testCase.verifyEqual(ice1_older.(fields{k}), value, ...
            sprintf('field %s moved without the flag', fields{k}));
      end
   end
   clear cleanup
end

function test_the_redistribution_channels_are_zero_when_off(testCase)
   % The redistribution channels must read zero on a default run;
   % matching another default run is not enough. The per-phase storage
   % closures consume them, so a nonzero value with the flag off would
   % move those closures for every existing diagnostic run.

   [base_opts, cleanup] = syntheticRunOpts(testCase);
   [ice1, ~] = icemodel.test.helpers.runSmbModel(base_opts);

   channels = { ...
      'mass_budget_vapor_redistribution_solid_mwe', ...
      'mass_budget_vapor_redistribution_liquid_mwe', ...
      'mass_budget_vapor_redistribution_unapplied_j_m2', ...
      'mass_budget_vapor_redistribution_solid_gross_mwe', ...
      'mass_budget_vapor_redistribution_liquid_gross_mwe', ...
      'mass_budget_vapor_redistribution_unapplied_gross_j_m2'};
   for n = 1:numel(channels)
      testCase.verifyEqual(ice1.(channels{n}), ...
         zeros(size(ice1.(channels{n}))), channels{n});
   end
   clear cleanup
end

function test_coupled_run_closes_the_vapor_identity(testCase)
   % The vapor closure identity must hold in a coupled run, not only in the
   % default one. The potential input equals the realized storage change plus
   % the overflow plus whatever the column could not apply:
   %
   %   potential = ro_liq * (Ls * solid + Lv * liquid + Lv * overflow)
   %               + unapplied
   %
   % No redistribution term appears: the transport runs after the vapor
   % budget takes its storage baseline, so the surface channels carry
   % surface exchange alone and the transport keeps its own accounting.

   [base_opts, cleanup] = syntheticRunOpts(testCase);
   [Ls, Lv, ro_liq] = icemodel.physicalConstant('Ls', 'Lv', 'ro_liq');

   [ice1, ~] = icemodel.test.helpers.runSmbModel( ...
      icemodel.resetopts(base_opts, 'use_coupled_vapor', true));

   % The transport runs after the vapor budget takes its storage baseline,
   % so the surface closure identity carries no redistribution term.
   accounted = ro_liq * ( ...
      Ls * ice1.mass_budget_vapor_solid_mwe ...
      + Lv * ice1.mass_budget_vapor_liquid_mwe ...
      + Lv * ice1.mass_budget_condensation_overflow_mwe) ...
      + ice1.mass_budget_unapplied_vapor_j_m2;
   potential = ice1.mass_budget_vapor_potential_j_m2;

   % The tolerance is subtraction roundoff on column-integrated checkpoints,
   % the same source the default-mode identity test allows for. The coupled
   % run integrates over every cell rather than the top one, so it collects
   % more of it.
   testCase.verifyEqual(potential, accounted, 'AbsTol', 1e-4);

   % A relative check as well, because the absolute tolerance above means
   % little without the scale it sits on.
   scale = max(abs(potential));
   testCase.verifyGreaterThan(scale, 0);
   testCase.verifyLessThan(max(abs(potential - accounted)) / scale, 1e-8);
   clear cleanup
end

function test_coupled_run_closes_the_per_phase_storage(testCase)
   % The per-phase storage closures must hold in a coupled run through the
   % production driver: the endpoint deltas equal the phase, surface
   % vapor, remesh, and redistribution increments together. This exercises
   % the driver's accepted-substep order end to end. Running the transport
   % before the vapor budget would count cross-phase transport in both the
   % vapor and the redistribution channels and break this closure by the
   % double-counted amount, and the endpoint gross bound would lose the
   % transport's contribution.

   [base_opts, cleanup] = syntheticRunOpts(testCase);
   [ice1, ~] = icemodel.test.helpers.runSmbModel( ...
      icemodel.resetopts(base_opts, 'use_coupled_vapor', true));

   solid_delta = ice1.mass_budget_solid_end_mwe ...
      - ice1.mass_budget_solid_start_mwe;
   liquid_delta = ice1.mass_budget_liquid_end_mwe ...
      - ice1.mass_budget_liquid_start_mwe;

   testCase.verifyEqual(solid_delta, ...
      ice1.mass_budget_phase_solid_mwe ...
      + ice1.mass_budget_vapor_solid_mwe ...
      + ice1.mass_budget_remesh_solid_mwe ...
      + ice1.mass_budget_vapor_redistribution_solid_mwe, 'AbsTol', 1e-10);
   testCase.verifyEqual(liquid_delta, ...
      ice1.mass_budget_phase_liquid_mwe ...
      + ice1.mass_budget_vapor_liquid_mwe ...
      + ice1.mass_budget_remesh_liquid_mwe ...
      + ice1.mass_budget_vapor_redistribution_liquid_mwe, 'AbsTol', 1e-10);

   % The endpoint gross bounds its signed delta with the transport
   % included in the substep storage-change context.
   testCase.verifyGreaterThanOrEqual( ...
      ice1.mass_budget_solid_storage_gross_mwe + 1e-12, abs(solid_delta));
   testCase.verifyGreaterThanOrEqual( ...
      ice1.mass_budget_liquid_storage_gross_mwe + 1e-12, ...
      abs(liquid_delta));
   clear cleanup
end

function test_coupled_run_changes_the_result(testCase)
   % A mode that changed nothing would pass every check above while doing
   % nothing. The coupled run must differ from the default one.

   [base_opts, cleanup] = syntheticRunOpts(testCase);
   [ice1_off, ~] = icemodel.test.helpers.runSmbModel(base_opts);
   [ice1_on, ~] = icemodel.test.helpers.runSmbModel( ...
      icemodel.resetopts(base_opts, 'use_coupled_vapor', true));

   testCase.verifyNotEqual(ice1_on.mass_budget_vapor_solid_mwe, ...
      ice1_off.mass_budget_vapor_solid_mwe);
   clear cleanup
end

function [base_opts, cleanup] = syntheticRunOpts(testCase)
   %SYNTHETICRUNOPTS Build one short diagnostic run on a synthetic workspace.

   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=24, dt_seconds=3600);
   cleanup = onCleanup(@() ...
      icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace));
   base_opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'icemodel', 2016, output_profile='diagnostic', solver=1);
   testCase.assertFalse(base_opts.use_coupled_vapor);
end
