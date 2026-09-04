function tests = test_surface_phase_partitioning
   %TEST_SURFACE_PHASE_PARTITIONING Verify dry/wet vapor partition behavior.
   tests = functiontests(localfunctions);
end

function test_dry_deposition_adds_ice(testCase)
   % Dry deposition on a cold surface should go directly into the ice
   % fraction rather than creating a liquid film.

   f_ice = 0.95;
   f_liq = 0.0;
   d_rof = 0.0;
   d_pevp = 1e-4;
   f_ice_min = 0.1;
   f_res_por = 0.02;

   [f_ice_new, f_liq_new, d_rof, d_vap_liq, d_vap_ice] = ...
      icemodel.surface.apply_surface_vapor_exchange( ...
      f_ice, f_liq, d_rof, d_pevp, 1, f_ice_min, f_res_por);

   testCase.verifyGreaterThan(f_ice_new(1), f_ice(1));
   testCase.verifyEqual(f_liq_new(1), f_liq(1), 'AbsTol', 0);
   testCase.verifyEqual(d_rof, 0, 'AbsTol', 0);
   testCase.verifyEqual(d_vap_liq, 0, 'AbsTol', 0);
   testCase.verifyGreaterThan(d_vap_ice, 0);
end

function test_wet_condensation_stays_liquid(testCase)
   % Once the surface is already wet, condensation should remain in the
   % liquid branch instead of increasing the ice fraction.

   f_ice = 0.90;
   f_liq = 0.05;
   d_rof = 0.0;
   d_pevp = 1e-4;
   f_ice_min = 0.1;
   f_res_por = 0.02;

   [f_ice_new, f_liq_new, d_rof, d_vap_liq, d_vap_ice, d_vap] = ...
      icemodel.surface.apply_surface_vapor_exchange( ...
      f_ice, f_liq, d_rof, d_pevp, 1, f_ice_min, f_res_por);

   testCase.verifyEqual(f_ice_new(1), f_ice(1), 'AbsTol', 0);
   testCase.verifyGreaterThan(f_liq_new(1), f_liq(1));
   testCase.verifyEqual(d_rof, 0, 'AbsTol', 0);
   testCase.verifyGreaterThan(d_vap_liq, 0);
   testCase.verifyEqual(d_vap_ice, 0, 'AbsTol', 0);

   % Liquid applies at Lv on the liquid-water basis, so the realized
   % exchange equals the demand up to the partition round trip.
   testCase.verifyEqual(d_vap, d_pevp, 'RelTol', 1e-12);
end

function test_dry_sublimation_removes_ice(testCase)
   % Dry sublimation should remove ice mass without creating spurious
   % liquid water.

   f_ice = 0.95;
   f_liq = 0.0;
   d_rof = 0.0;
   d_pevp = -1e-4;
   f_ice_min = 0.1;
   f_res_por = 0.02;

   [f_ice_new, f_liq_new] = icemodel.surface.apply_surface_vapor_exchange( ...
      f_ice, f_liq, d_rof, d_pevp, 1, f_ice_min, f_res_por);

   testCase.verifyLessThan(f_ice_new(1), f_ice(1));
   testCase.verifyEqual(f_liq_new(1), f_liq(1), 'AbsTol', 0);
end

function test_thin_surface_sublimates_existing_ice_before_remesh(testCase)
   % A top cell below the remesh floor must still spend its existing ice.

   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   f_ice = 0.05;
   f_liq = 0;
   d_pevp = -0.10;

   [f_ice_new, f_liq_new, d_rof, d_vap_liq, d_vap_ice, d_vap] = ...
      icemodel.surface.apply_surface_vapor_exchange( ...
      f_ice, f_liq, 0, d_pevp, 1, 0.1, 0.02);

   available_lwe = f_ice * ro_ice / ro_liq;

   testCase.verifyEqual(f_ice_new, 0, 'AbsTol', 0);
   testCase.verifyEqual(f_liq_new, f_liq, 'AbsTol', 0);
   testCase.verifyEqual(d_rof, 0, 'AbsTol', 0);
   testCase.verifyEqual(d_vap_liq, 0, 'AbsTol', 0);

   % This one-cell column reaches its storage floor before it satisfies the
   % demand. The ice increment reports the ice-fraction change, and the
   % realized exchange stops at the ice the cell held.
   testCase.verifyEqual(d_vap_ice, -f_ice, 'RelTol', 1e-14);
   testCase.verifyEqual(d_vap, -available_lwe, 'RelTol', 1e-14);
end

function test_cascade_conserves_areal_energy_on_nonuniform_mesh(testCase)
   % A wet top cell supplies liquid and ice before the remaining sublimation
   % reaches the thicker dry cell. The thickness ratio must preserve the
   % total areal latent energy.

   [Ls, Lv, ro_ice, ro_liq] = ...
      icemodel.physicalConstant('Ls', 'Lv', 'ro_ice', 'ro_liq');
   f_ice_min = 0.1;
   f_ice = [0.12; 0.90];
   f_liq = [0.04; 0];
   dz = [0.02; 0.10];
   f_res_por = 0.02;

   % The demand removes mobile liquid and spendable ice from cell 1. Cell 2
   % can supply the remaining sublimation.
   d_pevp = -0.10;

   [f_ice_new, f_liq_new, d_rof, d_vap_liq, d_vap_ice, d_vap] = ...
      icemodel.surface.apply_surface_vapor_exchange( ...
      f_ice, f_liq, 0, d_pevp, dz, f_ice_min, f_res_por);

   % Cell 1 loses liquid and stops at the retained ice limit. The remaining
   % sublimation removes ice from cell 2.
   testCase.verifyEqual(f_ice_new(1), f_ice_min, 'RelTol', 1e-14);
   testCase.verifyEqual(d_vap_ice(1), f_ice_min - f_ice(1), ...
      'RelTol', 1e-14);
   testCase.verifyLessThan(f_liq_new(1), f_liq(1));
   testCase.verifyLessThan(d_vap_liq(1), 0);
   testCase.verifyLessThan(d_vap_ice(2), 0);
   testCase.verifyEqual(f_liq_new(2), f_liq(2), 'AbsTol', 0);
   testCase.verifyEqual(d_vap_liq(2), 0, 'AbsTol', 0);
   testCase.verifyEqual(d_rof, 0, 'AbsTol', 0);

   % The applied liquid and ice changes satisfy the full surface demand.
   testCase.verifyEqual(Lv * ro_liq * sum(d_vap_liq .* dz) ...
      + Ls * ro_ice * sum(d_vap_ice .* dz), ...
      Lv * ro_liq * d_pevp * dz(1), 'RelTol', 1e-14);

   % D_VAP reports the total applied mass as a fraction of cell 1.
   testCase.verifyEqual(d_vap, ...
      (sum(d_vap_liq .* dz) ...
      + sum(d_vap_ice .* dz) * ro_ice / ro_liq) / dz(1), ...
      'RelTol', 1e-14);
end

function test_full_top_cell_rejects_solid_deposition_without_cascade(testCase)
   % Solid deposition a nearly solid top cell cannot store must stay
   % unrealized at the surface: deposited ice belongs at the surface, so
   % the rejected remainder never cascades into the cell below.

   [Ls, Lv, ro_ice, ro_liq] = ...
      icemodel.physicalConstant('Ls', 'Lv', 'ro_ice', 'ro_liq');
   f_ice = [0.995; 0.90];
   f_liq = [0; 0];
   dz = [0.04; 0.04];

   % A dry deposition demand larger than cell 1's representable storage.
   d_pevp = 0.02;

   [f_ice_new, f_liq_new, d_rof, d_vap_liq, d_vap_ice, d_vap] = ...
      icemodel.surface.apply_surface_vapor_exchange( ...
      f_ice, f_liq, 0, d_pevp, dz, 0.1, 0.02);

   % Cell 1 gains ice; cell 2 and all liquid outputs stay unchanged.
   testCase.verifyGreaterThan(f_ice_new(1), f_ice(1));
   testCase.verifyEqual(f_ice_new(2), f_ice(2), 'AbsTol', 0);
   testCase.verifyEqual(f_liq_new, f_liq, 'AbsTol', 0);
   testCase.verifyEqual(d_vap_liq, zeros(2, 1), 'AbsTol', 0);
   testCase.verifyEqual(d_vap_ice(2), 0, 'AbsTol', 0);
   testCase.verifyEqual(d_rof, 0, 'AbsTol', 0);

   % d_vap reports only the realized part, which is smaller than the full
   % dry-demand realization at Ls.
   testCase.verifyEqual(d_vap, ...
      sum(d_vap_ice .* dz) * ro_ice / ro_liq / dz(1), 'RelTol', 1e-14);
   testCase.verifyLessThan(d_vap, d_pevp * Lv / Ls);
end

function test_vappress_honors_satflag(testCase)
   % The phase flag should switch between ice and liquid saturation vapor
   % pressure curves at the same temperature.

   Tf = icemodel.physicalConstant('Tf');
   es_iceflag = icemodel.vapor.saturation_vapor_pressure(Tf - 5, false);
   es_waterflag = icemodel.vapor.saturation_vapor_pressure(Tf - 5, true);

   testCase.verifyNotEqual(es_iceflag, es_waterflag);
   testCase.verifyGreaterThan(es_waterflag, es_iceflag);
end

function test_scalar_demand_reaches_the_top_cell_alone(testCase)
   % The surface function passes one demand for the whole column. It must
   % land in the top cell and leave every interior cell untouched, which is
   % the same rule the merge look-ahead follows.

   f_ice = [0.90; 0.90; 0.90];
   f_liq = [0.05; 0.05; 0.05];

   [f_ice_new, f_liq_new, d_rof, d_vap_liq, d_vap_ice] = ...
      icemodel.surface.apply_surface_vapor_exchange( ...
      f_ice, f_liq, 0, 1e-4, ones(3, 1), 0.1, 0.02);

   testCase.verifyGreaterThan(f_liq_new(1), f_liq(1));
   testCase.verifyEqual(f_liq_new(2:3), f_liq(2:3), 'AbsTol', 0);
   testCase.verifyEqual(f_ice_new, f_ice, 'AbsTol', 0);
   testCase.verifyEqual(d_rof, 0, 'AbsTol', 0);

   % The applied-increment outputs contain one entry per cell. Wet deposition
   % enters only the top cell's liquid output.
   testCase.verifySize(d_vap_liq, [3, 1]);
   testCase.verifyEqual(d_vap_liq(1), f_liq_new(1) - f_liq(1), 'AbsTol', 0);
   testCase.verifyEqual(d_vap_liq(2:3), zeros(2, 1), 'AbsTol', 0);
   testCase.verifyEqual(d_vap_ice, zeros(3, 1), 'AbsTol', 0);
end

function test_wet_evaporation_crosses_from_liquid_to_ice_exactly(testCase)
   % A demand larger than the mobile liquid must spend that liquid at Lv and
   % only its remaining energy on ice at Ls. The surface function must preserve
   % the mixed-latent state and energy accounting.

   [Ls, Lv, ro_ice, ro_liq] = ...
      icemodel.physicalConstant('Ls', 'Lv', 'ro_ice', 'ro_liq');
   f_ice = 0.80;
   f_liq = 0.05;
   d_rof = 0.03;
   f_ice_min = 0.1;
   f_res_por = 0.02;
   [~, f_res] = icemodel.column.vapor_exchange_is_wet( ...
      f_ice, f_liq, f_res_por);
   d_pevp = -2 * (f_liq - f_res);

   [d_vap_liq, d_vap_ice, returned_f_res] = ...
      icemodel.surface.potential_surface_vapor_exchange( ...
      d_pevp, f_ice, f_liq, f_res_por);

   testCase.verifyEqual(returned_f_res, f_res, 'AbsTol', 0);
   testCase.verifyEqual(d_vap_liq, -(f_liq - f_res), 'RelTol', 1e-14);
   testCase.verifyLessThan(d_vap_ice, 0);
   testCase.verifyEqual(Lv * d_pevp, ...
      Lv * d_vap_liq + Ls * d_vap_ice, 'RelTol', 1e-14);

   [f_ice_actual, f_liq_actual, d_rof_actual, d_vap_liq_actual, ...
      d_vap_ice_actual, d_vap_actual] = ...
      icemodel.surface.apply_surface_vapor_exchange( ...
      f_ice, f_liq, d_rof, d_pevp, 1, f_ice_min, f_res_por);

   % Replay the legacy wet-crossing branch without calling either new
   % transfer helper. Liquid reaches its residual floor before the remaining
   % Lv-basis demand is converted to an ice-fraction change at Ls.
   d_aevp_oracle = -(f_liq - f_res);
   d_pevp_ice_oracle = d_pevp - d_aevp_oracle;
   d_psbl_oracle = d_pevp_ice_oracle ...
      * (Lv * ro_liq) / (Ls * ro_ice);
   f_liq_oracle = f_liq + d_aevp_oracle;
   f_ice_oracle = f_ice + d_psbl_oracle;
   d_rof_oracle = d_rof;
   d_vap_oracle = (f_liq_oracle - f_liq) ...
      + (f_ice_oracle - f_ice) * ro_ice / ro_liq;

   testCase.verifyEqual(f_ice_actual, f_ice_oracle, 'RelTol', 1e-14);
   testCase.verifyEqual(f_liq_actual, f_liq_oracle, 'AbsTol', 0);
   testCase.verifyEqual(d_rof_actual, d_rof_oracle, 'AbsTol', 0);

   % The per-phase outputs report each phase on its own fraction basis.
   testCase.verifyEqual(d_vap_liq_actual, d_aevp_oracle, 'RelTol', 1e-14);
   testCase.verifyEqual(d_vap_ice_actual, d_psbl_oracle, 'RelTol', 1e-14);
   testCase.verifyEqual(d_vap_actual, d_vap_oracle, 'RelTol', 1e-14);
end

function test_shared_transfer_reports_each_phase_on_the_lwe_basis(testCase)
   % APPLY_VAPOR_TRANSFER applies both phase increments and returns each
   % limited remainder.

   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   f_ice_min = 0.1;
   f_ice = [0.99; f_ice_min];
   f_liq = [0.0; 0.02];
   f_res = [0.0; 0.02];
   d_vap_liq = [0.02; -0.01];
   d_vap_ice = [0.02; -0.01];

   [f_ice_new, f_liq_new, d_liq_unapplied, d_ice_unapplied] = ...
      icemodel.column.apply_vapor_transfer(f_ice, f_liq, ...
      d_vap_liq, d_vap_ice, f_ice_min, f_res);

   liquid_applied = f_liq_new - f_liq;
   ice_applied = (f_ice_new - f_ice) * ro_ice / ro_liq;
   testCase.verifyEqual(liquid_applied + d_liq_unapplied, ...
      d_vap_liq, 'AbsTol', 1e-16);
   testCase.verifyEqual(ice_applied + d_ice_unapplied, ...
      d_vap_ice, 'AbsTol', 1e-16);
   testCase.verifyTrue(any(d_liq_unapplied ~= 0));
   testCase.verifyTrue(any(d_ice_unapplied ~= 0));
end
