function tests = test_surface_phase_partitioning
%TEST_SURFACE_PHASE_PARTITIONING Verify dry/wet vapor partition behavior.
   tests = functiontests(localfunctions);
end

function test_dry_deposition_adds_ice(testCase)

   [ro_ice, ro_liq, Ls, Lv] = icemodel.physicalConstant( ...
      'ro_ice', 'ro_liq', 'Ls', 'Lv');
   f_ice = 0.95;
   f_liq = 0.0;
   d_con = 0.0;
   d_pevp = 1e-4;

   [f_ice_new, f_liq_new, d_con_new, x_vap] = ICESUBL( ...
      f_ice, f_liq, d_con, ro_ice, ro_liq, Ls, Lv, d_pevp, false, 0.1, 0.02);

   testCase.verifyGreaterThan(f_ice_new(1), f_ice(1));
   testCase.verifyEqual(f_liq_new(1), f_liq(1), 'AbsTol', 0);
   testCase.verifyEqual(d_con_new(1), 0, 'AbsTol', 0);
   testCase.verifyEqual(x_vap, 0, 'AbsTol', 0);
end

function test_wet_condensation_stays_liquid(testCase)

   [ro_ice, ro_liq, Ls, Lv] = icemodel.physicalConstant( ...
      'ro_ice', 'ro_liq', 'Ls', 'Lv');
   f_ice = 0.90;
   f_liq = 0.05;
   d_con = 0.0;
   d_pevp = 1e-4;

   [f_ice_new, f_liq_new, d_con_new, x_vap] = ICESUBL( ...
      f_ice, f_liq, d_con, ro_ice, ro_liq, Ls, Lv, d_pevp, true, 0.1, 0.02);

   testCase.verifyEqual(f_ice_new(1), f_ice(1), 'AbsTol', 0);
   testCase.verifyGreaterThan(f_liq_new(1), f_liq(1));
   testCase.verifyEqual(d_con_new(1), 0, 'AbsTol', 0);
   testCase.verifyEqual(x_vap, 0, 'AbsTol', 0);
end

function test_dry_sublimation_removes_ice(testCase)

   [ro_ice, ro_liq, Ls, Lv] = icemodel.physicalConstant( ...
      'ro_ice', 'ro_liq', 'Ls', 'Lv');
   f_ice = 0.95;
   f_liq = 0.0;
   d_con = 0.0;
   d_pevp = -1e-4;

   [f_ice_new, f_liq_new] = ICESUBL( ...
      f_ice, f_liq, d_con, ro_ice, ro_liq, Ls, Lv, d_pevp, false, 0.1, 0.02);

   testCase.verifyLessThan(f_ice_new(1), f_ice(1));
   testCase.verifyEqual(f_liq_new(1), f_liq(1), 'AbsTol', 0);
end

function test_vappress_honors_satflag(testCase)

   Tf = icemodel.physicalConstant('Tf');
   es_iceflag = VAPPRESS(Tf - 5, Tf, false);
   es_waterflag = VAPPRESS(Tf - 5, Tf, true);

   testCase.verifyNotEqual(es_iceflag, es_waterflag);
   testCase.verifyGreaterThan(es_waterflag, es_iceflag);
end
