function tests = test_vapor_face_discretization
   %TEST_VAPOR_FACE_DISCRETIZATION Verify the shared vapor face rule.
   %
   % vapor_transport_terms owns the face interpolation, vapor mass
   % flux, positive donor-tangent energy term, deferred-correction flux, and
   % face donor latent heat. Its outputs must be conjugate. Energy across a
   % face must equal the mass flux times the donor phase's latent heat. A
   % column closure test can still balance a wrong face because adjacent
   % face errors can cancel.
   %
   % See also: icemodel.column.vapor_transport_terms
   tests = functiontests(localfunctions);
end

function test_face_energy_equals_latent_heat_times_face_mass(testCase)
   % The identity the whole conjugate discretization exists to hold:
   %
   %   k_vap * (T_north - T_south) / delz + q_deferred == L_face * U_vap
   %
   % at every interior face. The matrix part alone cannot reproduce the
   % flux, because it carries the donor tangent rather than the secant; the
   % deferred term restores the difference, so the sum must return the flux.

   [T, f_ice, f_liq, delz, fn, f_res_por] = faceFixture(6);
   [ro_vap, De, dro_vapdT] = nodeQuantities(T, f_liq);
   k_eff = icemodel.column.bulk_thermal_conductivity(T, f_ice, f_liq, 0);
   [~, k_vap_faces, q_vap_deferred, U_vap_faces, L_vap_faces] = ...
      icemodel.column.vapor_transport_terms( ...
      T, f_ice, f_liq, k_eff, ro_vap, dro_vapdT, De, delz, fn, f_res_por);

   % Rebuild the latent heat and the temperature difference the conductance
   % used, so the check reads the same faces the functions did.
   [L_faces, d_T] = faceLatentHeat(T, f_ice, f_liq, ro_vap, f_res_por);

   % The function's own returned face donor latent heat must match this
   % independent reconstruction exactly: both apply the same donor rule to
   % the same padded node arrays.
   testCase.verifyEqual(L_vap_faces, L_faces, 'AbsTol', 0);

   % Interior faces alone. Both boundary faces are set to zero by design and
   % carry no identity to check.
   interior = 2:numel(T);
   returned = k_vap_faces(interior) .* d_T(interior) ./ delz(interior) ...
      + q_vap_deferred(interior);
   expected = L_faces(interior) .* U_vap_faces(interior);

   scale = max(abs(expected));
   testCase.assertGreaterThan(scale, 0);
   testCase.verifyEqual(returned, expected, 'RelTol', 1e-12, ...
      'AbsTol', 1e-14 * scale);
end

function test_the_identity_holds_across_a_wet_dry_boundary(testCase)
   % The face where the two latent heats meet is the one a separate energy
   % rule would get wrong. The donor cell decides the latent heat, so a face
   % between a wet and a dry cell must carry the donor's, and the identity
   % must still close there.

   [T, f_ice, ~, delz, fn, f_res_por] = faceFixture(4);

   % Wet above, dry below, so face 3 spans the phase change.
   f_liq = [0.05; 0.05; 0; 0];
   [ro_vap, De, dro_vapdT] = nodeQuantities(T, f_liq);
   k_eff = icemodel.column.bulk_thermal_conductivity(T, f_ice, f_liq, 0);
   [~, k_vap_faces, q_vap_deferred, U_vap_faces, L_vap_faces] = ...
      icemodel.column.vapor_transport_terms( ...
      T, f_ice, f_liq, k_eff, ro_vap, dro_vapdT, De, delz, fn, f_res_por);
   [L_faces, d_T] = faceLatentHeat(T, f_ice, f_liq, ro_vap, f_res_por);

   % The spanning face takes one of the two latent heats, not a blend, and
   % the function's own output must match the independent reconstruction.
   [Ls, Lv] = icemodel.physicalConstant('Ls', 'Lv');
   testCase.verifyTrue(L_faces(3) == Ls || L_faces(3) == Lv);
   testCase.verifyEqual(L_vap_faces(3), L_faces(3), 'AbsTol', 0);

   returned = k_vap_faces(3) * d_T(3) / delz(3) + q_vap_deferred(3);
   expected = L_faces(3) * U_vap_faces(3);
   testCase.verifyEqual(returned, expected, 'RelTol', 1e-12);
end

function test_both_boundary_faces_carry_no_vapor_energy(testCase)
   % The bottom is closed and the surface exchange enters through the SEB's
   % Qe. A vapor term at either boundary would count the surface exchange
   % twice, which is what the Neumann boundary decision exists to prevent.
   % Both the matrix part and the deferred flux must be zero there.

   [T, f_ice, f_liq, delz, fn, f_res_por] = faceFixture(5);
   [ro_vap, De, dro_vapdT] = nodeQuantities(T, f_liq);
   k_eff = icemodel.column.bulk_thermal_conductivity(T, f_ice, f_liq, 0);
   [~, k_vap_faces, q_vap_deferred, U_vap_faces] = ...
      icemodel.column.vapor_transport_terms( ...
      T, f_ice, f_liq, k_eff, ro_vap, dro_vapdT, De, delz, fn, f_res_por);

   testCase.verifyEqual(k_vap_faces(1), 0);
   testCase.verifyEqual(k_vap_faces(numel(T) + 1), 0);
   testCase.verifyEqual(q_vap_deferred(1), 0);
   testCase.verifyEqual(q_vap_deferred(numel(T) + 1), 0);
   testCase.verifyEqual(U_vap_faces(1), 0);
   testCase.verifyEqual(U_vap_faces(numel(T) + 1), 0);
end

function test_an_isothermal_column_moves_no_vapor(testCase)
   % Equal temperatures and one phase give equal saturation densities, so
   % every face flux is zero. This is the reproduction case the acceptance
   % policy names: an isothermal column must leave the surface path acting
   % alone. The deferred term must also vanish, so the energy side moves
   % nothing either.

   [~, f_ice, f_liq, delz, fn, f_res_por] = faceFixture(5);
   Tf = icemodel.physicalConstant('Tf');
   T = (Tf - 5) * ones(5, 1);
   [ro_vap, De, dro_vapdT] = nodeQuantities(T, f_liq);
   k_eff = icemodel.column.bulk_thermal_conductivity(T, f_ice, f_liq, 0);
   [~, ~, q_vap_deferred, returned] = ...
      icemodel.column.vapor_transport_terms( ...
      T, f_ice, f_liq, k_eff, ro_vap, dro_vapdT, De, delz, fn, f_res_por);
   expected = zeros(6, 1);
   testCase.verifyEqual(returned, expected, 'AbsTol', 0);
   testCase.verifyEqual(q_vap_deferred, expected, 'AbsTol', 0);
end

function test_the_matrix_part_is_the_positive_donor_tangent(testCase)
   % The matrix part carries the donor cell's tangent slope. Saturation
   % vapor density rises with temperature, so the tangent is positive at
   % every face, in every state; that positivity is what keeps the
   % assembled system diagonally dominant. A face that spans no temperature
   % difference still gets a finite, nonzero matrix coefficient from the
   % tangent.

   [~, f_ice, f_liq, delz, fn, f_res_por] = faceFixture(4);
   Tf = icemodel.physicalConstant('Tf');

   % Faces 2 and 4 span no temperature difference; face 3 does.
   T = [Tf - 6; Tf - 6; Tf - 2; Tf - 2];
   [ro_vap, De, dro_vapdT] = nodeQuantities(T, f_liq);
   k_eff = icemodel.column.bulk_thermal_conductivity(T, f_ice, f_liq, 0);
   [~, k_vap_faces] = icemodel.column.vapor_transport_terms( ...
      T, f_ice, f_liq, k_eff, ro_vap, dro_vapdT, De, delz, fn, f_res_por);

   testCase.verifyTrue(all(isfinite(k_vap_faces)));
   testCase.verifyTrue(all(k_vap_faces >= 0));
   testCase.verifyGreaterThan(abs(k_vap_faces(2)), 0);
end

function test_the_face_diffusivity_is_the_harmonic_mean(testCase)
   % Patankar (1980) Eq. 4.9. The harmonic mean makes a face between a
   % conducting cell and an insulating one behave like the insulating one,
   % which is the property an arithmetic mean loses.

   % Two nodes span one open interior face. Give vapor density a unit jump so
   % the returned mass flux exposes the face diffusivity directly.
   T = [267; 266];
   f_ice = 0.5 * ones(2, 1);
   f_liq = zeros(2, 1);
   ro_vap = [1; 2];
   dro_vapdT = ones(2, 1);
   De = [1e-9; 1e-5];
   k_eff = ones(2, 1);
   delz = ones(3, 1);
   fn = [0.5; 0.5; 0.5];
   f_res_por = 0.02;

   [~, ~, ~, U_vap_faces] = ...
      icemodel.column.vapor_transport_terms( ...
      T, f_ice, f_liq, k_eff, ro_vap, dro_vapdT, De, delz, fn, f_res_por);
   De_face = 1 / (0.5 / De(1) + 0.5 / De(2));
   returned = -U_vap_faces(2) / (ro_vap(2) - ro_vap(1));
   testCase.verifyEqual(returned, De_face, 'RelTol', 1e-15);

   % The insulating node dominates: the face sits nearer the small value than
   % the arithmetic mean would put it.
   testCase.verifyLessThan(returned, mean(De));

   % Both boundary mass fluxes are closed even though their conductivity is
   % defined for the thermal equation.
   testCase.verifyEqual(U_vap_faces([1, 3]), zeros(2, 1));
end

function test_interface_conductivity_adds_the_harmonic_bulk_term(testCase)
   % The helper must return the combined face conductivity expected by the
   % assembler: harmonic non-vapor conduction plus the vapor matrix term.

   [T, f_ice, f_liq, delz, fn, f_res_por] = faceFixture(5);
   [ro_vap, De, dro_vapdT] = nodeQuantities(T, f_liq);
   k_eff = linspace(0.5, 2.5, numel(T))';

   [k_eff_faces, k_vap_faces] = ...
      icemodel.column.vapor_transport_terms( ...
      T, f_ice, f_liq, k_eff, ro_vap, dro_vapdT, De, delz, fn, f_res_por);

   % Rebuild only the ordinary conduction term independently.
   expected = 1.0 ./ ((1.0 - fn) ./ [k_eff(1); k_eff] ...
      + fn ./ [k_eff; k_eff(end)]);

   testCase.verifyEqual(k_eff_faces - k_vap_faces, expected, ...
      'RelTol', 1e-15);
end

function test_the_deferred_term_carries_an_isothermal_wet_dry_face(testCase)
   % Two neighbours at equal temperature but different phase hold
   % different saturation vapor densities, because the ice and water
   % curves differ. Mass crosses the face while d_T is zero, so no
   % conductance times a temperature difference can carry the energy. The
   % deferred term must carry all of it: the face energy must equal
   % L_face * U_vap with the matrix part contributing nothing.

   [~, f_ice, ~, delz, fn, f_res_por] = faceFixture(4);
   Tf = icemodel.physicalConstant('Tf');

   % Isothermal, wet over dry, so face 3 spans the phase change at equal T.
   T = (Tf - 5) * ones(4, 1);
   f_liq = [0.05; 0.05; 0; 0];
   [ro_vap, De, dro_vapdT] = nodeQuantities(T, f_liq);
   k_eff = icemodel.column.bulk_thermal_conductivity(T, f_ice, f_liq, 0);
   [~, k_vap_faces, q_vap_deferred, U_vap_faces, L_vap_faces] = ...
      icemodel.column.vapor_transport_terms( ...
      T, f_ice, f_liq, k_eff, ro_vap, dro_vapdT, De, delz, fn, f_res_por);
   [L_faces, d_T] = faceLatentHeat(T, f_ice, f_liq, ro_vap, f_res_por);

   % The function's own returned face donor latent heat must match this
   % independent reconstruction exactly.
   testCase.verifyEqual(L_vap_faces(3), L_faces(3), 'AbsTol', 0);

   % Mass crosses the face at zero temperature difference.
   testCase.assertGreaterThan(abs(U_vap_faces(3)), 0);
   testCase.assertEqual(d_T(3), 0);

   % The matrix part moves nothing there; the deferred term moves L * U.
   matrix_part = k_vap_faces(3) * d_T(3) / delz(3);
   required = L_faces(3) * U_vap_faces(3);
   testCase.verifyEqual(matrix_part, 0);
   testCase.verifyGreaterThan(abs(q_vap_deferred(3)), 0);
   testCase.verifyEqual(matrix_part + q_vap_deferred(3), required, ...
      'RelTol', 1e-12);
end

function test_the_assembled_coefficients_keep_diagonal_dominance(testCase)
   % The dominance property lives in the assembled coefficients, so this
   % asserts on them, not on the face terms alone. In the worst region the
   % verification sweep found, low-density snow with a wet half over a dry
   % half, the secant slope is negative and its magnitude reaches 21.85
   % times the face conduction conductance. The split keeps that sign in
   % the deferred source: every assembled off-diagonal must stay
   % nonnegative and every row diagonally dominant, while the converged
   % face flux still carries the Bergeron-Findeisen direction against the
   % temperature gradient.

   Tf = icemodel.physicalConstant('Tf');
   JJ = 20;
   f_res_por = 0.02;
   [dz, delz, ~, ~, fn] = icemodel.column.control_volume_mesh( ...
      JJ * 0.04, 0.04);
   dz = dz(1:JJ);
   delz = delz(1:JJ + 1);
   fn = fn(1:JJ + 1);

   % The worst region the V4 sweep found: low-density snow, shallow
   % gradient, a wet upper half over a dry lower half.
   T = (Tf - 10) + linspace(0, 0.5, JJ)';
   f_ice = 0.30 * ones(JJ, 1);
   f_liq = zeros(JJ, 1);
   f_liq(1:JJ / 2) = 0.05;

   [ro_vap, dro_vapdT] = ...
      icemodel.vapor.saturation_vapor_density(T, f_liq);
   [~, De] = icemodel.vapor.vapor_thermal_conductivity( ...
      T, f_liq, dro_vapdT);
   k_eff = icemodel.column.bulk_thermal_conductivity(T, f_ice, f_liq, 0);
   [k_eff_faces, k_vap_faces, q_vap_deferred, U_vap_faces, L_vap_faces] = ...
      icemodel.column.vapor_transport_terms( ...
      T, f_ice, f_liq, k_eff, ro_vap, dro_vapdT, De, delz, fn, f_res_por);

   % The exact face flux still runs against the temperature gradient at the
   % wet/dry faces: that is the Bergeron-Findeisen direction the design must
   % keep. Losing it would mean conjugacy was removed rather than treated.
   [L_faces, d_T] = faceLatentHeat(T, f_ice, f_liq, ro_vap, f_res_por);

   % The function's own returned face donor latent heat must match this
   % independent reconstruction exactly, at every face this worst-region
   % fixture exercises.
   testCase.verifyEqual(L_vap_faces, L_faces, 'AbsTol', 0);
   interior = 2:JJ;
   q_total = k_vap_faces(interior) .* d_T(interior) ./ delz(interior) ...
      + q_vap_deferred(interior);
   testCase.assertTrue(any( ...
      q_total .* d_T(interior) < 0 & d_T(interior) ~= 0), ...
      'the Bergeron-Findeisen branch does not occur');

   % The identity holds face by face in the worst region too.
   expected = L_faces(interior) .* U_vap_faces(interior);
   testCase.verifyEqual(q_total, expected, 'RelTol', 1e-12);

   % Assemble the coupled system the way solve_column_enthalpy does and
   % assert on its coefficients. The column sits outside the melt zone, so
   % the melt-zone transform leaves the conductances unscaled.
   f_wat = icemodel.column.water_fraction(f_ice, f_liq);
   [~, dHdT, dFdT] = icemodel.column.bulk_enthalpy( ...
      T, f_ice, f_liq, f_wat, ro_vap, dro_vapdT);
   [aN, aP, aS] = icemodel.column.assemble_enthalpy_system( ...
      T, f_ice, f_liq, dHdT, dFdT, dro_vapdT, zeros(JJ, 1), ...
      zeros(JJ, 1), zeros(JJ, 1), k_eff_faces, delz, dz, 900, ...
      Tf - 10, 0, 0, 1, q_vap_deferred);

   testCase.verifyTrue(all(aN >= 0));
   testCase.verifyTrue(all(aS >= 0));
   testCase.verifyTrue(all(aP > 0));

   % Diagonal dominance itself: the diagonal exceeds the off-diagonal sum
   % by the positive storage coefficient in every row.
   testCase.verifyTrue(all(aP - aN - aS > 0));

   % The bead's acceptance sweep, compacted to its corners and midpoints:
   % temperature gradients 0.5 to 25 K, surface depressions 0 to 20 K,
   % liquid fractions 0 to 0.10, ice fractions 0.30 to 0.85, wet upper
   % half over dry lower half. The matrix part is a positive tangent at
   % every face by construction, so no state in the sweep may produce a
   % negative off-diagonal or lose dominance. A 2 K margin keeps every
   % state below the melt zone: the enthalpy transform rescales the
   % coefficients inside it, and the faces the bead measured are the
   % sub-freezing ones.
   dt = 900;
   solver_dt = 0.1;
   solver = 1;
   tol = 1e-4;
   maxiter = 100;
   n_states = 0;
   for grad = [0.5, 5, 25]
      for depression = [0, 10, 20]
         for f_liq_wet = [0, 0.05, 0.10]
            for f_ice_s = [0.30, 0.85]
               n_states = n_states + 1;
               T = (Tf - 2 - depression - grad) + linspace(0, grad, JJ)';
               f_ice = f_ice_s * ones(JJ, 1);
               f_liq = zeros(JJ, 1);
               f_liq(1:JJ / 2) = f_liq_wet;

               [ro_vap, dro_vapdT] = ...
                  icemodel.vapor.saturation_vapor_density(T, f_liq);
               [~, De] = icemodel.vapor.vapor_thermal_conductivity( ...
                  T, f_liq, dro_vapdT);
               k_eff = icemodel.column.bulk_thermal_conductivity( ...
                  T, f_ice, f_liq, 0);
               [k_eff_faces, k_vap_faces, q_vap_deferred, U_vap_faces] = ...
                  icemodel.column.vapor_transport_terms( ...
                  T, f_ice, f_liq, k_eff, ro_vap, dro_vapdT, De, delz, ...
                  fn, f_res_por);

               % Every face quantity must remain finite throughout the
               % acceptance matrix, and the implicit vapor share must stay
               % nonnegative.
               face_values = [k_eff_faces; k_vap_faces; ...
                  q_vap_deferred; U_vap_faces];
               testCase.assertTrue(all(isfinite(face_values)), sprintf( ...
                  'nonfinite face quantity in acceptance state %d', n_states));
               testCase.assertTrue(all(k_vap_faces >= 0), sprintf( ...
                  'negative vapor matrix term in acceptance state %d', ...
                  n_states));

               f_wat = icemodel.column.water_fraction(f_ice, f_liq);
               [~, dHdT, dFdT] = icemodel.column.bulk_enthalpy( ...
                  T, f_ice, f_liq, f_wat, ro_vap, dro_vapdT);
               [aN, aP, aS] = icemodel.column.assemble_enthalpy_system( ...
                  T, f_ice, f_liq, dHdT, dFdT, dro_vapdT, zeros(JJ, 1), ...
                  zeros(JJ, 1), zeros(JJ, 1), k_eff_faces, delz, dz, dt, ...
                  T(1), 0, 0, 1, q_vap_deferred);

               testCase.assertTrue(all(aN >= 0) && all(aS >= 0) ...
                  && all(aP > 0) && all(aP - aN - aS > 0), sprintf( ...
                  'dominance lost at grad=%g depression=%g f_liq=%g f_ice=%g', ...
                  grad, depression, f_liq_wet, f_ice_s));

               % Exercise the production nonlinear solve on the same state.
               % The Dirichlet surface equals the top-node temperature and
               % there are no external sources. A 0.1 s substep isolates
               % numerical health: these coefficient-stress corners include
               % deliberately cold liquid states that the adaptive driver
               % would not ask to advance by the 900 s assembly interval.
               [T_solve, f_ice_solve, f_liq_solve, k_eff_solve, ...
                  U_vap_solve, ~, ok_solve, ~, ~, ~] = ...
                  icemodel.column.solve_column_enthalpy( ...
                  T(1), T, f_ice, f_liq, 0, 0, zeros(JJ, 1), ...
                  zeros(JJ, 1), dz, delz, fn, solver_dt, solver, tol, ...
                  maxiter, ...
                  1, false, 10, false, f_res_por);
               testCase.assertTrue(ok_solve, sprintf( ...
                  'column solve failed in acceptance state %d', n_states));
               accepted_values = [T_solve; f_ice_solve; f_liq_solve; ...
                  k_eff_solve; U_vap_solve];
               testCase.assertTrue(all(isfinite(accepted_values)), sprintf( ...
                  'nonfinite accepted solve state %d', n_states));
            end
         end
      end
   end
   testCase.verifyEqual(n_states, 54);
end

function test_the_donor_phase_matches_the_mass_applier_in_the_band(testCase)
   % In the band where the residual floor and the fixed storage threshold
   % disagree, the mass applier classifies a cell wet while
   % icemodel.vapor.latent_enthalpy_switch still says dry. The face latent
   % heat must follow the applier's predicate,
   % icemodel.column.vapor_exchange_is_wet, so the energy the solve moves
   % and the mass the applier moves use one latent heat.

   [Ls, Lv, Tf] = icemodel.physicalConstant('Ls', 'Lv', 'Tf');
   f_res_por = 0.02;

   % The band cell the bead measured: wet to the applier, dry to the
   % storage switch.
   f_ice_band = 0.90;
   f_liq_band = 0.0135;
   testCase.assertTrue(icemodel.column.vapor_exchange_is_wet( ...
      f_ice_band, f_liq_band, f_res_por));
   testCase.assertEqual(icemodel.vapor.latent_enthalpy_switch(f_liq_band), Ls);

   % Two cells, the band cell warmer and on top so it donates the flux.
   T = [Tf - 3; Tf - 5];
   f_ice = f_ice_band * ones(2, 1);
   f_liq = [f_liq_band; 0];
   delz = [0.04; 0.04; 0.04];
   fn = [0.5; 0.5; 0.5];
   [ro_vap, De, dro_vapdT] = nodeQuantities(T, f_liq);

   % The band cell holds the higher vapor density, so it is the donor.
   testCase.assertGreaterThan(ro_vap(1), ro_vap(2));

   k_eff = icemodel.column.bulk_thermal_conductivity(T, f_ice, f_liq, 0);
   [~, k_vap_faces, ~, ~, L_vap_faces] = ...
      icemodel.column.vapor_transport_terms( ...
      T, f_ice, f_liq, k_eff, ro_vap, dro_vapdT, De, delz, fn, f_res_por);

   % The function's own returned face donor latent heat carries the band
   % cell's value directly.
   testCase.verifyEqual(L_vap_faces(2), Lv, 'AbsTol', 0);

   % Recover the latent heat the face carries from the matrix part too: the
   % donor tangent and face diffusivity are known, so L divides out. This is
   % an independent cross-check on the returned value above.
   De_faces = 1.0 ./ ((1.0 - fn) ./ [De(1); De] ...
      + fn ./ [De; De(end)]);
   returned = k_vap_faces(2) / (De_faces(2) * dro_vapdT(1));
   testCase.verifyEqual(returned, Lv, 'RelTol', 1e-12);
end

function test_the_surface_fraction_scales_with_the_supplying_phase(testCase)
   % The same energy demand moves less mass as ice than as liquid, because
   % sublimation costs Ls and evaporation costs Lv. A wet surface must return
   % the demand unchanged, and a dry one must return it smaller by Lv / Ls.

   [Ls, Lv] = icemodel.physicalConstant('Ls', 'Lv');
   d_pevp = -1e-4;
   f_res_por = 0.02;

   % A wet top cell and a dry one, on the predicate that owns the decision.
   f_ice = 0.85;
   f_liq_wet = 0.05;
   f_liq_dry = 0.0;
   testCase.assertTrue(icemodel.column.vapor_exchange_is_wet( ...
      f_ice, f_liq_wet, f_res_por));
   testCase.assertFalse(icemodel.column.vapor_exchange_is_wet( ...
      f_ice, f_liq_dry, f_res_por));

   [wet_liq, wet_ice, wet_f_res] = ...
      icemodel.surface.potential_surface_vapor_exchange( ...
      d_pevp, f_ice, f_liq_wet, f_res_por);
   [dry_liq, dry_ice, dry_f_res] = ...
      icemodel.surface.potential_surface_vapor_exchange( ...
      d_pevp, f_ice, f_liq_dry, f_res_por);

   % With Lv the liquid demand passes through; the dry ice increment scales
   % by Lv/Ls. Each state returns the floor used for phase selection.
   testCase.verifyEqual(wet_liq, d_pevp, 'RelTol', 1e-14);
   testCase.verifyEqual(wet_ice, 0);
   testCase.verifyEqual(dry_liq, 0);
   testCase.verifyEqual(dry_ice / wet_liq, Lv / Ls, 'RelTol', 1e-14);
   testCase.verifyLessThan(abs(dry_ice), abs(wet_liq));
   testCase.verifyGreaterThan(wet_f_res, 0);
   testCase.verifyGreaterThan(dry_f_res, 0);
end

function test_the_surface_fraction_sums_over_substeps(testCase)
   % The driver accumulates this fraction per substep rather than converting
   % one step total, because the phase can change within a step. The helper is
   % linear in the demand, so a sum of substeps taken at one phase equals the
   % single conversion of their total. That is what makes the per-substep
   % accumulation safe for a step whose phase never changes.

   f_ice = 0.85;
   f_liq = 0.0;
   f_res_por = 0.02;
   d_pevp = [-3e-5; -1e-5; -6e-5];

   summed_liq = 0;
   summed_ice = 0;
   for k = 1:numel(d_pevp)
      [d_liq, d_ice] = icemodel.surface.potential_surface_vapor_exchange( ...
         d_pevp(k), f_ice, f_liq, f_res_por);
      summed_liq = summed_liq + d_liq;
      summed_ice = summed_ice + d_ice;
   end
   [returned_liq, returned_ice] = ...
      icemodel.surface.potential_surface_vapor_exchange( ...
      sum(d_pevp), f_ice, f_liq, f_res_por);

   testCase.verifyEqual(summed_liq, returned_liq, 'AbsTol', 0);
   testCase.verifyEqual(summed_ice, returned_ice, 'RelTol', 1e-14);
end

function test_the_surface_flux_is_a_unit_conversion(testCase)
   % The phase correction happens in potential_surface_vapor_exchange, so
   % this carries no latent heat. It converts a fraction to the face mass
   % flux and nothing else, which keeps the driver's accumulator a fraction.

   ro_liq = icemodel.physicalConstant('ro_liq');
   d_vap_sfc = -1e-4;
   dz_top = 0.04;
   dt = 900;

   returned = icemodel.surface.surface_vapor_mass_flux(d_vap_sfc, dz_top, dt);
   expected = d_vap_sfc * ro_liq * dz_top / dt;

   testCase.verifyEqual(returned, expected, 'RelTol', 1e-14);
end

function [L_faces, d_T] = faceLatentHeat(T, f_ice, f_liq, ro_vap, f_res_por)
   %FACELATENTHEAT Rebuild the donor-cell latent heat and face temperature drop.

   [Ls, Lv] = icemodel.physicalConstant('Ls', 'Lv');

   JJ = numel(T);
   T_nodes = [T(1); T; T(JJ)];
   ro_vap_nodes = [ro_vap(1); ro_vap; ro_vap(JJ)];
   f_ice_nodes = [f_ice(1); f_ice; f_ice(JJ)];
   f_liq_nodes = [f_liq(1); f_liq; f_liq(JJ)];

   d_ro_vap = ro_vap_nodes(2:JJ+2) - ro_vap_nodes(1:JJ+1);

   % The conductance multiplies the north-minus-south difference, because the
   % flux is positive downward. Return that orientation, not its negative.
   d_T = T_nodes(1:JJ+1) - T_nodes(2:JJ+2);

   % The donor phase follows the mass applier's predicate, so this helper
   % agrees with the production rule by construction.
   wet_nodes = icemodel.column.vapor_exchange_is_wet( ...
      f_ice_nodes, f_liq_nodes, f_res_por);
   L_nodes = Ls * ones(JJ + 2, 1);
   L_nodes(wet_nodes) = Lv;
   L_north = L_nodes(1:JJ+1);
   L_south = L_nodes(2:JJ+2);
   donor_is_north = d_ro_vap <= 0;
   L_faces = L_south;
   L_faces(donor_is_north) = L_north(donor_is_north);
end

function [ro_vap, De, dro_vapdT] = nodeQuantities(T, f_liq)
   %NODEQUANTITIES Evaluate the accepted-state vapor quantities at the nodes.

   % Evaluate the same density, tangent, and diffusivity pair the column
   % solver uses for one accepted Picard state.
   [ro_vap, dro_vapdT] = ...
      icemodel.vapor.saturation_vapor_density(T, f_liq);
   [~, De] = icemodel.vapor.vapor_thermal_conductivity( ...
      T, f_liq, dro_vapdT);
end

function [T, f_ice, f_liq, delz, fn, f_res_por] = faceFixture(JJ)
   %FACEFIXTURE Return one dry column with a temperature gradient.

   Tf = icemodel.physicalConstant('Tf');
   [~, delz, ~, ~, fn] = icemodel.column.control_volume_mesh(JJ * 0.04, 0.04);
   delz = delz(1:JJ + 1);
   fn = fn(1:JJ + 1);

   % A gradient, so every interior face carries a real flux to check.
   T = (Tf - 8) + linspace(0, 6, JJ)';
   f_ice = 0.6 * ones(JJ, 1);
   f_liq = zeros(JJ, 1);

   % The residual pore fraction the donor-phase predicate needs.
   f_res_por = 0.02;
end
