function tests = test_vapor_face_discretization
   %TEST_VAPOR_FACE_DISCRETIZATION Verify the shared vapor face rule.
   %
   % The coupled vapor mode moves mass with vapor_face_quantities and energy
   % with vapor_face_conductance. Both build on
   % icemodel.column.vapor_face_diffusivity. The design requires the two to be
   % conjugate: the energy the enthalpy solve moves across a face must equal
   % the mass that crosses it times the latent heat of the phase that supplied
   % it. Nothing else in the suite tests that identity at the face. A column
   % closure test can balance while a face is wrong, because the errors of an
   % adjacent pair cancel in the sum.
   %
   % See also: icemodel.column.vapor_face_quantities,
   %  icemodel.column.vapor_face_conductance,
   %  icemodel.column.vapor_face_diffusivity
   tests = functiontests(localfunctions);
end

function test_face_energy_equals_latent_heat_times_face_mass(testCase)
   % The identity the whole conjugate discretization exists to hold:
   %
   %   k_vap * (T_north - T_south) / delz == L_face * U_vap
   %
   % at every interior face. The conductance is the mass flux rewritten in
   % temperature, so recovering the flux from it must return the flux.

   [T, f_liq, delz, fn] = faceFixture(6);
   [ro_vap, De, dro_vapdT] = nodeQuantities(T, f_liq);

   U_vap_faces = icemodel.column.vapor_face_quantities( ...
      ro_vap, ro_vap(1), De, De(1), delz, fn);
   k_vap_faces = icemodel.column.vapor_face_conductance( ...
      T, f_liq, ro_vap, dro_vapdT, De, fn);

   % Rebuild the latent heat and the temperature difference the conductance
   % used, so the check reads the same faces the functions did.
   [L_faces, d_T] = faceLatentHeat(T, f_liq, ro_vap);

   % Interior faces alone. Both boundary faces are set to zero by design and
   % carry no identity to check.
   interior = 2:numel(T);
   returned = k_vap_faces(interior) .* d_T(interior) ./ delz(interior);
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

   [T, ~, delz, fn] = faceFixture(4);

   % Wet above, dry below, so face 3 spans the phase change.
   f_liq = [0.05; 0.05; 0; 0];
   [ro_vap, De, dro_vapdT] = nodeQuantities(T, f_liq);

   U_vap_faces = icemodel.column.vapor_face_quantities( ...
      ro_vap, ro_vap(1), De, De(1), delz, fn);
   k_vap_faces = icemodel.column.vapor_face_conductance( ...
      T, f_liq, ro_vap, dro_vapdT, De, fn);
   [L_faces, d_T] = faceLatentHeat(T, f_liq, ro_vap);

   % The spanning face takes one of the two latent heats, not a blend.
   [Ls, Lv] = icemodel.physicalConstant('Ls', 'Lv');
   testCase.verifyTrue(L_faces(3) == Ls || L_faces(3) == Lv);

   returned = k_vap_faces(3) * d_T(3) / delz(3);
   expected = L_faces(3) * U_vap_faces(3);
   testCase.verifyEqual(returned, expected, 'RelTol', 1e-12);
end

function test_both_boundary_faces_carry_no_vapor_energy(testCase)
   % The bottom is closed and the surface exchange enters through the SEB's
   % Qe. A conductance at either boundary would count the surface exchange
   % twice, which is what the Neumann boundary decision exists to prevent.

   [T, f_liq, ~, fn] = faceFixture(5);
   [ro_vap, De, dro_vapdT] = nodeQuantities(T, f_liq);

   k_vap_faces = icemodel.column.vapor_face_conductance( ...
      T, f_liq, ro_vap, dro_vapdT, De, fn);

   testCase.verifyEqual(k_vap_faces(1), 0);
   testCase.verifyEqual(k_vap_faces(numel(T) + 1), 0);
end

function test_an_isothermal_column_moves_no_vapor(testCase)
   % Equal temperatures give equal saturation densities, so every face flux
   % is zero. This is the reproduction case the acceptance policy names: an
   % isothermal column must leave the surface path acting alone.

   [~, f_liq, delz, fn] = faceFixture(5);
   Tf = icemodel.physicalConstant('Tf');
   T = (Tf - 5) * ones(5, 1);
   [ro_vap, De] = nodeQuantities(T, f_liq);

   returned = icemodel.column.vapor_face_quantities( ...
      ro_vap, ro_vap(1), De, De(1), delz, fn);
   expected = zeros(6, 1);

   testCase.verifyEqual(returned, expected, 'AbsTol', 0);
end

function test_the_tangent_replaces_an_undefined_secant(testCase)
   % Where a face spans no temperature difference the secant slope divides by
   % zero. Its limit is the tangent, so the conductance must fall back to the
   % donor node's derivative and stay finite.

   [~, f_liq, ~, fn] = faceFixture(4);
   Tf = icemodel.physicalConstant('Tf');

   % Faces 2 and 4 span no temperature difference; face 3 does.
   T = [Tf - 6; Tf - 6; Tf - 2; Tf - 2];
   [ro_vap, De, dro_vapdT] = nodeQuantities(T, f_liq);

   k_vap_faces = icemodel.column.vapor_face_conductance( ...
      T, f_liq, ro_vap, dro_vapdT, De, fn);

   testCase.verifyTrue(all(isfinite(k_vap_faces)));

   % A tangent-slope face is nonzero: ro_vap has a real derivative there even
   % though the neighbours share a temperature.
   testCase.verifyGreaterThan(abs(k_vap_faces(2)), 0);
end

function test_the_face_diffusivity_is_the_harmonic_mean(testCase)
   % Patankar (1980) Eq. 4.9. The harmonic mean makes a face between a
   % conducting cell and an insulating one behave like the insulating one,
   % which is the property an arithmetic mean loses.

   % Two column nodes under a conducting north boundary, so the two interior
   % faces span the conducting/insulating pair in both directions and the
   % third face is the closed bottom, where the deepest node repeats.
   De = [1e-9; 1e-5];
   De_N = 1e-5;
   fn = [0.5; 0.5; 0.5];

   returned = icemodel.column.vapor_face_diffusivity(De, De_N, fn);
   expected = [ ...
      1 / (0.5 / 1e-5 + 0.5 / 1e-9); ...
      1 / (0.5 / 1e-9 + 0.5 / 1e-5); ...
      1e-5];

   testCase.verifyEqual(returned, expected, 'RelTol', 1e-15);

   % The insulating node dominates: the face sits nearer the small value than
   % the arithmetic mean would put it.
   testCase.verifyLessThan(returned(1), mean([De_N; De(1)]));

   % The bottom face repeats the deepest node, so the harmonic mean of a
   % value with itself returns that value.
   testCase.verifyEqual(returned(3), De(2), 'RelTol', 1e-15);
end

function test_both_flux_paths_share_one_face_diffusivity(testCase)
   % The mass flux and the energy conductance must call the same rule. If one
   % gained a weighting factor the other did not, the identity above would
   % fail while every individual function still looked correct.

   [T, f_liq, delz, fn] = faceFixture(5);
   [ro_vap, De] = nodeQuantities(T, f_liq);

   [~, returned] = icemodel.column.vapor_face_quantities( ...
      ro_vap, ro_vap(1), De, De(1), delz, fn);

   % Rebuild the same faces directly from the shared rule. The mass path
   % passes the surface diffusivity as the north value, which this fixture
   % sets to De(1).
   expected = icemodel.column.vapor_face_diffusivity(De, De(1), fn);

   testCase.verifyEqual(returned, expected, 'AbsTol', 0);
end

function test_an_isothermal_wet_dry_face_breaks_conjugacy(testCase)
   % Known defect, bead icemodel-55x. Two neighbours at equal temperature but
   % different phase hold different saturation vapor densities, because the
   % ice and water curves differ. Mass crosses the face while d_T is zero.
   % The conductance falls back to the tangent slope. The assembly multiplies
   % it by the zero temperature difference, so the face carries no energy at
   % all. The identity energy = L * mass fails there.
   %
   % The flux is not expressible as a conductance times a temperature
   % difference at such a face. It belongs in the source term, which is the
   % deferred correction bead icemodel-ziq already owns.
   %
   % This test pins the defect so it cannot be lost. Invert it when the fix
   % lands: the energy the face carries must then equal L_face * U_vap.

   [~, ~, delz, fn] = faceFixture(4);
   Tf = icemodel.physicalConstant('Tf');

   % Isothermal, wet over dry, so face 3 spans the phase change at equal T.
   T = (Tf - 5) * ones(4, 1);
   f_liq = [0.05; 0.05; 0; 0];
   [ro_vap, De, dro_vapdT] = nodeQuantities(T, f_liq);

   U_vap_faces = icemodel.column.vapor_face_quantities( ...
      ro_vap, ro_vap(1), De, De(1), delz, fn);
   k_vap_faces = icemodel.column.vapor_face_conductance( ...
      T, f_liq, ro_vap, dro_vapdT, De, fn);
   [L_faces, d_T] = faceLatentHeat(T, f_liq, ro_vap);

   % Mass crosses the face.
   testCase.assertGreaterThan(abs(U_vap_faces(3)), 0);
   testCase.assertEqual(d_T(3), 0);

   % The energy the assembly carries is zero, and the identity wants L * U.
   carried = k_vap_faces(3) * d_T(3) / delz(3);
   required = L_faces(3) * U_vap_faces(3);
   testCase.verifyEqual(carried, 0);
   testCase.verifyGreaterThan(abs(required), 0);
end

function test_a_wet_dry_face_can_lose_diagonal_dominance(testCase)
   % Pin the known limitation bead icemodel-ziq treats. The conductance goes
   % negative where the vapor flux opposes the temperature gradient, and the
   % enthalpy assembly adds it to an always-positive conduction conductance.
   % In low-density snow the negative part exceeds that conductance, so the
   % coupled matrix is not diagonally dominant there.
   %
   % This test documents the state rather than requiring it. It fails if the
   % negative branch stops occurring, which would mean the conjugacy the
   % design rests on was removed.
   %
   % It does not fail when the treatment lands. The remedy keeps this negative
   % conductance and changes how icemodel.column.assemble_enthalpy_system
   % consumes it, which this test never calls. Bead icemodel-ziq carries the
   % duty to replace this test with one on the assembled coefficient.

   Tf = icemodel.physicalConstant('Tf');
   JJ = 20;
   [~, ~, ~, ~, fn] = icemodel.column.control_volume_mesh(JJ * 0.04, 0.04);
   fn = fn(1:JJ + 1);

   % The worst region the V4 sweep found: low-density snow, shallow gradient,
   % a wet upper half over a dry lower half.
   T = (Tf - 10) + linspace(0, 0.5, JJ)';
   f_ice = 0.30 * ones(JJ, 1);
   f_liq = zeros(JJ, 1);
   f_liq(1:JJ / 2) = 0.05;

   [ro_vap, De] = icemodel.column.accepted_vapor_quantities(T, f_liq);
   [~, dro_vapdT] = icemodel.vapor.saturation_vapor_density(T, f_liq);
   k_vap = icemodel.column.vapor_face_conductance( ...
      T, f_liq, ro_vap, dro_vapdT, De, fn);

   % The assembly's face conductance, built the same way it builds it.
   k_eff = icemodel.column.bulk_thermal_conductivity(T, f_ice, f_liq);
   k_nodes = [k_eff(1); k_eff; k_eff(JJ)];
   k_faces = 1.0 ./ ((1.0 - fn) ./ k_nodes(1:JJ + 1) ...
      + fn ./ k_nodes(2:JJ + 2));

   negative = k_vap < 0;
   testCase.assertTrue(any(negative), ...
      'the Bergeron-Findeisen branch does not occur');

   % The worst ratio exceeds one, which is the loss of dominance itself.
   returned = max(abs(k_vap(negative)) ./ k_faces(negative));
   testCase.verifyGreaterThan(returned, 1);
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

   wet = icemodel.surface.potential_surface_vapor_exchange( ...
      d_pevp, f_ice, f_liq_wet, f_res_por);
   dry = icemodel.surface.potential_surface_vapor_exchange( ...
      d_pevp, f_ice, f_liq_dry, f_res_por);

   % With Lv the latent heats cancel and the demand passes through.
   testCase.verifyEqual(wet, d_pevp, 'RelTol', 1e-14);
   testCase.verifyEqual(dry / wet, Lv / Ls, 'RelTol', 1e-14);
   testCase.verifyLessThan(abs(dry), abs(wet));
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

   summed = 0;
   for k = 1:numel(d_pevp)
      summed = summed + icemodel.surface.potential_surface_vapor_exchange( ...
         d_pevp(k), f_ice, f_liq, f_res_por);
   end
   returned = icemodel.surface.potential_surface_vapor_exchange( ...
      sum(d_pevp), f_ice, f_liq, f_res_por);

   testCase.verifyEqual(summed, returned, 'RelTol', 1e-14);
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

function [L_faces, d_T] = faceLatentHeat(T, f_liq, ro_vap)
   %FACELATENTHEAT Rebuild the donor-cell latent heat and face temperature drop.

   JJ = numel(T);
   T_nodes = [T(1); T; T(JJ)];
   ro_vap_nodes = [ro_vap(1); ro_vap; ro_vap(JJ)];
   f_liq_nodes = [f_liq(1); f_liq; f_liq(JJ)];

   d_ro_vap = ro_vap_nodes(2:JJ+2) - ro_vap_nodes(1:JJ+1);

   % The conductance multiplies the north-minus-south difference, because the
   % flux is positive downward. Return that orientation, not its negative.
   d_T = T_nodes(1:JJ+1) - T_nodes(2:JJ+2);

   L_nodes = icemodel.vapor.latent_enthalpy_switch(f_liq_nodes);
   L_north = L_nodes(1:JJ+1);
   L_south = L_nodes(2:JJ+2);
   donor_is_north = d_ro_vap <= 0;
   L_faces = L_south;
   L_faces(donor_is_north) = L_north(donor_is_north);
end

function [ro_vap, De, dro_vapdT] = nodeQuantities(T, f_liq)
   %NODEQUANTITIES Evaluate the accepted-state vapor quantities at the nodes.

   [ro_vap, De] = icemodel.column.accepted_vapor_quantities(T, f_liq);

   % The tangent slope the conductance falls back on. Take it from the same
   % saturation relation the densities came from, the way
   % icemodel.column.solve_column_enthalpy does.
   [~, dro_vapdT] = icemodel.vapor.saturation_vapor_density(T, f_liq);
end

function [T, f_liq, delz, fn] = faceFixture(JJ)
   %FACEFIXTURE Return one dry column with a temperature gradient.

   Tf = icemodel.physicalConstant('Tf');
   [~, delz, ~, ~, fn] = icemodel.column.control_volume_mesh(JJ * 0.04, 0.04);
   delz = delz(1:JJ + 1);
   fn = fn(1:JJ + 1);

   % A gradient, so every interior face carries a real flux to check.
   T = (Tf - 8) + linspace(0, 6, JJ)';
   f_liq = zeros(JJ, 1);
end
