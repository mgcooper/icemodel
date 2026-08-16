function tests = test_vapor_mass_transfer
   %TEST_VAPOR_MASS_TRANSFER Pin the vapor transport kernel.
   %
   % icemodel.column.vapor_mass_transfer solves Fick's law on the saturation
   % vapor density and grows grains from the result. The model keeps only the
   % grain radius today, so the flux and the mass source leave the kernel
   % untested by any consumer. These cases pin the sign conventions, the
   % boundary treatment, the conservation property of the divergence, and the
   % grain-growth branches, so a later change to the discretization has to
   % move a test rather than pass unnoticed.
   %
   % See also: icemodel.column.vapor_mass_transfer
   tests = functiontests(localfunctions);
end

function test_isothermal_column_moves_no_interior_vapor(testCase)
   % With one temperature everywhere and the surface at that temperature, the
   % saturation vapor density is uniform, so every face gradient is zero. This
   % is the reproduction case the coupled mode must match: no interior flux,
   % and no mass source in any cell.

   [T, f_ice, f_liq, radius, dz, delz, fn] = kernelFixture(6);
   Ts = T(1);

   [~, U_vap_faces, dm_vap] = icemodel.column.vapor_mass_transfer( ...
      T, Ts, f_ice, f_liq, radius, dz, delz, fn, 900);

   testCase.verifyEqual(U_vap_faces, zeros(numel(T) + 1, 1), 'AbsTol', 1e-18);
   testCase.verifyEqual(dm_vap, zeros(numel(T), 1), 'AbsTol', 1e-18);
end

function test_flux_sign_follows_the_temperature_gradient(testCase)
   % Positive flux is downward, into the column. A surface warmer than the
   % column holds more vapor, so vapor moves down and the top face flux is
   % positive. Reversing the gradient reverses the sign.

   [T, f_ice, f_liq, radius, dz, delz, fn] = kernelFixture(6);

   % Warm surface over a cold column drives vapor downward.
   [~, U_warm] = icemodel.column.vapor_mass_transfer( ...
      T, T(1) + 5, f_ice, f_liq, radius, dz, delz, fn, 900);
   testCase.verifyGreaterThan(U_warm(1), 0);

   % A cold surface over the same column drives vapor upward, out of it.
   [~, U_cold] = icemodel.column.vapor_mass_transfer( ...
      T, T(1) - 5, f_ice, f_liq, radius, dz, delz, fn, 900);
   testCase.verifyLessThan(U_cold(1), 0);
end

function test_bottom_face_is_closed_and_divergence_conserves_mass(testCase)
   % The bottom boundary is zero flux, so the column exchanges vapor with the
   % surface alone. The mass source is the flux divergence, so the mass the
   % cells gain must equal the mass that crossed the top face.

   [T, f_ice, f_liq, radius, dz, delz, fn] = kernelFixture(6);
   T = T + linspace(0, 3, numel(T))';

   [~, U_vap_faces, dm_vap] = icemodel.column.vapor_mass_transfer( ...
      T, T(1) + 4, f_ice, f_liq, radius, dz, delz, fn, 900);

   testCase.verifyEqual(U_vap_faces(end), 0);

   % dm_vap is a volumetric rate, so scale by the cell thickness before
   % summing. A closed bottom makes that sum the top face flux exactly.
   testCase.verifyEqual(sum(dm_vap .* dz), U_vap_faces(1), 'RelTol', 1e-12);
end

function test_dry_grain_growth_follows_the_vapor_flux(testCase)
   % Jordan Eq. 33 grows a dry grain from the vapor flux magnitude, so a
   % column with a gradient must grow and a column without one must not.

   [T, f_ice, f_liq, radius, dz, delz, fn] = kernelFixture(6);
   f_liq = zeros(size(f_liq));

   still = icemodel.column.vapor_mass_transfer( ...
      T, T(1), f_ice, f_liq, radius, dz, delz, fn, 900);
   testCase.verifyEqual(still, radius, 'AbsTol', 0);

   grown = icemodel.column.vapor_mass_transfer( ...
      T, T(1) + 5, f_ice, f_liq, radius, dz, delz, fn, 900);
   testCase.verifyGreaterThan(grown(1), radius(1));

   % Jordan's model is monotonic: grains never shrink, whatever the sign of
   % the flux, because growth uses its magnitude.
   shrunk = icemodel.column.vapor_mass_transfer( ...
      T, T(1) - 5, f_ice, f_liq, radius, dz, delz, fn, 900);
   testCase.verifyGreaterThanOrEqual(min(shrunk - radius), 0);
end

function test_wet_grain_growth_uses_the_liquid_branches(testCase)
   % Jordan Eqs. 34a and 34b grow a wet grain from the liquid fraction rather
   % than the vapor flux, and cap the rate at f_liq = 0.09. An isothermal
   % column isolates those branches, because the vapor flux is then zero.

   [T, f_ice, ~, radius, dz, delz, fn] = kernelFixture(3);
   Ts = T(1);

   % Below the cap, more liquid grows the grain faster.
   low = icemodel.column.vapor_mass_transfer( ...
      T, Ts, f_ice, 0.02 * ones(3, 1), radius, dz, delz, fn, 900);
   higher = icemodel.column.vapor_mass_transfer( ...
      T, Ts, f_ice, 0.08 * ones(3, 1), radius, dz, delz, fn, 900);
   testCase.verifyGreaterThan(higher(1), low(1));

   % At and above the cap the rate does not depend on the liquid fraction.
   at_cap = icemodel.column.vapor_mass_transfer( ...
      T, Ts, f_ice, 0.09 * ones(3, 1), radius, dz, delz, fn, 900);
   above_cap = icemodel.column.vapor_mass_transfer( ...
      T, Ts, f_ice, 0.30 * ones(3, 1), radius, dz, delz, fn, 900);
   testCase.verifyEqual(above_cap, at_cap, 'AbsTol', 0);
end

function test_grain_radius_stops_at_the_parameter_maximum(testCase)
   % r_max is the cap the spectral model expects, so growth must clamp there
   % rather than run away.
   %
   % The step below is ten days, far longer than any the model takes. Wet
   % growth runs at g2 = 4e-12, so a grain starting one percent below the
   % 2.5 mm cap needs a step that long to reach it. The case exists to reach
   % the clamp, not to represent a step the solver would take.

   r_max = icemodel.parameterLookup('r_max');
   [T, f_ice, ~, ~, dz, delz, fn] = kernelFixture(3);
   radius = repmat(0.99 * r_max, 3, 1);

   returned = icemodel.column.vapor_mass_transfer( ...
      T, T(1) + 20, f_ice, 0.30 * ones(3, 1), radius, dz, delz, fn, 8.64e5);

   testCase.verifyEqual(returned, repmat(r_max, 3, 1), 'AbsTol', 0);

   % A realistic step leaves the grain below the cap, which shows the clamp
   % is not simply pinning every result to r_max.
   unclamped = icemodel.column.vapor_mass_transfer( ...
      T, T(1) + 20, f_ice, 0.30 * ones(3, 1), radius, dz, delz, fn, 900);
   testCase.verifyLessThan(unclamped(1), r_max);
   testCase.verifyGreaterThan(unclamped(1), radius(1));
end

function test_precomputed_node_quantities_reproduce_the_self_computed(testCase)
   % The coupled path hands the kernel node quantities the caller already
   % holds, so the kernel does not evaluate its own. Passing the same values
   % the kernel would have computed must change nothing at all. This pins
   % that equivalence, not where the caller got the values.

   [T, f_ice, f_liq, radius, dz, delz, fn] = kernelFixture(6);
   T = T + linspace(0, 3, numel(T))';
   Ts = T(1) + 4;

   [r_self, U_self, dm_self] = icemodel.column.vapor_mass_transfer( ...
      T, Ts, f_ice, f_liq, radius, dz, delz, fn, 900);

   ro_vap = icemodel.vapor.saturation_vapor_density(T, f_liq);
   De = icemodel.vapor.vapor_diffusivity(T);
   [r_given, U_given, dm_given] = icemodel.column.vapor_mass_transfer( ...
      T, Ts, f_ice, f_liq, radius, dz, delz, fn, 900, ro_vap, De);

   testCase.verifyEqual(r_given, r_self, 'AbsTol', 0);
   testCase.verifyEqual(U_given, U_self, 'AbsTol', 0);
   testCase.verifyEqual(dm_given, dm_self, 'AbsTol', 0);
end

function test_neumann_top_replaces_the_saturated_ghost_node(testCase)
   % In coupled mode the top face carries the turbulent exchange the surface
   % energy balance already computed. The ghost node is not evaluated, so the
   % surface temperature cannot influence the result.

   [T, f_ice, f_liq, radius, dz, delz, fn] = kernelFixture(6);
   T = T + linspace(0, 3, numel(T))';
   ro_vap = icemodel.vapor.saturation_vapor_density(T, f_liq);
   De = icemodel.vapor.vapor_diffusivity(T);
   dt = 900;

   % The caller supplies a liquid-water volume fraction; the kernel forms the
   % kilogram basis itself.
   d_vap_sfc = -1e-9;
   U_top = icemodel.surface.surface_vapor_mass_flux(d_vap_sfc, dz(1), dt);

   [~, U_faces, dm_vap] = icemodel.column.vapor_mass_transfer( ...
      T, T(1), f_ice, f_liq, radius, dz, delz, fn, dt, ...
      ro_vap, De, d_vap_sfc);

   % The converted flux is the top face, exactly.
   testCase.verifyEqual(U_faces(1), U_top, 'AbsTol', 0);

   % Two very different surface temperatures must give one answer, which is
   % what proves the ghost node is gone rather than merely overwritten after
   % contributing somewhere else.
   [~, U_warm] = icemodel.column.vapor_mass_transfer( ...
      T, T(1) + 40, f_ice, f_liq, radius, dz, delz, fn, dt, ...
      ro_vap, De, d_vap_sfc);
   testCase.verifyEqual(U_warm, U_faces, 'AbsTol', 0);

   % The column still conserves: a closed bottom makes the mass the cells
   % gain equal the mass that crossed the surface face.
   testCase.verifyEqual(sum(dm_vap .* dz), U_top, 'RelTol', 1e-12);
end

function test_neumann_top_leaves_the_interior_faces_diffusive(testCase)
   % Only face 1 changes. Below it the transport stays diffusive, so the
   % interior faces must match the flag-off run on the same column.

   [T, f_ice, f_liq, radius, dz, delz, fn] = kernelFixture(6);
   T = T + linspace(0, 3, numel(T))';
   ro_vap = icemodel.vapor.saturation_vapor_density(T, f_liq);
   De = icemodel.vapor.vapor_diffusivity(T);

   % Pad the flag-off run with node 1, which is what the coupled path does,
   % so the two runs differ in the top face alone.
   [~, U_padded] = icemodel.column.vapor_mass_transfer( ...
      T, T(1), f_ice, f_liq, radius, dz, delz, fn, 900, ro_vap, De, 0);
   [~, U_coupled] = icemodel.column.vapor_mass_transfer( ...
      T, T(1), f_ice, f_liq, radius, dz, delz, fn, 900, ro_vap, De, -3e-6);

   testCase.verifyEqual(U_coupled(2:end), U_padded(2:end), 'AbsTol', 0);
   testCase.verifyNotEqual(U_coupled(1), U_padded(1));
end

function [T, f_ice, f_liq, radius, dz, delz, fn] = kernelFixture(JJ)
   %KERNELFIXTURE Return one small cold column on the production mesh.

   Tf = icemodel.physicalConstant('Tf');
   [dz, delz, ~, ~, fn] = ...
      icemodel.column.control_volume_mesh(JJ * 0.04, 0.04);
   dz = dz(1:JJ);
   delz = delz(1:JJ + 1);
   fn = fn(1:JJ + 1);
   T = (Tf - 5) * ones(JJ, 1);
   f_ice = 0.85 * ones(JJ, 1);
   f_liq = 0.01 * ones(JJ, 1);
   radius = 5e-4 * ones(JJ, 1);
end
