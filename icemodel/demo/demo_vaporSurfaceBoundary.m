%[text] # The surface boundary of the coupled vapor model
%[text] This demo explains why the interior vapor transport operator closes its top face, and shows that the surface still exchanges vapor mass and latent energy on every substep.
%[text] "Closed" describes where the surface exchange is computed and applied. It does not mean the surface blocks vapor.
%[text] The short version:
%[text] - Interior faces (2..JJ) carry diffusive vapor transport: a matrix conductance $k\_{\\mathrm{vap}}$ plus a deferred correction, and the conjugate mass flux $U\_{\\mathrm{vap}}$.
%[text] - Face 1 is the snow-atmosphere interface. Its vapor exchange is turbulent, not pore-diffusive, so the SEB computes it (`Qe`) and the surface exchange step applies the mass ($d\_{\\mathrm{pevp}}$). The interior operator carries no $k\_{\\mathrm{vap}}$, no $U\_{\\mathrm{vap}}$, and no deferred flux at face 1, so the exchange enters the model exactly once.
%[text] - This is the standard Patankar treatment of a flux boundary condition: at a boundary, the known (or linearized) boundary flux replaces the face conductance. \
%%
%[text] ## 1. A small column and its faces
%[text] Build a short column on the production mesh and evaluate the coupled face terms. Face indexing: face `j` sits on top of node `j`, so face 1 is the surface, faces 2..JJ are interior, and face JJ+1 is the bottom.
JJ = 8;
dz0 = 0.04;
[dz, delz, ~, ~, fn] = icemodel.column.control_volume_mesh(JJ * dz0, dz0);
dz = dz(1:JJ);
delz = delz(1:JJ+1);
fn = fn(1:JJ+1);

[Tf, Ls, Lv, ro_ice, ro_liq] = ...
   icemodel.physicalConstant('Tf', 'Ls', 'Lv', 'ro_ice', 'ro_liq');

% A cold, dry column with a temperature gradient that drives upward vapor
% motion (warm below, cold above), typical of winter firn.
T = (Tf - 12) + linspace(0, 6, JJ)';
f_ice = 0.55 * ones(JJ, 1);
f_liq = zeros(JJ, 1);
f_res_por = 0.02;

% Node quantities. K_EFF is the vapor-free node conductivity. The coupled
% scheme omits k_vap from node conductivities because vapor transport is
% computed on faces. U_vap is a face quantity. Computing face k_vap from the
% same De, donor latent heat, and ro_vap difference couples the energy and mass
% transfers. The matrix plus deferred energy transfer equals L times the
% requested mass transfer. When no storage limit binds, the applied mass
% matches the request.
%
% ro_vap is exponential in T, so the exact face flux needs the secant ro_vap
% difference across the face. Putting node k_vap into k_eff would instead use
% the node tangent d(ro_vap)/dT. The face construction combines the tangent
% matrix term with a deferred correction and keeps the coefficients positive.
%
% A needle-probe comparison uses the nodal conductivity
% k_eff + (1 - f_ice - f_liq) .* k_vap.
[ro_vap, dro_vapdT] = icemodel.vapor.saturation_vapor_density(T, f_liq);
[~, De] = icemodel.vapor.vapor_thermal_conductivity(T, f_liq, dro_vapdT);
k_eff = icemodel.column.bulk_thermal_conductivity(T, f_ice, f_liq, 0);

[k_eff_faces, k_vap_faces, q_vap_deferred_faces, U_vap_faces, ...
   L_vap_faces] = icemodel.column.vapor_transport_terms(T, f_ice, ...
   f_liq, k_eff, ro_vap, dro_vapdT, De, delz, fn, f_res_por);

face_table = table((1:JJ+1)', k_eff_faces, k_vap_faces, U_vap_faces, ...
   'VariableNames', {'face', 'k_eff_face', 'k_vap_face', 'U_vap'});
disp(face_table)
%%
%[text] Face 1 and face JJ+1 carry zero vapor conductance and zero vapor mass flux; every interior face carries both. The combined face conductance at face 1 is the vapor-free $k\_{\\mathrm{eff}}$`(1)` alone (the harmonic mean of the endpoint with itself):
assert(k_vap_faces(1) == 0 && k_vap_faces(JJ+1) == 0)
assert(U_vap_faces(1) == 0 && U_vap_faces(JJ+1) == 0)
assert(q_vap_deferred_faces(1) == 0 && q_vap_deferred_faces(JJ+1) == 0)
assert(abs(k_eff_faces(1) - k_eff(1)) < 1e-14)
assert(all(k_vap_faces(2:JJ) > 0))
%%
%[text] ## 2. What face 1 actually carries
%[text] The enthalpy solve connects the surface temperature $T\_{\\mathrm{sfc}}$ to node 1 through the half-cell conductance
%[text] $a\_1 = \\frac{k\_{\\mathrm{eff}}(1)}{\\Delta z\_1}, \\qquad \\Delta z\_1 = \\frac{dz(1)}{2}$
%[text] and the SEB enters the same node through the Robin linearization $F\_\\mathrm{c} + F\_\\mathrm{p} \* T\_\\mathrm{sfc}$. $F\_\\mathrm{c}$ and $F\_\\mathrm{p}$ linearize the full surface energy balance, including the turbulent latent flux $Q\_\\mathrm{e}$. Eliminating $T\_{\\mathrm{sfc}}$ between the surface balance and the node-1 equation gives the closure the couplers use:
%[text] $T\_{\\mathrm{sfc}} = \\frac{F\_\\mathrm{c} + a\_1 T\_1}{a\_1 - F\_\\mathrm{p}}$
%[text] So the surface exchanges latent energy on every iteration, through $Q\_\\mathrm{e}$ inside $F\_\\mathrm{c}$/$F\_\\mathrm{p}$, never through a face conductance. $T\_{\\mathrm{sfc}}$ is a real unknown on the face, eliminated algebraically.
assert(abs(delz(1) - dz(1) / 2) < 1e-14)
a1 = k_eff_faces(1) / delz(1);

% conductive_heat_flux is the same half-cell link, written as a flux. The
% two agree exactly because both use the vapor-free k_eff(1).
T_sfc = Tf - 15;
Qc = icemodel.surface.conductive_heat_flux(k_eff, T, dz, T_sfc);
assert(abs(Qc - a1 * (T(1) - T_sfc)) < 1e-10)
fprintf('a1 = %.3f W m-2 K-1, Qc = %.2f W m-2\n', a1, Qc)
%%
%[text] ## 3. Why no diffusive vapor conductance at face 1
%[text] ### 3a. The transport physics changes at the surface
%[text] $k\_{\\mathrm{vap}}$ parameterizes molecular diffusion of saturated vapor through stagnant pore air. Above the surface the air is not stagnant pore space: the surface layer is turbulent, and the exchange follows the bulk aerodynamic law the SEB computes. Compare the two conductances for the air between the surface and a 2 m measurement height:
kappa = 0.4;
z_ref = 2.0;
z0 = 1e-3;
wspd = (0.5:0.5:10)';

% Neutral bulk aerodynamic conductance [m s-1] (illustrative form).
g_turb = kappa^2 * wspd ./ (log(z_ref / z0))^2;

% Molecular diffusive conductance of the same air layer [m s-1].
g_diff = De(1) / z_ref;

figure
semilogy(wspd, g_turb, '-', 'LineWidth', 1.5)
hold on
yline(g_diff, '--', 'molecular diffusion', 'LineWidth', 1.5)
xlabel('wind speed [m s^{-1}]')
ylabel('surface-air vapor conductance [m s^{-1}]')
title('Turbulent vs molecular surface-air exchange')
legend('bulk aerodynamic (SEB)', 'Location', 'southeast')

fprintf(['At 5 m s-1 wind, turbulent exchange is %.0f times the ', ...
   'molecular rate.\n'], g_turb(wspd == 5) / g_diff)
%[text] A diffusive face-1 conductance would model the atmosphere as stagnant pore air. It would use the wrong transport law and the wrong magnitude. The SEB bulk formula uses the real turbulent transfer coefficient.
%%
%[text] ### 3b. One exchange with the atmosphere, applied once
%[text] `Qe` already IS the surface-to-atmosphere latent exchange. An open face 1 diffusing to an air value would exchange with the outside in parallel with `Qe` and count the same physical flux twice. The mass side mirrors this: the surface exchange step applies the demand derived from the converged `Qe` ($d\_{\\mathrm{pevp}}$), so a nonzero $U\_{\\mathrm{vap}}$`(1)` would move the same mass twice. Note the scope: this argument rules out a second SURFACE-ATMOSPHERE path. It does not by itself rule out a vapor term on the node-1-to-skin link, which is the next section's question.
%%
%[text] ### 3c. But vapor also diffuses between node 1 and the skin
%[text] The sharpest form of the question: for heat, the model keeps BOTH the turbulent skin-atmosphere sensible heat flux `Qh` and the diffusive node-1-to-skin conductive flux `Qc` - they are segments on opposite sides of the interface, not double counting. Vapor is physically parallel: pore vapor diffuses from node 1 toward the skin while the skin exchanges turbulently with the air. An exact skin balance would include that pore-vapor term, and it would NOT double count `Qe`. So why is the half-cell link $a\_1$ vapor-free?
%[text] The reason is mass. Conduction carries energy without mass, so the massless skin can receive `Qc` from below and hand `Qh` to the air; nothing else is implied. Vapor carries energy WITH mass. The model's skin has no mass store, so a node-1-to-skin vapor flux has no reservoir to land in: every kilogram the surface exchanges must come from cell 1 either way. The coupled scheme therefore routes the whole surface vapor exchange through cell 1 directly ($d\_{\\mathrm{pevp}}$), and the conjugacy rule (every face that moves vapor energy moves L times that mass, see section 6) then forces the energy to take the same route. A $k\_{\\mathrm{vap}}$ inside $a\_1$ would move latent energy across a face that can move no mass, recreating the energy-without-mass inconsistency the coupled design exists to remove.
%[text] Nothing is lost by this routing. Vapor that sublimates near node 1 and deposits at the skin is redistribution INSIDE the top half cell, below grid resolution; the cell-1 enthalpy accounts the phase change wherever the mass budget puts it. The one thing the choice does affect is the $T\_{\\mathrm{sfc}}$ diagnosis, because $a\_1$ sets how tightly the skin tracks node 1. That effect is percent-level:
[k_vap_node, ~] = icemodel.vapor.vapor_thermal_conductivity( ...
   T, f_liq, dro_vapdT);
f_air = 1 - f_ice - f_liq;
fprintf('vapor contribution at node 1 / k_eff(1) = %.1f%%\n', ...
   100 * f_air(1) * k_vap_node(1) / k_eff(1))
%[text] And in the top millimeters the stagnant-saturated-pore assumption behind $k\_{\\mathrm{vap}}$ is at its weakest (wind pumping ventilates the surface), while the SEB already prescribes the skin's vapor state as saturated at $T\_{\\mathrm{sfc}}$. Putting a pore-diffusion closure and the SEB closure on the same thin region would stack two models of one process. If skin-interior vapor exchange ever needs to be resolved, the consistent extension is a true surface node with mass storage, not a lone $k\_{\\mathrm{vap}}$ in $a\_1$.
%%
%[text] ## 4. The mass side: a Neumann boundary by operator splitting
%[text] With both boundary faces closed, the requested interior face transfers conserve column mass. The dz-weighted (areal) increments sum to zero exactly:
dt = 900;
d_vap_faces = U_vap_faces(2:JJ) * dt / ro_liq;
d_north = -d_vap_faces ./ dz(1:JJ-1);
d_south = d_vap_faces ./ dz(2:JJ);
d_interior = zeros(JJ, 1);
d_interior(1:JJ-1) = d_interior(1:JJ-1) + d_north;
d_interior(2:JJ) = d_interior(2:JJ) + d_south;
fprintf('areal sum of interior increments: %.3e (exactly zero)\n', ...
   sum(d_interior .* dz))
assert(abs(sum(d_interior .* dz)) < 1e-18)
%%
%[text] When both sides of each interior transfer are applied, the boundary mass flux is the only term that changes the column total. It arrives as the SEB-derived demand: `Qe -> d_pevp -> apply_surface_vapor_exchange`. A cell storage limit can reject one side of an interior transfer, so the applied transport increments can have a nonzero total; the transport budget records that change. The surface boundary is a Neumann (specified-flux) condition implemented as an operator split: interior transport with closed boundary faces, plus one boundary source.
%%
%[text] ## 5. Patankar's treatment of the boundary
%[text] Patankar (1980), *Numerical Heat Transfer and Fluid Flow*:
%[text] - Section 4.2-3, "The Interface Conductivity" (Eqs. 4.8-4.9, about pp. 44-45): the harmonic-mean face conductance applies between two grid points inside the domain. `vapor_transport_terms` uses exactly this rule, fn-weighted, for interior faces. \
%[text] - Section 4.2-6, "Boundary Conditions" (about pp. 46-48): at a boundary, integrate over the half control volume and insert the boundary flux directly. For a given flux, or a flux written with a transfer coefficient as $q\_B = h(T\_\\infty - T\_B)$, no interface conductivity is constructed at the boundary face. The flux replaces the conductance-times-gradient product. \
%[text] The SEB linearization is Patankar's third-kind (transfer-coefficient) condition with `Fc + Fp * T_sfc` in place of $h(T\_\\infty - T\_B)$: the radiative, sensible, and latent exchanges all fold into `Fc` and `Fp`. The half-cell link $a\_1 = k(1)/(dz/2)$ is the conductance between the boundary point and the adjacent grid point. It uses the cell's own conductivity, with no harmonic mean, because only one material sits between the surface and node 1. And it is the VAPOR-FREE conductivity because the vapor part of the surface exchange is already inside the boundary flux (`Qe`); pore diffusion ends where the pore space ends.
%%
%[text] ## 6. Consistency with the discretized PDE
%[text] The solve discretizes one enthalpy conservation law,
%[text] $\\frac{\\partial H}{\\partial t} = -\\frac{\\partial}{\\partial z}\\left(q\_{cond} + q\_{vap}\\right)$
%[text] with $q\_{vap} = -L\\,D\_e\\,\\frac{\\partial \\rho\_v}{\\partial z}$ on interior faces and the SEB as the top boundary flux. The interior face terms are conjugate: the matrix term plus deferred correction reconstructs the face latent heat times the requested mass flux. When no storage limit binds, the split step applies that requested mass.
T_pad = [T(1); T; T(JJ)];
d_T = T_pad(1:JJ+1) - T_pad(2:JJ+2);
Q_vap_reconstructed = k_vap_faces .* d_T ./ delz + q_vap_deferred_faces;
Q_vap_conjugate = L_vap_faces .* U_vap_faces;
assert(max(abs(Q_vap_reconstructed - Q_vap_conjugate)) < 1e-12)
fprintf('max |Q_vap - L*U_vap| over all faces: %.2e W m-2\n', ...
   max(abs(Q_vap_reconstructed - Q_vap_conjugate)))
%%
%[text] Natural follow-up questions:
%[text] - Does interior transport change the top node temperature? Yes. Face 2 carries $k\_{\\mathrm{vap}}$ and the deferred flux into the node-1 equation, and the transported mass changes $f\_{\\mathrm{ice}}$/$f\_{\\mathrm{liq}}$, which changes the next substep's properties.
%[text] - Does it change the top node $k\_{\\mathrm{eff}}$? No node $k\_{\\mathrm{eff}}$ contains vapor in the coupled scheme; vapor is a face quantity. The archived `ice2.k_eff` is the vapor-free node conductivity.
%[text] - Which faces include $k\_{\\mathrm{vap}}$? Faces 2..JJ. Face 1 (the surface) and face JJ+1 (the bottom) do not. \
%%
%[text] ## 7. `Qe` as the top boundary condition
%[text] The converged `Qe` defines the potential surface energy demand. `potential_surface_vapor_demand` converts `Qe` to `d_pevp` once per accepted substep. `apply_surface_vapor_exchange` applies the phase changes after the coupler converges; its applied result gives the realized mass flux when a storage limit rejects part of the demand.
%[text] `Fc`, `Fp`, and `Qe` update during every outer coupler iteration. The phase fractions update once per accepted substep, outside the Picard solve. Compare the potential and applied surface-budget terms to quantify vapor demand that the column could not accept.

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
