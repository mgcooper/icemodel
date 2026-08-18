function [d_vap, dm_vap, U_vap_faces] = couple_vapor_transport( ...
      ro_vap, De, dz, delz, fn, dt)
   %COUPLE_VAPOR_TRANSPORT Redistribute vapor mass between the cells.
   %
   %  [d_vap, dm_vap, U_vap_faces] = icemodel.column.couple_vapor_transport( ...
   %     ro_vap, De, dz, delz, fn, dt)
   %
   % Returns the mass each cell gains or loses to its neighbours over one
   % substep, as a liquid-water volume fraction.
   %
   % This is the interior transport alone. The surface exchange is not part
   % of it, and the top face is closed here. The two obey different
   % invariants, and mixing them loses mass.
   %
   % The surface exchange conserves ENERGY. The surface energy balance fixes
   % Qe, and the mass follows from whichever phase supplies it: the same Qe
   % sublimates less mass than it evaporates. That is why d_pevp is an energy
   % demand and why icemodel.surface.apply_surface_vapor_exchange may spend
   % part of one demand on liquid and the rest on ice.
   %
   % Interior transport conserves MASS. Fick's law moves a definite mass
   % between cells, and the energy is whatever the phase changes need. A cell
   % that receives vapor receives that mass whether it lands in liquid or in
   % ice.
   %
   % Routing interior transport through the energy-demand path would let a
   % wet cell whose liquid margin runs out spend the remainder at a different
   % latent heat. The mass it applied would then not be the mass that
   % arrived. Nothing would record the difference, because nothing was
   % rejected.
   %
   % The divergence is linear in the faces, so closing the top face here
   % separates the two exactly. The caller applies this on a mass basis and
   % the surface demand on its own energy basis.
   %
   % Conservation. Both boundaries are closed, so the volume-weighted sum
   % over the column is zero: sum(d_vap .* dz) == 0, and likewise
   % sum(dm_vap .* dz) == 0. The flux divergence telescopes to
   % U_vap_faces(1) - U_vap_faces(JJ+1), which both closed boundaries make
   % zero. This moves mass between cells and creates none. The unweighted
   % sum is zero only on a uniform grid.
   %
   % Inputs
   %   ro_vap   - Saturation vapor density at the nodes [kg m-3] (JJ x 1),
   %              from the accepted solve.
   %   De       - Effective vapor diffusivity at the nodes [m2 s-1] (JJ x 1),
   %              from the accepted solve.
   %   dz       - Control-volume thicknesses [m] (JJ x 1).
   %   delz     - Distances between adjacent node centers [m] (JJ+1 x 1).
   %   fn       - Interface interpolation weights (JJ+1 x 1).
   %   dt       - Substep length [s].
   %
   % Outputs
   %   d_vap  - Vapor gained per cell as a liquid-water volume fraction [-]
   %            (JJ x 1). sum(d_vap .* dz) is zero, because both boundaries
   %            are closed.
   %   dm_vap - Volumetric mass source rate [kg m-3 s-1] (JJ x 1).
   %   U_vap_faces - Face mass flux [kg m-2 s-1] (JJ+1 x 1), positive
   %            downward, both boundary faces zero. The caller accumulates
   %            its magnitudes so grain growth consumes the fluxes the
   %            column transported.
   %
   % See also: icemodel.column.vapor_face_quantities,
   %  icemodel.column.vapor_mass_transfer,
   %  icemodel.surface.apply_surface_vapor_exchange,
   %  icemodel.surface.potential_surface_vapor_tendency
   %
   %#codegen

   persistent ro_liq
   if isempty(ro_liq)
      ro_liq = icemodel.physicalConstant('ro_liq');
   end

   JJ = numel(ro_vap);

   % Interior faces come from the one shared face rule. The mass this moves
   % and the energy the solve moves are built from the same quantities.
   % Padding the top with node 1 closes that face, which leaves the surface
   % exchange to its own path.
   U_vap_faces = icemodel.column.vapor_face_quantities( ...
      ro_vap, ro_vap(1), De, De(1), delz, fn);

   % Net flux into each control volume [kg m-3 s-1]. Both boundaries are
   % closed, so sum(dm_vap .* dz) is zero.
   dm_vap = (U_vap_faces(1:JJ) - U_vap_faces(2:JJ+1)) ./ dz;

   % Express the mass as the liquid-water volume fraction the applier reads.
   d_vap = dm_vap * dt / ro_liq;
end
