function [f_ice, f_liq, d_vap_faces, vapor_solid, vapor_liquid, ledger] = ...
      couple_vapor_step(T, f_ice, f_liq, f_ice_solve, f_liq_solve, ...
      U_vap_faces, d_vap_faces, vapor_solid, vapor_liquid, dz, dt, ...
      f_ice_min, f_res_por, ledger, use_mass_budget)
   %COUPLE_VAPOR_STEP Apply accepted interior vapor transport once.
   %
   %  [f_ice, f_liq, d_vap_faces, vapor_solid, vapor_liquid, ledger] = ...
   %     icemodel.column.couple_vapor_step(T, f_ice, f_liq, ...
   %     f_ice_solve, f_liq_solve, U_vap_faces, d_vap_faces, ...
   %     vapor_solid, vapor_liquid, dz, dt, f_ice_min, f_res_por, ...
   %     ledger, use_mass_budget)
   %
   % Call this once per accepted substep after the surface budget closes.
   % At this point, d_liq records changes between the checkpoint and the
   % surface mass-balance call, so earlier transport would be read as melt or
   % refreezing. The vapor storage baseline covers only surface exchange, so
   % earlier transport would be scored as surface exchange. Shortfall here
   % belongs to redistribution, not surface closure.
   % U_VAP_FACES is the mass flux accepted by the enthalpy solve [kg m-2
   % s-1], positive downward, with both boundary faces closed. Its divergence
   % supplies the signed liquid-water-equivalent node increments.
   %
   % Each interior face is split into donor-phase removal and matching
   % receiver addition so cross-phase faces retain the donor latent heat. The
   % phase decision uses the solve-state fractions before surface exchange;
   % storage limits use the current state after surface exchange. Face 1 is
   % closed in this interior call; the driver records surface exchange there.
   % Gross face magnitudes feed the grain-growth budget.
   %
   % D_VAP_FACES records gross face throughput for the forcing step [m w.e.].
   % VAPOR_SOLID and VAPOR_LIQUID carry the per-phase storage increments that
   % the follow-on remesh budget consumes.
   %
   % Inputs
   %   T              - Column temperature at the accepted substep [K].
   %   f_ice, f_liq   - Column phase fractions after surface exchange [-].
   %   f_ice_solve    - Ice fraction at the converged solve state [-].
   %   f_liq_solve    - Liquid fraction at the converged solve state [-].
   %   U_vap_faces    - Accepted interior vapor mass flux [kg m-2 s-1].
   %   d_vap_faces    - Running gross face exchange [m w.e.] (JJ+1 x 1).
   %   vapor_solid    - Substep solid storage-change context [m w.e.].
   %   vapor_liquid   - Substep liquid storage-change context [m w.e.].
   %   dz             - Control-volume thickness [m].
   %   dt             - Substep length [s].
   %   f_ice_min      - Minimum retained ice fraction [-].
   %   f_res_por      - Residual liquid-water fraction per pore volume [-].
   %   ledger         - Forcing-step mass and energy ledger.
   %   use_mass_budget - True when the diagnostic ledger is being built.
   %
   % Outputs
   %   f_ice, f_liq   - Phase fractions after the interior exchange [-].
   %   d_vap_faces    - Running record with this substep's face magnitudes.
   %   vapor_solid    - Solid storage-change context with transport
   %                    increments [m w.e.].
   %   vapor_liquid   - Liquid storage-change context with transport
   %                    increments [m w.e.].
   %   ledger         - Ledger with redistribution increments and shortfall.
   %
   % See also: icemodel.column.apply_vapor_transfer,
   %  icemodel.column.accumulate_redistribution_budget
   %
   %#codegen

   persistent ro_liq
   if isempty(ro_liq)
      ro_liq = icemodel.physicalConstant('ro_liq');
   end

   JJ = numel(f_ice);

   % Gross throughput accumulates magnitudes so reversing substeps add.
   d_vap_faces = d_vap_faces + abs(U_vap_faces) * dt / ro_liq;

   % Preserve the donor phase on every face. A node-wise net divergence cannot
   % do this across a wet/dry face because incoming and outgoing donors can
   % differ. The solve-state fractions determine phase routing; the current
   % fractions determine capacity.
   [wet, f_res] = icemodel.column.vapor_exchange_is_wet( ...
      f_ice_solve, f_liq_solve, f_res_por);
   d_vap_liq_nodes = zeros(JJ, 1);
   d_vap_ice_nodes = zeros(JJ, 1);
   d_vap_face = U_vap_faces(2:JJ) * dt / ro_liq;
   donor_is_north = d_vap_face >= 0;
   donor_wet = wet(1:JJ-1);
   wet_south = wet(2:JJ);
   donor_wet(~donor_is_north) = wet_south(~donor_is_north);
   d_liq_north = -d_vap_face ./ dz(1:JJ-1) .* donor_wet;
   d_liq_south = d_vap_face ./ dz(2:JJ) .* donor_wet;
   d_ice_north = -d_vap_face ./ dz(1:JJ-1) .* ~donor_wet;
   d_ice_south = d_vap_face ./ dz(2:JJ) .* ~donor_wet;
   d_vap_liq_nodes(1:JJ-1) = d_vap_liq_nodes(1:JJ-1) + d_liq_north;
   d_vap_liq_nodes(2:JJ) = d_vap_liq_nodes(2:JJ) + d_liq_south;
   d_vap_ice_nodes(1:JJ-1) = d_vap_ice_nodes(1:JJ-1) + d_ice_north;
   d_vap_ice_nodes(2:JJ) = d_vap_ice_nodes(2:JJ) + d_ice_south;

   if use_mass_budget
      [solid_r, liquid_r] = icemodel.column.integrate_column_budget( ...
         T, f_ice, f_liq, dz);
   end

   [f_ice, f_liq, d_vap_liq_unapplied, d_vap_ice_unapplied] = ...
      icemodel.column.apply_vapor_transfer(f_ice, f_liq, ...
      d_vap_liq_nodes, d_vap_ice_nodes, f_ice_min, f_res);

   if use_mass_budget
      [solid_after, liquid_after] = ...
         icemodel.column.integrate_column_budget(T, f_ice, f_liq, dz);
      vapor_solid = vapor_solid + (solid_after - solid_r);
      vapor_liquid = vapor_liquid + (liquid_after - liquid_r);

      % Convert both phase shortfalls to the redistribution ledger's
      % ice-fraction energy basis.
      d_sbl_err_cpl = icemodel.column.vapor_shortfall_ice_equivalent( ...
         d_vap_liq_unapplied, d_vap_ice_unapplied);
      ledger = icemodel.column.accumulate_redistribution_budget( ...
         ledger, solid_r, liquid_r, T, f_ice, f_liq, dz, d_sbl_err_cpl);
   end
end
