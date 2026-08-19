function ledger = accumulate_redistribution_budget( ...
      ledger, solid_r, liquid_r, T, f_ice, f_liq, dz, d_sbl_err_cpl)
   %ACCUMULATE_REDISTRIBUTION_BUDGET Record interior vapor transport energy.
   %
   %  ledger = icemodel.column.accumulate_redistribution_budget( ...
   %     ledger, solid_r, liquid_r, T, f_ice, f_liq, dz, d_sbl_err_cpl)
   %
   % Call this once per accepted substep, right after
   % icemodel.column.apply_vapor_transfer realizes the node-wise transport in
   % the cells. SOLID_R and LIQUID_R are the storage totals from just before that
   % step. D_SBL_ERR_CPL is the phase-resolved shortfall converted to an
   % equivalent ice-fraction energy record.
   %
   % Interior transport conserves mass: it moves vapor between cells and
   % creates none. It does not conserve the per-phase storage totals. A
   % cell that gives up vapor from a liquid film loses liquid storage, and
   % a cell that takes it into ice gains solid storage, so the solid and
   % liquid totals move in opposite directions while their sum holds.
   %
   % The transport runs after icemodel.column.accumulate_vapor_budget takes
   % its storage baseline, so no interior term reaches the surface vapor
   % channels, and the surface closure identity holds without a
   % redistribution correction:
   %
   %   vapor_potential = ro_liq * (Ls * vapor_solid + Lv * vapor_liquid
   %                     + Lv * condensation_overflow) + unapplied_vapor
   %
   % The per-phase storage closures use the increments recorded here. The
   % solid storage delta closes against the phase, surface-vapor, remesh, and
   % redistribution-solid terms. The liquid storage delta closes against the
   % corresponding liquid terms.
   %
   % The shortfall channel records increments rejected by per-cell limits. A
   % binding clamp breaks transport mass closure: the donor applies its
   % increment while the receiver cannot. The shortfall keeps this mismatch
   % visible. Cross-phase transfer populates the phase channels. A bound clamp
   % populates the shortfall channel. Unconstrained same-phase redistribution
   % can have nonzero face flux while all six channels remain zero.
   %
   % Inputs
   %   ledger            - Forcing-step ledger.
   %   solid_r, liquid_r - Storage totals [m w.e.] from before the
   %                       redistribution step.
   %   T, f_ice, f_liq   - Column state after the redistribution step.
   %   dz                - Control-volume thickness [m].
   %   d_sbl_err_cpl     - Per-cell unapplied transport [-] (JJ x 1), in
   %                       ice-fraction units from
   %                       icemodel.column.couple_vapor_step.
   %
   % Output
   %   ledger            - Ledger with the per-phase redistribution
   %                       increments and the interior shortfall
   %                       accumulated.
   %
   % See also: icemodel.column.apply_vapor_transfer,
   %  icemodel.column.accumulate_vapor_budget
   %
   %#codegen

   persistent Ls ro_ice
   if isempty(Ls)
      [Ls, ro_ice] = icemodel.physicalConstant('Ls', 'ro_ice');
   end

   [solid_after, liquid_after] = ...
      icemodel.column.integrate_column_budget(T, f_ice, f_liq, dz);

   redistribution_solid = solid_after - solid_r;
   redistribution_liquid = liquid_after - liquid_r;

   % Add rather than assign: a forcing step can take several accepted
   % substeps and each one redistributes.
   ledger.mass_budget_vapor_redistribution_solid_mwe = ...
      ledger.mass_budget_vapor_redistribution_solid_mwe ...
      + redistribution_solid;
   ledger.mass_budget_vapor_redistribution_liquid_mwe = ...
      ledger.mass_budget_vapor_redistribution_liquid_mwe ...
      + redistribution_liquid;

   % Keep the absolute totals so exchanges of opposite sign within one
   % forcing step cannot cancel to a signed zero.
   ledger.mass_budget_vapor_redistribution_solid_gross_mwe = ...
      ledger.mass_budget_vapor_redistribution_solid_gross_mwe ...
      + abs(redistribution_solid);
   ledger.mass_budget_vapor_redistribution_liquid_gross_mwe = ...
      ledger.mass_budget_vapor_redistribution_liquid_gross_mwe ...
      + abs(redistribution_liquid);

   % The shortfall the per-cell limits rejected, scaled cell by cell before
   % the sum on the same energy basis the surface unapplied channel uses.
   % This lives in its own channel pair so the surface closure identity
   % never carries an interior term.
   unapplied_weighted = ro_ice * Ls * d_sbl_err_cpl(:) .* dz(:);
   ledger.mass_budget_vapor_redistribution_unapplied_j_m2 = ...
      ledger.mass_budget_vapor_redistribution_unapplied_j_m2 ...
      + sum(unapplied_weighted);

   % Magnitudes per cell before the sum, so a rejecting cell and a starved
   % cell of opposite sign in one substep cannot cancel to a signed zero.
   ledger.mass_budget_vapor_redistribution_unapplied_gross_j_m2 = ...
      ledger.mass_budget_vapor_redistribution_unapplied_gross_j_m2 ...
      + sum(abs(unapplied_weighted));
end
