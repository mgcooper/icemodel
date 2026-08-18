function [ledger, vapor_solid, vapor_liquid] = accumulate_vapor_budget( ...
      ledger, solid_p, liquid_p, T, f_ice, f_liq, dz, d_pevp, d_rof, d_sbl_err)
   %ACCUMULATE_VAPOR_BUDGET Add one substep's surface vapor-exchange increments.
   %
   %  [ledger, vapor_solid, vapor_liquid] = ...
   %     icemodel.column.accumulate_vapor_budget( ...
   %     ledger, solid_p, liquid_p, T, f_ice, f_liq, dz, ...
   %     d_pevp, d_rof, d_sbl_err)
   %
   % Call this once per accepted substep, right after
   % icemodel.column.budget_surface_mass_balance applies the vapor exchange.
   % This keeps the storage change apart from the potential input and from the
   % two amounts the top control volume cannot accept. Those two amounts are
   % condensation above the pore capacity, and deposition the control volume
   % has no energy to apply.
   %
   % Inputs
   %   ledger              - Forcing-step ledger with the phase channels
   %                         already accumulated.
   %   solid_p, liquid_p   - Post-phase storage [m w.e.] from
   %                         icemodel.column.accumulate_phase_budget.
   %   T, f_ice, f_liq     - Column state after the vapor exchange.
   %   dz                  - Control-volume thickness [m].
   %   d_pevp              - Potential vapor-driven top-layer liquid change [-].
   %   d_rof               - Condensation overflow in liquid fraction [-].
   %   d_sbl_err           - Signed unapplied vapor-driven ice change per
   %                         cell [-].
   %
   % Outputs
   %   ledger              - Ledger with the accumulated vapor budgets.
   %   vapor_solid,        - Signed substep vapor increments [m w.e.], needed to
   %   vapor_liquid          close gross storage after remeshing.
   %
   % See also: icemodel.column.accumulate_phase_budget,
   %  icemodel.column.budget_surface_mass_balance
   %
   %#codegen

   persistent Ls Lv ro_ice ro_liq
   if isempty(Ls)
      [Ls, Lv, ro_ice, ro_liq] = icemodel.physicalConstant( ...
         'Ls', 'Lv', 'ro_ice', 'ro_liq');
   end

   % Compute the post-budget_surface_mass_balance total solid and liquid mass in
   % mwe.
   [solid_v, liquid_v] = icemodel.column.integrate_column_budget( ...
      T, f_ice, f_liq, dz);

   % Subtract the pre-budget_surface_mass_balance totals from the post totals on
   % one fixed storage basis, so the increment holds no change of reference
   % density. This call and the call that produced solid_p and liquid_p in
   % accumulate_phase_budget must use the same densities.
   vapor_solid = solid_v - solid_p;
   vapor_liquid = liquid_v - liquid_p;

   % The potential input and the rejected exchanges are energies, not storage.
   % d_pevp is a top-layer fraction scaled by its latent heat. d_sbl_err is
   % one fraction per cell, so integrate it over the column. With a top-cell
   % tendency only the first entry is nonzero, and the integral is then the
   % top-cell term alone.
   % The condensation overflow is a liquid depth the column does not store.
   %
   % Scale each cell before the sum, not the sum afterwards. Multiplication
   % is not associative in floating point. Summing the dz-weighted fractions
   % first and scaling once would move the last bits of a value the
   % surface-only path already produces.
   vapor_potential = ro_liq * Lv * d_pevp * dz(1);

   % A scalar d_sbl_err is the top cell's rejection, which is what the
   % surface-only path produces and what a caller from before the column
   % path still passes. Integrating it against every dz would spread one
   % cell's rejection over the whole column and inflate both channels from
   % the dz(1) term to a sum(dz) one. Weight the scalar form by dz(1) and
   % keep the per-cell sum for the vector form.
   if isscalar(d_sbl_err)
      unapplied_weighted = ro_ice * Ls * d_sbl_err * dz(1);
   else
      unapplied_weighted = ro_ice * Ls * d_sbl_err(:) .* dz(:);
   end
   unapplied_vapor = sum(unapplied_weighted);

   % Take the magnitude per cell before summing. The surface applier
   % returns one entry per cell with only the top entry nonzero, so today
   % this equals the magnitude of the net; the per-cell form stays so the
   % ledger integrates one shape if a future path fills more entries.
   % Interior transport never reaches this channel: it runs after this
   % budget closes and records its shortfall in the redistribution's own
   % unapplied pair.
   unapplied_vapor_gross = sum(abs(unapplied_weighted));

   % Compute condensation overflow in mwe.
   condensation_overflow = d_rof * dz(1);

   % Add into the forcing-step ledger rather than assigning, because a
   % forcing step can take several accepted substeps and each contributes its
   % own vapor exchange.
   ledger.mass_budget_vapor_solid_mwe = ...
      ledger.mass_budget_vapor_solid_mwe + vapor_solid;
   ledger.mass_budget_vapor_liquid_mwe = ...
      ledger.mass_budget_vapor_liquid_mwe + vapor_liquid;
   ledger.mass_budget_vapor_potential_j_m2 = ...
      ledger.mass_budget_vapor_potential_j_m2 + vapor_potential;
   ledger.mass_budget_unapplied_vapor_j_m2 = ...
      ledger.mass_budget_unapplied_vapor_j_m2 + unapplied_vapor;

   % d_pevp and d_sbl_err are per-substep. d_rof is reset once per forcing
   % step in newtimestep and accumulated across substeps, so it already holds
   % the step total. ASSIGN it, do not add: adding would count the same
   % overflow once per substep.
   ledger.mass_budget_condensation_overflow_mwe = condensation_overflow;

   % Retain absolute substep values so sublimation and deposition within one
   % forcing step cannot cancel to a signed zero.
   ledger.mass_budget_vapor_solid_gross_mwe = ...
      ledger.mass_budget_vapor_solid_gross_mwe + abs(vapor_solid);
   ledger.mass_budget_vapor_liquid_gross_mwe = ...
      ledger.mass_budget_vapor_liquid_gross_mwe + abs(vapor_liquid);
   ledger.mass_budget_vapor_potential_gross_j_m2 = ...
      ledger.mass_budget_vapor_potential_gross_j_m2 + abs(vapor_potential);

   % Overflow only ever grows, so the running total's magnitude already equals
   % the sum of the per-substep magnitudes.
   ledger.mass_budget_condensation_overflow_gross_mwe = ...
      abs(condensation_overflow);
   ledger.mass_budget_unapplied_vapor_gross_j_m2 = ...
      ledger.mass_budget_unapplied_vapor_gross_j_m2 + unapplied_vapor_gross;
end
