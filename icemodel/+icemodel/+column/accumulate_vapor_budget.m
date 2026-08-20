function budget = accumulate_vapor_budget( ...
      budget, T, f_ice, f_liq, dz, d_pevp, d_rof, d_sbl_err)
   %ACCUMULATE_VAPOR_BUDGET Add one substep's surface vapor-exchange increments.
   %
   %  budget = icemodel.column.accumulate_vapor_budget( ...
   %     budget, T, f_ice, f_liq, dz, d_pevp, d_rof, d_sbl_err)
   %
   % icemodel.column.budget_surface_mass_balance calls this once per
   % accepted substep, right after it applies the vapor exchange.
   % This keeps the storage change apart from the potential input and from the
   % two amounts the top control volume cannot accept. Those two amounts are
   % condensation above the pore capacity, and deposition the control volume
   % has no energy to apply.
   %
   % The pre-exchange baseline comes from budget.substep (solid_p,
   % liquid_p, written by accumulate_phase_budget); the post-exchange
   % storage written back to budget.substep (solid_v, liquid_v) is the
   % baseline the interior-transport budget consumes.
   %
   % Inputs
   %   budget              - Forcing-step budget with the phase channels and
   %                         post-phase baselines already accumulated.
   %   T, f_ice, f_liq     - Column state after the vapor exchange.
   %   dz                  - Control-volume thickness [m].
   %   d_pevp              - Potential vapor-driven top-layer liquid change [-].
   %   d_rof               - Condensation overflow in liquid fraction [-].
   %   d_sbl_err           - Signed unapplied vapor-driven ice change per
   %                         cell [-].
   %
   % Outputs
   %   budget              - Budget with the accumulated vapor increments and
   %                         the post-exchange baselines in budget.substep.
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
   vapor_solid = solid_v - budget.substep.solid_p;
   vapor_liquid = liquid_v - budget.substep.liquid_p;

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
   % cell's rejection over the whole column and inflate the channel from
   % the dz(1) term to a sum(dz) one. Weight the scalar form by dz(1) and
   % keep the per-cell sum for the vector form.
   if isscalar(d_sbl_err)
      unapplied_vapor = ro_ice * Ls * d_sbl_err * dz(1);
   else
      unapplied_vapor = sum(ro_ice * Ls * d_sbl_err(:) .* dz(:));
   end

   % Add into the forcing-step budget rather than assigning, because a
   % forcing step can take several accepted substeps and each contributes its
   % own vapor exchange.
   budget.mass_budget_vapor_solid_mwe = ...
      budget.mass_budget_vapor_solid_mwe + vapor_solid;
   budget.mass_budget_vapor_liquid_mwe = ...
      budget.mass_budget_vapor_liquid_mwe + vapor_liquid;
   budget.mass_budget_vapor_potential_j_m2 = ...
      budget.mass_budget_vapor_potential_j_m2 + vapor_potential;
   budget.mass_budget_unapplied_vapor_j_m2 = ...
      budget.mass_budget_unapplied_vapor_j_m2 + unapplied_vapor;

   % d_pevp and d_sbl_err are per-substep. d_rof is reset once per forcing
   % step in newtimestep and accumulated across substeps, so it already holds
   % the step total. ASSIGN it, do not add: adding would count the same
   % overflow once per substep.
   budget.mass_budget_condensation_overflow_mwe = d_rof * dz(1);

   % Hand the post-exchange baselines to the interior-transport budget.
   budget.substep.solid_v = solid_v;
   budget.substep.liquid_v = liquid_v;
end
