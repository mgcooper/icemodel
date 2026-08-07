function [ledger, vapor_solid, vapor_liquid] = accumulate_vapor_budget( ...
      ledger, solid_p, liquid_p, T, f_ice, f_liq, dz, d_pevp, d_rof, d_sbl_err)
   %ACCUMULATE_VAPOR_BUDGET Add one substep's surface vapor-exchange channels.
   %
   %  [ledger, vapor_solid, vapor_liquid] = ...
   %     icemodel.column.accumulate_vapor_budget( ...
   %     ledger, solid_p, liquid_p, T, f_ice, f_liq, dz, ...
   %     d_pevp, d_rof, d_sbl_err)
   %
   % Call this once per accepted substep, immediately after
   % icemodel.column.budget_surface_mass_balance applies the vapor exchange, so
   % the realized storage change is separated from the potential input and from
   % the two channels the top control volume could not accept.
   %
   % Inputs
   %   ledger              - Forcing-step ledger carrying the phase channels.
   %   solid_p, liquid_p   - Post-phase storage [m w.e.] from
   %                         icemodel.column.accumulate_phase_budget.
   %   T, f_ice, f_liq     - Column state after the vapor exchange.
   %   dz                  - Control-volume thickness [m].
   %   d_pevp              - Potential vapor-driven top-layer liquid change [-].
   %   d_rof               - Condensation overflow in liquid fraction [-].
   %   d_sbl_err           - Signed unapplied vapor-driven ice change [-].
   %
   % Outputs
   %   ledger              - Ledger with the vapor channels accumulated.
   %   vapor_solid,        - Signed substep vapor increments [m w.e.], needed to
   %   vapor_liquid          close storage throughput after remeshing.
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

   % Realized exchange is the storage change across the vapor update alone.
   [solid_v, liquid_v] = icemodel.column.integrate_column_budget(T, f_ice, f_liq, dz);
   vapor_solid = solid_v - solid_p;
   vapor_liquid = liquid_v - liquid_p;

   % The potential input and the rejected channels are energies, not storage:
   % d_pevp and d_sbl_err are top-layer fractions scaled by their latent heats,
   % and the condensation overflow is a liquid depth the column did not store.
   vapor_potential = ro_liq * Lv * d_pevp * dz(1);
   unapplied_vapor = ro_ice * Ls * d_sbl_err * dz(1);

   % d_rof is not like its neighbours. d_pevp and d_sbl_err are per-substep
   % quantities, but d_rof is reset once per forcing step in newtimestep and
   % then accumulated across substeps, so by the time it reaches here it is
   % already the running step total. Adding it on every accepted substep would
   % count the same overflow once per substep. The ledger is also per forcing
   % step, so the current total IS the ledger value.
   condensation_overflow = d_rof * dz(1);

   ledger.mass_budget_vapor_solid_mwe = ...
      ledger.mass_budget_vapor_solid_mwe + vapor_solid;
   ledger.mass_budget_vapor_liquid_mwe = ...
      ledger.mass_budget_vapor_liquid_mwe + vapor_liquid;
   ledger.mass_budget_vapor_potential_j_m2 = ...
      ledger.mass_budget_vapor_potential_j_m2 + vapor_potential;
   ledger.mass_budget_condensation_overflow_mwe = condensation_overflow;
   ledger.mass_budget_unapplied_vapor_j_m2 = ...
      ledger.mass_budget_unapplied_vapor_j_m2 + unapplied_vapor;

   % Retain absolute substep activity so sublimation and deposition within one
   % forcing step cannot cancel to a signed zero.
   ledger.mass_budget_vapor_solid_throughput_mwe = ...
      ledger.mass_budget_vapor_solid_throughput_mwe + abs(vapor_solid);
   ledger.mass_budget_vapor_liquid_throughput_mwe = ...
      ledger.mass_budget_vapor_liquid_throughput_mwe + abs(vapor_liquid);
   ledger.mass_budget_vapor_potential_throughput_j_m2 = ...
      ledger.mass_budget_vapor_potential_throughput_j_m2 + abs(vapor_potential);
   % Overflow only ever grows, so the running total's magnitude already equals
   % the sum of the per-substep magnitudes.
   ledger.mass_budget_condensation_overflow_throughput_mwe = ...
      abs(condensation_overflow);
   ledger.mass_budget_unapplied_vapor_throughput_j_m2 = ...
      ledger.mass_budget_unapplied_vapor_throughput_j_m2 + abs(unapplied_vapor);
end
