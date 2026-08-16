function ledger = accumulate_redistribution_budget( ...
      ledger, solid_r, liquid_r, T, f_ice, f_liq, dz)
   %ACCUMULATE_REDISTRIBUTION_BUDGET Record interior vapor transport energy.
   %
   %  ledger = icemodel.column.accumulate_redistribution_budget( ...
   %     ledger, solid_r, liquid_r, T, f_ice, f_liq, dz)
   %
   % Call this once per accepted substep, right after
   % icemodel.column.apply_vapor_transport moves vapor between
   % cells. SOLID_R and LIQUID_R are the storage totals from just before that
   % step.
   %
   % Interior transport conserves mass: it moves vapor between cells and
   % creates none. It does not conserve energy-weighted storage. A cell that
   % gives up vapor from ice releases it at Ls. A cell that takes it into
   % liquid stores it at Lv. The column's Ls-and-Lv-weighted storage therefore
   % moves even though its mass does not. Nothing supplies that energy as a
   % potential, because the enthalpy solve already carries the vapor latent
   % heat through its conduction term.
   %
   % The vapor closure identity compares the surface potential against the
   % energy-weighted storage change, so it needs this term subtracted:
   %
   %   vapor_potential = ro_liq * (Ls * vapor_solid + Lv * vapor_liquid
   %                     + Lv * condensation_overflow)
   %                     + unapplied_vapor - vapor_redistribution
   %
   % The channel stays zero unless the coupled vapor mode is on, so the
   % identity is unchanged for a default run.
   %
   % Inputs
   %   ledger            - Forcing-step ledger.
   %   solid_r, liquid_r - Storage totals [m w.e.] from before the
   %                       redistribution step.
   %   T, f_ice, f_liq   - Column state after the redistribution step.
   %   dz                - Control-volume thickness [m].
   %
   % Output
   %   ledger            - Ledger with the redistribution energy accumulated.
   %
   % See also: icemodel.column.apply_vapor_transport,
   %  icemodel.column.accumulate_vapor_budget
   %
   %#codegen

   persistent Ls Lv ro_liq
   if isempty(Ls)
      [Ls, Lv, ro_liq] = icemodel.physicalConstant('Ls', 'Lv', 'ro_liq');
   end

   [solid_after, liquid_after] = ...
      icemodel.column.integrate_column_budget(T, f_ice, f_liq, dz);

   redistribution = ro_liq * (Ls * (solid_after - solid_r) ...
      + Lv * (liquid_after - liquid_r));

   % Add rather than assign: a forcing step can take several accepted
   % substeps and each one redistributes.
   ledger.mass_budget_vapor_redistribution_j_m2 = ...
      ledger.mass_budget_vapor_redistribution_j_m2 + redistribution;

   % Keep the absolute total so exchanges of opposite sign within one forcing
   % step cannot cancel to a signed zero.
   ledger.mass_budget_vapor_redistribution_gross_j_m2 = ...
      ledger.mass_budget_vapor_redistribution_gross_j_m2 ...
      + abs(redistribution);
end
