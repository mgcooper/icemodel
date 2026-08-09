function increments = ablationLedgerIncrements(ledger)
   %ABLATIONLEDGERINCREMENTS Per-interval ablation terms from the mass ledger.
   %
   %  increments = ...
   %     icemodel.verification.helpers.ablationLedgerIncrements(ledger)
   %
   % The comparator scores these and the runner plots them. Defining them once
   % keeps the scored quantity and the plotted curve the same quantity; each
   % caller still applies its own prefix and rebasing convention.
   %
   % solid_balance is signed: it falls when refreezing exceeds melt, so it must
   % not be read as geometric surface rise.
   %
   % surface_loss is the mass a top-cell removal exported. A merge gives the
   % joined cell the mean of the pair, so this over-counts what the removed
   % cell held and is not a surface mass flux.
   %
   % solid_vapor_loss is sublimation minus deposition. Runoff cannot see it,
   % because sublimated ice never becomes liquid, so adding it to runoff gives
   % a continuous ablation proxy. Liquid vapor exchange stays out: evaporation
   % removes pore water that runoff already counted as melt that left, and
   % condensation adds pore liquid rather than ice.
   %
   % Inputs
   %  ledger - model rows carrying the mass_budget_* channels
   %
   % Outputs
   %  increments - struct with solid_balance, surface_loss, solid_vapor_loss

   increments = struct( ...
      'solid_balance', -(ledger.mass_budget_phase_solid_mwe ...
      + ledger.mass_budget_vapor_solid_mwe), ...
      'surface_loss', ledger.mass_budget_top_export_solid_mwe ...
      + ledger.mass_budget_top_export_liquid_mwe, ...
      'solid_vapor_loss', -ledger.mass_budget_vapor_solid_mwe);
end
