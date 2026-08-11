function increments = ablationLedgerIncrements(ledger)
   %ABLATIONLEDGERINCREMENTS Per-interval ablation terms from the mass ledger.
   %
   %  increments = ...
   %     icemodel.verification.helpers.ablationLedgerIncrements(ledger)
   %
   % The comparator scores these terms and the runner plots them. Each caller
   % applies its own prefix and rebasing convention.
   %
   % solid_balance is signed: it falls when refreezing exceeds melt, so it must
   % not be read as geometric surface rise.
   %
   % surface_loss is the mass a top-cell removal exported. A merge gives the
   % joined cell the mean of the pair, so this over-counts what the removed
   % cell held and is not a surface mass flux.
   %
   % solid_vapor_loss is sublimation minus deposition. Runoff does not include
   % it, because sublimated ice never becomes liquid. Adding it to runoff
   % therefore gives a continuous ablation proxy. This term excludes liquid
   % vapor exchange: evaporation removes pore water that runoff already counted
   % as melt that left, and condensation adds pore liquid, not ice.
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
