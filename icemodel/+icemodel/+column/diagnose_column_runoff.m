function ice1 = diagnose_column_runoff(ice1, ice2, opts)
   %DIAGNOSE_COLUMN_RUNOFF Diagnose cumulative runoff from column mass changes.
   %
   % Runoff is a postprocessed water budget for the column reservoir, not a
   % direct modeled flux (the column never drains during runtime). Per step the
   % reservoir gains melt and condensation and loses refreezing, evaporation,
   % and runoff. Condensation above the top cell's pore capacity never
   % enters the reservoir, so it runs off here with no residence time.
   %
   % This function credits refreezing only up to the liquid supplied within a
   % trailing opts.tlag window. That residence time stops the column's whole
   % accumulated meltwater from refreezing when the melt season ends. It also
   % limits how much of a day's melt can refreeze overnight.
   %
   % Runoff is a running total, not a cumulative sum. An hour whose evaporation
   % exceeds its melt and overflow lowers the total.
   %
   %#codegen

   dz = opts.dz_thermal;
   df_liq = ice2.df_liq;
   n_steps = size(df_liq, 2);

   % Partition melt/freeze phase change only. In budget_surface_mass_balance,
   % df_liq is assigned before vapor exchange, so melt and refreezing here do
   % not include evaporation or condensation.
   melt = zeros(n_steps, 1);
   freeze = zeros(n_steps, 1);
   for n = 1:n_steps
      melt(n) = sum(dz(1) .* df_liq(df_liq(:, n) > 0, n));
      freeze(n) = sum(-dz(1) .* df_liq(df_liq(:, n) < 0, n));
   end

   % Liquid vapor exchange, signed positive into the reservoir. In
   % budget_surface_mass_balance, df_evp differences f_liq (after melt/freeze
   % phase change), so it represents condensation and evaporation only;
   % sublimation and deposition act directly on f_ice at each substep update.
   vapor_liquid = liquidVaporSupply(ice2, dz, n_steps);

   % Overflow is condensation the top cell's pore capacity rejected.
   overflow = overflowSupply(ice1, dz, n_steps);

   % Supply is liquid that arrives and can later refreeze. Evaporation removes
   % liquid that would otherwise have run off, so it reduces runoff rather than
   % the liquid available for refreezing.
   supply = melt + max(vapor_liquid, 0.0);

   runoff = zeros(n_steps, 1);
   for n = 1 + opts.tlag:n_steps
      supplysum = sum(supply(n - opts.tlag:n));
      potfreeze = max(min(freeze(n), supplysum), 0.0);
      netrunoff = melt(n) + vapor_liquid(n) + overflow(n) - potfreeze;
      runoff(n, 1) = max(runoff(n - 1) + netrunoff, 0.0);
   end

   ice1.melt = cumsum(melt);           % cumulative melt
   ice1.runoff = runoff;               % cumulative runoff
   ice1.freeze = cumsum(freeze);       % cumulative freeze

   % Cumulative mass removed by remeshing, in mwe. df_lyr is the
   % water-equivalent fraction each merge exported, so scaling by the cell
   % thickness gives mass. This is not surface-displacement: only removals of
   % the top cell translate the grid downward.
   if isfield(ice2, 'df_lyr')
      ice1.dlayer = transpose(cumsum(sum(dz(:) .* ice2.df_lyr)));
   end
end

function values = liquidVaporSupply(ice2, dz, n_steps)
   %LIQUIDVAPORSUPPLY Return signed liquid vapor exchange per step [m].

   % Initialize an empty column.
   values = zeros(n_steps, 1);

   % Every icemodel output profile includes df_evp. Reduced test payloads can
   % omit it, so return a zero column when the field is missing.
   if ~isfield(ice2, 'df_evp')
      return
   end

   % Otherwise integrate the liquid vapor exchange over the column in mwe.
   values = reshape(sum(dz(1) .* ice2.df_evp, 1), [], 1);
end

function values = overflowSupply(ice1, dz, n_steps)
   %OVERFLOWSUPPLY Return condensation overflow per step [m].
   %
   % ice1.df_rof is a top-cell liquid fraction, matching how df_liq and df_evp
   % are stored, so it is scaled by the same control-volume thickness to compute
   % runoff in mwe.

   % Initialize a column of zeros.
   values = zeros(n_steps, 1);

   % Return a column of zeros to callers who don't have df_rof in their output.
   if ~isfield(ice1, 'df_rof')
      return
   end

   % Otherwise, compute the overflow from the top cell in mwe.
   values = dz(1) .* reshape(double(ice1.df_rof), [], 1);
end
