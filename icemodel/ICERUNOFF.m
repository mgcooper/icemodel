function ice1 = ICERUNOFF(ice1,ice2,opts)
   %ICERUNOFF compute runoff from the ice column
   %
   %#codegen

   dz = opts.dz_thermal;
   df_liq = ice2.df_liq;

   % partition runoff into melt/freeze
   melt = zeros(size(df_liq,2), 1);
   freeze = zeros(size(df_liq,2), 1);
   for n = 1:size(df_liq,2)
      melt(n) = sum(dz(1).*df_liq(df_liq(:,n)>0, n));
      freeze(n) = sum(-dz(1).*df_liq(df_liq(:,n)<0, n));
   end

   tlag = opts.tlagcolumn;
   runoff = zeros(size(melt));

   for n = 1+tlag:length(melt)
      meltsumlag = sum(melt(n-tlag:n));
      potrunoff = melt(n);
      potfreeze = min(freeze(n),meltsumlag);
      potfreeze = max(potfreeze,0.0);
      if meltsumlag > 0.0
         netrunoff = potrunoff-potfreeze;
         runoff(n,1) = max(runoff(n-1)+netrunoff, 0.0);
      else
         runoff(n,1) = runoff(n-1) + potrunoff;
      end
   end
   ice1.melt = cumsum(melt);           % cumulative melt
   ice1.runoff = runoff;               % cumulative runoff
   ice1.freeze = cumsum(freeze);       % cumulative freeze

   % compute cumulative layer change if it's included in the output
   if isfield(ice2, 'df_lyr')
      ice1.dlayer = transpose(cumsum(sum(dz(:) .* ice2.df_lyr)));
   end
end
