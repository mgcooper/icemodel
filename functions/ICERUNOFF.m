function ice1 = ICERUNOFF(ice1,ice2,opts)
   %ICERUNOFF compute runoff from the ice column

   dz = opts.dz_thermal;
   df_liq = ice2.df_liq;

   % partition runoff into melt/freeze
   melt = zeros(size(df_liq,2), 1);
   frze = zeros(size(df_liq,2), 1);
   for n = 1:size(df_liq,2)
      melt(n) = sum(dz(1).*df_liq(df_liq(:,n)>0, n));
      frze(n) = sum(-dz(1).*df_liq(df_liq(:,n)<0, n));
   end

   tlag = opts.tlagcolumn;
   runoff = zeros(size(melt));

   for n = 1+tlag:length(melt)
      meltsumlag = sum(melt(n-tlag:n));
      potrunoff = melt(n);
      potfreeze = min(frze(n),meltsumlag);
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
   ice1.freeze = cumsum(frze);         % cumulative freeze

   % compute cumulative drainage if it's included in the output
   if isfield(ice2,'df_drn')
      ice1.drain = transpose(cumsum(sum(dz(:) .* ice2.df_drn)));
   end

   % compute runoff and drain directly from df_liq (note: df_liq is the change
   % in liquid water fraction, which is the change in liquid water thickness
   % divided by the layer thickness dz, so multiplying by dz here cancels dz
   % and you get the change in liquid water thickness of each layer,
   % regardless of whatever density is/was for each df_liq)
end
