function ice1 = SRFRUNOFF(ice1,ro_liq,Ls,Lf,dt)
   %SRFRUNOFF compute surface runoff

   % init arrays
   mlt = zeros(size(ice1.Qe));
   frz = zeros(size(ice1.Qe));
   sub = zeros(size(ice1.Qe));
   con = zeros(size(ice1.Qe));

   % timesteps with surface sublimation / deposition and melt / freeze
   isub = ice1.Qe < 0;
   icon = ice1.Qe > 0;
   imlt = ice1.Qm > 0;
   ifrz = ice1.Qf > 0;

   sub(isub) = -ice1.Qe(isub) ./ (Ls .* ro_liq) .* dt;
   con(icon) = ice1.Qe(icon) ./ (Ls .* ro_liq) .* dt;
   mlt(imlt) = ice1.Qm(imlt) ./ (Lf .* ro_liq) .* dt;
   frz(ifrz) = ice1.Qf(ifrz) ./ (Lf .* ro_liq) .* dt;

   rof = mlt + con;

   ice1.melt = cumsum(mlt);
   ice1.freeze = cumsum(frz);
   ice1.subl = cumsum(sub);
   ice1.cond = cumsum(con);
   ice1.runoff = cumsum(rof);
end
