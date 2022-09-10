function ice1 = SRFRUNOFF(enbal,ro_liq,Ls,Lf,dt)

% init arrays
   mlt      =  zeros(size(enbal.Qe));
   frz      =  zeros(size(enbal.Qe));
   sub      =  zeros(size(enbal.Qe));
   con      =  zeros(size(enbal.Qe));
   
% timesteps with surface sublimation / deposition and melt / freeze
   isub     =  enbal.Qe < 0;
   icon     =  enbal.Qe > 0;
   imlt     =  enbal.Qm > 0;
   ifrz     =  enbal.Qf > 0;

   sub(isub)   = -enbal.Qe(isub) ./ (Ls .* ro_liq) .* dt;
   con(icon)   =  enbal.Qe(icon) ./ (Ls .* ro_liq) .* dt;
   mlt(imlt)   =  enbal.Qm(imlt) ./ (Lf .* ro_liq) .* dt;
   frz(ifrz)   =  enbal.Qf(ifrz) ./ (Lf .* ro_liq) .* dt;

   rof         =  mlt + con;
   
   ice1.melt   = tocolumn(cumsum(mlt));
   ice1.freeze = tocolumn(cumsum(frz));
   ice1.subl   = tocolumn(cumsum(sub));
   ice1.cond   = tocolumn(cumsum(con));
   ice1.runoff = tocolumn(cumsum(rof));