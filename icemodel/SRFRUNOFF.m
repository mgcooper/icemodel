function ice1 = SRFRUNOFF(ice1,ro_liq,Ls,Lf,dt)
   %SRFRUNOFF compute surface runoff
   %
   %#codegen

   % init arrays
   mlt = zeros(size(ice1.Qe));
   frz = zeros(size(ice1.Qe));
   sub = zeros(size(ice1.Qe));
   con = zeros(size(ice1.Qe));

   % timesteps with surface sublimation / deposition and melt / freeze
   isub = ice1.Qe < 0;
   icon = ice1.Qe > 0;
   imlt = ice1.Qm > 0;
   ifrz = ice1.Qm < 0;

   sub(isub) = ice1.Qe(isub) ./ (Ls .* ro_liq) .* dt .* -1.0;
   con(icon) = ice1.Qe(icon) ./ (Ls .* ro_liq) .* dt;
   mlt(imlt) = ice1.Qm(imlt) ./ (Lf .* ro_liq) .* dt;
   frz(ifrz) = ice1.Qm(ifrz) ./ (Lf .* ro_liq) .* dt .* -1.0;

   rof = mlt + con;

   ice1.melt = cumsum(mlt);
   ice1.freeze = cumsum(frz);
   ice1.subl = cumsum(sub);
   ice1.cond = cumsum(con);
   ice1.runoff = cumsum(rof);
end

%{
function rof = SURFRUNOFF(ice1,opts)

   % 14 Oct 2023, This is an old version SURFRUNOFF (not SRFRUNOFF) when
   melt/freeze/cond were computed in (I think) ICEABLATION and added to ice1 on
   each timestep. This could be revived to experiment with a simple
   surface-based parameterization of refreezing.
   mlt = ice1.melt;
   frz = ice1.freeze;
   con = ice1.cond;
   rof = zeros(size(mlt));

   for n = 1+opts.tlag:length(mlt)
      mltsumlag = sum(mlt(n-opts.tlag:n) + con(n-opts.tlag:n));
      potrunoff = con(n) + mlt(n);
      potfreeze = min(frz(n),mltsumlag);
      potfreeze = max(potfreeze,0.0);

      if isfield('opts','skinfreeze') && opts.skinfreeze == true
         if mltsumlag > 0.0
            netrunoff = potrunoff-potfreeze;
            rof(n,1) = max(rof(n-1)+netrunoff,0.0);
         else
            rof(n,1) = rof(n-1) + potrunoff;
         end
      else
         rof(n,1) = rof(n-1) + potrunoff;
      end
   end
end
%}
