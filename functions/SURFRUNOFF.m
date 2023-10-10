function rof = SURFRUNOFF(ice1,opts)

   tlag = opts.tlagsurf;
   mlt = ice1.melt;
   frz = ice1.freeze;
   con = ice1.cond;
   rof = zeros(size(mlt));


   for n = 1+tlag:length(mlt)
      meltsumlag = sum(mlt(n-tlag:n) + con(n-tlag:n));
      potrunoff = con(n) + mlt(n);
      potfreeze = min(frz(n),meltsumlag);
      potfreeze = max(potfreeze,0.0);

      if isfield('opts','skinfreeze') && opts.skinfreeze == true
         if meltsumlag > 0.0
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

