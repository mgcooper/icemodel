function Tsfc0 = MELTTEMP(Tsfc, Tf)
   %MELTTEMP set surface temperature to melt temp

   %if (depth_swe>0.0 && Tsfc>Tf)
   if Tsfc > Tf
      Tsfc0 = Tf;
   else
      Tsfc0 = Tsfc;
   end
end
