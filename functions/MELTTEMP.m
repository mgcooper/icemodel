function Tsfc = MELTTEMP(Tsfc, Tf)
   %MELTTEMP Set surface temperature to melt temp

   %if (depth_swe>0.0 && Tsfc>Tf)
   if Tsfc > Tf
      Tsfc = Tf;
   end
end
