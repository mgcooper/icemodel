function Ts = MELTTEMP(Ts, Tf)
   %MELTTEMP Set surface temperature to melt temp

   % if (depth_swe > 0.0 && Ts > Tf)
   if Ts > Tf
      Ts = Tf;
   end
end
