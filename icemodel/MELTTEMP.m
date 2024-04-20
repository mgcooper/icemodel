function Ts = MELTTEMP(Ts, Tf)
   %MELTTEMP Set surface temperature to melt temperature
   %
   %  Ts = MELTTEMP(Ts, Tf)
   %
   % See also: MFENERGY, ENBALANCE

   % if (depth_swe > 0.0 && Ts > Tf)
   if Ts > Tf
      Ts = Tf;
   end
end
