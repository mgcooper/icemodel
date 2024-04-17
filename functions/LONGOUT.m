function Qle = LONGOUT(Ts, emiss, SB)
   %LONGOUT compute outgoing longwave emitted by the surface
   %
   %  Qle = LONGOUT(Ts, emiss, SB)
   %
   % See also: LONGIN, ENBALANCE, SEBFLUX

   Qle = -emiss * SB * Ts ^ 4; % [W m-2]
end
