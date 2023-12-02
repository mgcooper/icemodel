function Qle = LONGOUT(Ts, emiss, SB)
   %LONGOUT compute outgoing longwave emitted by the surface
   Qle = -emiss * SB * Ts ^ 4; % [W m-2]
end
