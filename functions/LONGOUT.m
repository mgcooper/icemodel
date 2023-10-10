function Qle = LONGOUT(Tsfc,emiss,SB)
%LONGOUT compute outgoing longwave emitted by the surface
Qle = - emiss*SB*Tsfc^4; % [W m-2]