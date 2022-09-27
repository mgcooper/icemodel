%--------------------------------------------------------------------------
%   COMPUTE THE SOLAR RADIATION PENETRATING THE ICE SURFACE
%--------------------------------------------------------------------------
function Qsip = SOLARPEN(Qsi,chi)
    Qsip = (1.0 - chi) * Qsi;