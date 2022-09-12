%--------------------------------------------------------------------------
%   COMPUTE THE SENSIBLE HEAT FLUX
%--------------------------------------------------------------------------

function Qh = SENSIBLE(De,S,Tair,Tsfc,cv_air)
%--------------------------------------------------------------------------

% for speed:
    Qh      =   cv_air * De * S * (Tair - Tsfc);
    
% for clarity:
%     Qh      =   ro_air * Cp_air * De_h * stability * (Tair - Tsfc);     
%   [W m-2] =   [kg m-3] * [J kg-1 K-1] * [m s-1] * [-] * [K]  	  

% to implement scalar coeffs hoffman uses D_h
%   Qh = ro_air * Cp * D_h * stability * (Tair - Tsfc)   