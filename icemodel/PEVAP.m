function [d_pevp, pevp] = PEVAP(Qe, Lv, ro_liq, dt, dz)
   %PEVAP Potential vapor tendency as top-layer liquid fraction change and flux
   %
   % Definition:
   % d_pevp = (Qe / (Lv * ro_liq)) * (dt / dz)
   %        = pevp * (dt / dz)
   % pevp   = Qe / (Lv * ro_liq)
   %
   % This uses latent heat flux [W m-2] to:
   % 1) convert heat flux to a water-equivalent mass flux [kg m-2 s-1]
   % 2) divide by ro_liq [kg m-3] to get pevp [m s-1]
   % 3) integrate pevp over dt [s] and normalize by dz [m] to get d_pevp [-]
   %
   % Units:
   % [m s-1] = [J s-1 m-2] / ([J kg-1] * [kg m-3])
   % [-] = [J s-1 m-2] / ([J kg-1] * [kg m-3]) * [s] / [m]
   %
   % pevp > 0 and d_pevp > 0 : condensation potential
   % pevp < 0 and d_pevp < 0 : evaporation/sublimation potential
   %
   % Note: To budget sublimation, outputs require conversion to ice frac change:
   % d_psbl = d_pevp * (Lv * ro_liq) / (Ls * ro_ice)
   %
   % See also: ICESUBL
   %
   %#codegen
   pevp = Qe / (Lv * ro_liq);
   d_pevp = pevp * dt / dz;
end
