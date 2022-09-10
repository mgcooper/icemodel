function H = TOTALHEAT(T,f_liq,f_ice,cv_liq,cv_ice,roLf,Tf,H_vap,dz,dt)
% this is enthalpy (J/m3/K) converted to W/m2
   H  = ((cv_ice.*f_ice+cv_liq.*f_liq).*(T-Tf)+roLf.*f_liq+H_vap).*dz./dt;