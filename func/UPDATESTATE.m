function [  k_eff,                                                      ...
            ro_sno,                                                     ...
            cp_sno,                                                     ...
            drovdT ] =  UPDATESTATE(T,f_ice,f_liq,ro_ice,ro_liq,ro_air, ...
                        cp_ice,cp_liq,k_liq,Ls,Rv,Tf)

% consituent bulk densities                                       kg/m3
   bd_ice      =  f_ice.*ro_ice;
   bd_liq      =  f_liq.*ro_liq;
   bd_air      =  (1.0-f_liq-f_ice).*ro_air;
  
% combined bulk density  (pt)
   ro_sno      =  (bd_ice + bd_liq + bd_air);         % (eq 3)    kg/m3
   
% constituent bulk volume-specific heat capacities                J/m3/K
   cv_ice      =  cp_ice.*bd_ice;
   cv_liq      =  cp_liq.*bd_liq;
 % cv_air      =  cp_air.*bd_air;
    
% combined bulk mass-specific heat capacity (ct)
   cp_sno      =  (cv_ice + cv_liq) ./ ro_sno;        % (eq 55)   J/kg/K
    
% effective thermal conductivity    
   k_eff       =  GETGAMMA(T,f_liq,f_ice,ro_ice,k_liq,Ls,Rv,Tf);

% change in saturation vapor density with temperature             kg/m3
%	[~,drovdT]  =  VAPORHEAT(T,f_liq,f_ice,Tf,Rv,Ls);
    
% energy coefficient                                              W/m2/K
%   aP0         =   bd_sno.*cp_sno.*dz./dt;