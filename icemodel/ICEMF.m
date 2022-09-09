%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [  T,                                                          ...
            f_ice,                                                      ...
            f_liq,                                                      ...
            d_liq,                                                      ...
            d_drn,                                                      ...
            lcflag ]    =  ICEMF(T,f_ice,f_liq,ro_ice,ro_liq,cp_ice,    ...
                           Lf,Ls,Lv,Tf,TL,fcp,xf_liq,Sc,Sp,JJ_therm,    ...
                           f_min,fopts,dz_therm,dt_new,Qe,liqflag,      ...
                           ro_iwe,d_liq,d_drn)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% xf_liq is from the prior met-step, so each sub-step d_liq is overridden,
% and the total d_liq is relative to each met-step, as it should be,
% whereas the other values in here are updated at each substep. However,
% what happens when a layr is removed within a sub-stepping?

% need to test putting this after subl 
% d_drn gets updated if layers are combined      
   d_liq       =  d_liq + f_liq - xf_liq; % test
   xf_liq      =  f_liq;
   
% these have dimensions of f_liq and f_ice i.e. h_liq/h and h_ice/h
   if liqflag == true
      f_liq(1) =  min(max(f_liq(1)+Qe/(Lv*ro_liq)*dt_new/dz_therm,0),1);
   else
      f_ice(1) =  min(max(f_ice(1)+Qe/(Ls*ro_ice)*dt_new/dz_therm,f_min),1);
   end

% update the top layer temperature
   T(1)        =  Tf-sqrt(((f_liq(1)+f_ice(1).*ro_iwe)./f_liq(1)-1))./fcp;
   
% combine layers if any layer is <f_min, or if this step's sublimation
% would reduce any layer to <f_min (predict the need to combine next step)
   lyrmrg      =  (f_ice+Qe/(Ls*ro_ice)*dt_new/dz_therm)<=f_min;
   lcflag      =  false(size(f_ice));
    
% if lyrmerge is all false and no layers are < h_min, return (no combine)
if any(lyrmrg) % && ~any(h_ice<=h_min)

% This is a counter that keeps the ji index on track with the original 1:JJ
%   column, since layers are removed within the loop
    ii = 0;
    
    for j = 1:JJ_therm
        
        ji = j + ii;

% Combine layers
        if (f_ice(ji)+Qe/(Ls*ro_ice)*dt_new/dz_therm)<=f_min || lyrmrg(ji)==true
            
            [j1,j2]     =  LAYERINDS(ji,f_ice);

% prep for layer combination
            lcflag(j1)  =  true;
    
        [   f_liq(j2),                                                  ...
            f_ice(j2),                                                  ...
            T(j2),                                                      ...
            Sc(j2),                                                     ...
            Sp(j2) ]    =  COMBINEHEAT(f_liq,f_ice,T,Sc,Sp,Tf,TL,       ...
                           cp_ice,Lf,ro_ice,ro_liq,fcp,j1,j2,fopts);
            
% Remove the combined layers and add new layers to the bottom.
            f_liq(j1)   =  [];
            f_ice(j1)   =  [];
            T(j1)       =  [];
            Sc(j1)      =  [];
            Sp(j1)      =  [];
            lyrmrg(j1)  =  [];
            f_ice       =  vertcat(f_ice,f_ice(end));
            f_liq       =  vertcat(f_liq,f_liq(end));
            T           =  vertcat(T,T(end));
            Sc          =  vertcat(Sc,Sc(end));
            Sp          =  vertcat(Sp,Sp(end));    
            lyrmrg      =  vertcat(lyrmrg,lyrmrg(end));
            ii          =  ii - 1;
        else
            %ii         =   ii + 0; % for reference, obvi not needed
            %cflag(j)   =   false;
        end     
    end
%     d_liq =  d_liq + xf_liq-f_liq; % test
    d_drn = xf_liq-f_liq;
end
