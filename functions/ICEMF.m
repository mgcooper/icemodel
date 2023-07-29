function [  T,                                                          ...
            f_ice,                                                      ...
            f_liq,                                                      ...
            d_liq,                                                      ...
            d_drn,                                                      ...
            x_err,                                                      ...
            lcflag ] =  ICEMF(T,f_ice,f_liq,ro_ice,ro_liq,cp_ice,       ...
                        Lf,Ls,Lv,Tf,TL,fcp,xf_liq,Sc,Sp,JJ_therm,       ...
                        f_min,fopts,dz_therm,dt_new,Qe,liqflag,         ...
                        ro_iwe,d_liq,d_drn,flmin,metiter)
%ICEMF compute ice melt-freeze and combine layers if necessary

% Compute d_liq first. xf_liq is the value of f_liq returned by this function at
% the end of the last step. The only process that affects f_liq between its
% assignment to xf_liq and here is melt/freeze in ICEENBAL. Therefore, d_liq
% represents the change in f_liq due to melt/freeze. Evaporaton/condensation is
% accounted for next, and at present is not d_drn gets updated  
% if layers are combined   
d_liq  = d_liq + f_liq - xf_liq; % d_liq < 0 means refreezing, d_liq > 0 melt
xf_liq = f_liq;

% want: 
% d_liq > 0 = melt, done
% d_evp > 0 = cond, done
% d_drn > 0 = runoff due to extra condensation, done
% d_com > 0 = runoff in layer combo ... 
% 
% d_com is tricky ... but recall its always the top layer that is removed, so
% if f_liq - xf_liq < 0, it means the water content is lower after layer combo,
% so some water was expelled as runoff, but in the melt/freeze case, 
% f_liq - xf_liq < 0 means refreezing
% if f_liq - xf_liq > 0, it means layer combination increased the water content,
% ... pick up on this 

% update the lower melt-zone boundary
f_liq_min = (f_liq(1)+f_ice(1).*ro_iwe).*flmin; % 0.001;

% evaporation/condensation/sublimation
[  f_ice, ...
   f_liq, ...
   d_drn, ...
   x_err ]  =  ICESUBL(f_ice, f_liq, d_drn, ro_ice, ro_liq, Ls, Lv, ...
               f_min, f_liq_min, dz_therm, dt_new, Qe, liqflag, metiter);

% need:
d_evp = f_liq - xf_liq;
% if evap, d_evp < 0
% if cond, d_evp > 0

% % budget evap / subl
% d_evp = f_liq - xf_liq;
% d_sub = f_ice - xf_ice;

% update the top layer temperature
T(1) = Tf-sqrt(((f_liq(1)+f_ice(1).*ro_iwe)./f_liq(1)-1))./fcp;

% combine layers if any layer is <f_min, or if this step's sublimation
% would reduce any layer to <f_min (predict the need to combine next step)
lyrmrg = f_ice <= f_min | (f_ice+Qe/(Ls*ro_ice)*dt_new/dz_therm) <= f_min;
lcflag = false(size(f_ice));
 
% if lyrmerge is all false and no layers are < h_min, return (no combine)
if any(lyrmrg) % && ~any(h_ice<=h_min)

% This is a counter that keeps the ji index on track with the original 1:JJ
%   column, since layers are removed within the loop
   ii = 0;
    
   for j = 1:JJ_therm
      
      ji = j + ii;

      % Combine layers
      if lyrmrg(ji)==true
         
         [j1,j2] = LAYERINDS(ji,f_ice);

% prep for layer combination
         lcflag(j1) = true;
    
      [  f_liq(j2),  ...
         f_ice(j2),  ...
         T(j2),      ...
         Sc(j2),     ...
         Sp(j2) ]    = COMBINEHEAT(f_liq,f_ice,T,Sc,Sp,Tf,TL, ...
                       cp_ice,Lf,ro_ice,ro_liq,fcp,j1,j2,fopts);
            
% Remove the combined layers and add new layers to the bottom.
         f_ice       = vertcat(f_ice(~lyrmrg),f_ice(end));
         f_liq       = vertcat(f_liq(~lyrmrg),f_liq(end));
         T           = vertcat(T(~lyrmrg),T(end));
         Sc          = vertcat(Sc(~lyrmrg),Sc(end));
         Sp          = vertcat(Sp(~lyrmrg),Sp(end));    
         lyrmrg      = vertcat(lyrmrg(~lyrmrg),lyrmrg(end));
         ii          = ii - 1;
      else
         %ii         = ii + 0; % for reference, not needed
         %cflag(j)   = false;
      end
   end
   
   % need:
   % d_com = f_liq - xf_liq, but f_liq = xf_liq + d_evp, so:
   % d_com = f_liq - xf_liq - d_evp;

   % this says: drainage = extra condensation + change in f_liq due to
   % condensation and 
   d_drn = d_drn + xf_liq - f_liq + d_evp; % = d_drn - d_com
   
   % think this is wrong
   % d_drn = d_drn + f_liq - xf_liq - d_evp; % = d_drn + d_com
   
   % original:
   % d_drn = d_drn + xf_liq - f_liq;
   
   % test 8
   % d_drn = d_drn + f_liq - xf_liq;
end

% % f_liq budget:
% 1 d_liq = f_liq - xf_liq
% 2 xf_liq = f_liq
% 3 f_liq = f_liq + d_evap (evap / cond in ICESUBL)
% 4 d_drn = d_drn + x_evap (extra condensation that runs off in ICESUBL)
% 5 f_liq = f_liq + d_comb (extra water that runs off in COMBINEHEAT)
% 6 d_drn = d_drn + xf_liq - f_liq
% 
% Combining terms, step 6 can be written:
% d_drn = x_evap + f_liq - (f_liq _ d_evap + d_comb)
% d_drn = x_evap - d_evap - d_comb
% d_drn = x_evap - d_evap - d_comb
% 
% if x_evap is non-zero, it is guaranteed to be positive, and it means
% condensation occurred, therefore d_evap is positive, and d_evap 


% % f_liq budget:
% 1 d_liq = f_liq - xf_liq
% 2 xf_liq = f_liq
% 3 f_liq = f_liq + d_evap (evap / cond in ICESUBL)
% 4 d_drn = d_drn + x_evap (extra condensation that runs off in ICESUBL)
% 5 f_liq = f_liq + d_comb (extra water that runs off in COMBINEHEAT)
% 6 d_drn = d_drn + xf_liq - f_liq
% 
% Combining terms, step 6 can be written:
% d_drn = x_evap + f_liq - (f_liq _ d_evap + d_comb)
% d_drn = x_evap - d_evap - d_comb
% d_drn = x_evap - d_evap - d_comb
% 
% if x_evap is non-zero, it is guaranteed to be positive, and it means
% condensation occurred, therefore d_evap is positive, and d_evap 

% assert(f_liq(1) + f_ice(1) <= 1.0)
% printf(f_liq(1), 10)
% printf(f_ice(1), 10)
% printf(T(1), 10)
% printf(f_liq_min(1), 10)

% xf_liq is from the prior met-step, so each sub-step d_liq is overridden,
% and the total d_liq is relative to each met-step, as it should be,
% whereas the other values in here are updated at each substep. However,
% what happens when a layr is removed within a sub-stepping?