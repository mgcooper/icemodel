%--------------------------------------------------------------------------
%   CALCULATE SURFACE MELT/FREEZE, SUBLIMATION, AND ABLATION RATE
%--------------------------------------------------------------------------

function [  surf_mlt,                                                   ...
            surf_frz,                                                   ...
            surf_sub,                                                   ...
            isubl,                                                      ...
            f_ice    ] =   ICEABLATION2(Qm,Qe,Qf,surf_mlt,surf_frz,     ...
                           surf_sub,ro_ice,ro_liq,Lf,Ls,dz,dt,opts,f_ice)
%--------------------------------------------------------------------------
% All the liquid fluxes are added to surf_runoff. Each individual flux is
% saved as a cumulative sum. Runoff is calculated as a post-process step
% (see SURF_RUNOFF), which accounts for condensation that would have
% refrozen and melt available for runoff
   
   imeltfreeze =  0.0;
   
% negative = surface sublimation, positive = deposition
   isubl       =  Qe / (Ls * ro_liq) * dt;                   % [m w.e.]
   surf_sub    =  surf_sub + isubl;
      
% reduce/increase the ice content due to sublimation in the top layer
   f_ice(1)    =  f_ice(1) + isubl*ro_liq/ro_ice/dz(1);
   
% surface melt / freeze
if (Qm>0.0 && opts.skinmodel==true) || (Qm>0.0 && opts.skinmelt==true)
      
   imeltfreeze =  Qm / (Lf * ro_liq) * dt;
   surf_mlt    =  surf_mlt + imeltfreeze;

elseif (Qf > 0.0 && opts.skinfreeze==true)
   
   imeltfreeze =  Qf / (Lf * ro_liq) * dt;
   surf_frz    =  surf_frz + imeltfreeze;
end

% decrease the ice surface and add imelt to surf_runoff
   f_ice(1)    =  f_ice(1) + imeltfreeze*ro_liq/ro_ice/dz(1);

% this would instead add the melt to the upper grid cell
   % f_liq(1) =  f_liq(1) + imelt/dz;

end

% note: 
% isubl  = -Qe/(Ls*row)*dt [m w.e.]
% hsubl  = isubl*row/roi   [m i.e.]
% fsubl  = hsubl/ht = isubl*row/roi/ht; (see f_ice update in code above)

