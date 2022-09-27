function liqflux = LIQFLUX(frac_liq,frac_ice,dt)
    
% compute the pore space fraction and residual water fraction
    liqresid    =   0.05;               % residual water 
    frac_por    =   1.0 - frac_ice;     % pore fraction
    frac_res    =   liqresid*frac_por;  % residual water fraction

% compute the flux
    if frac_liq > frac_res

% compute the available capacity
        availCap    =   frac_por - frac_res;

% compute the saturated hydraulic conductivity (Colbeck 1972)
        k_snow      =   3.41875e-5 * exp(15.9*(1-frac_ice));

% if grain size is known, try Shimizu 1970 (ro_snow = mLayerTheta?)
%       k_snow      =   7.7e-4*grain_sz^2*9.8/1.787*exp(-0.0078*mLayerTheta*1000);

% compute the relative saturation
        if (availCap > 0.0)
            relSat  =   (frac_liq-frac_res)/availCap;
        else
            relSat  =   1.0;
        end

% calculate the flux out of the snowpack or layer
        mw          =   3.0;
        liqflux     =   dt*k_snow*relSat^mw;
        
% if grain size is known, use the van Genuchten equation
%         n         =   (15.68*exp(-0.46*2*grain_sz))+1;
%         m         =   1 - (1/n);
%         liqflux   =   dt*k_snow*relSat^0.5*((1-((1-relSat^(1/m))^m))^2);
    else
% flow does not occur
        liqflux         =   0.0;
    end

% liq_res     =  0.02;                % resid. water volume / pore volume [-]
% frac_por    =  1 - frac_ice;        % pore fraction                     [-]
% frac_res    =  liq_res * frac_por;  % residual water fraction           [-]
% liq_sat     =  frac_liq / frac_por; % water saturation     [-]
% ava_cap     =  frac_por - frac_res; % available capacity   [-]
% rel_sat     = (frac_liq - frac_res) / ava_cap  % relative saturation  [-]

% clarity on the analogues with colbeck:
% Swi   = liq_res  = frac_res/frac_por = h_res/h_por = h_res/(h_tot-h_ice)
% Sw    = liq_sat  = frac_liq/frac_por = h_liq/h_por = h_liq/(h_tot-h_ice)
% phi   = porosity = frac_por          = h_por/h_tot = (h_tot-h_ice)/h_tot
% Sstar = rel_sat  = (frac_liq-frac_res)/ava_cap = (h_liq-h_res)/ava_cap



% NOTE: didn't clarify this
% %--------------------------------------------------------------------------
% function [  fracliquid,                                                 ...
%             mLayerVolFracLiq,                                           ...
%             mLayerVolFracIce,                                           ...
%             retent ]            =   UPDATE_STATE(Tk,mLayerTheta,        ...
%                                     fc_param,Tf,ro_water,retent)
% %--------------------------------------------------------------------------
%     
%     iden_ice            =   917;
% 
% % fraction of liquid water
%     fracliquid          =   1.0/(1.0 + (fc_param * (Tf - min(Tk,Tf)))^2.0);
% 
% % calculate volumetric fractions of liquid and ice content
%     mLayerVolFracLiq    =   fracliquid * mLayerTheta;
%     mLayerVolFracIce    =   (1.0 - fracliquid) * mLayerTheta * (ro_water/iden_ice);
% 
% % calculate the retention
%     retent              =   min(0.02,max(retent,0.75*mLayerVolFracLiq));
% 
% end


