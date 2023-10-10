function liqflux = LIQFLUX(frac_liq,frac_ice,dt)
   %LIQFLUX Compute the liquid water flux.

   % march 2023, it should be straightforward to compute this after each
   % timestep, for each layer, and whatever gets to the base, either release
   % immediately, or put in a conceptual aquifer that drains slowly. but the
   % main thing is that at each step, this moves water around, so f_liq changes,
   % and total heat changes, and if the timestep is small enough, it could go at
   % the end of each step as an update state type thing, but if that fails, it
   % will need to go inside the solver

   % compute the pore space fraction and residual water fraction
   liqresid = 0.05;               % residual water
   frac_por = 1.0 - frac_ice;     % pore fraction
   frac_res = liqresid*frac_por;  % residual water fraction

   % compute the flux
   if frac_liq > frac_res

      % compute the available capacity
      availCap = frac_por - frac_res;

      % compute the saturated hydraulic conductivity (Colbeck 1972)
      k_snow = 3.41875e-5 * exp(15.9*(1-frac_ice));

      % if grain size is known, try Shimizu 1970 (ro_snow = mLayerTheta?)
      % k_snow = 7.7e-4*grain_sz^2*9.8/1.787*exp(-0.0078*mLayerTheta*1000);

      % compute the relative saturation
      if (availCap > 0.0)
         relSat = (frac_liq-frac_res)/availCap;
      else
         relSat = 1.0;
      end

      % calculate the flux out of the snowpack or layer
      liqflux = dt*k_snow*relSat^3.0;

      % if grain size is known, use the van Genuchten equation
      % n = (15.68*exp(-0.46*2*grain_sz))+1;
      % m = 1 - (1/n);
      % liqflux = dt*k_snow*relSat^0.5*((1-((1-relSat^(1/m))^m))^2);
   else
      % flow does not occur
      liqflux = 0.0;
   end
end
% liq_res = 0.02;                % resid. water volume / pore volume [-]
% frac_por = 1 - frac_ice;        % pore fraction                     [-]
% frac_res = liq_res * frac_por;  % residual water fraction           [-]
% liq_sat = frac_liq / frac_por; % water saturation     [-]
% ava_cap = frac_por - frac_res; % available capacity   [-]
% rel_sat = (frac_liq - frac_res) / ava_cap  % relative saturation  [-]

% clarity on the analogues with colbeck:
% Swi = liq_res = frac_res/frac_por = h_res/h_por = h_res/(h_tot-h_ice)
% Sw = liq_sat = frac_liq/frac_por = h_liq/h_por = h_liq/(h_tot-h_ice)
% phi = porosity = frac_por = h_por/h_tot = (h_tot-h_ice)/h_tot
% Sstar = rel_sat = (frac_liq-frac_res)/ava_cap = (h_liq-h_res)/ava_cap


% % NOTE: didn't clarify this, but updated march 2023 with fracwater method
% function [fracliquid,mLayerVolFracLiq,mLayerVolFracIce,retent] = ...
%    UPDATE_STATE(T,mLayerTheta,fcp,Tf,ro_liq,retent)
%
% % fraction of liquid water
%    fracliquid = 1.0/(1.0 + (fcp * (Tf - min(T,Tf)))^2.0);
%
% % % if using f_wat method, should be:
% %    Tdep = Tf-min(T,Tf);
% %    fracwater = f_liq + f_ice.*ro_iwe;                  % frac_wat_old
% %    fracliquid = fracwater./(1.0 + fcp^2.*Tdep.*Tdep);   % eq 67, Jordan
% % But note: mLayerTheta = fracwater, nearly certain
%
% % calculate volumetric fractions of liquid and ice content
%    mLayerVolFracLiq = fracliquid * mLayerTheta;
%    mLayerVolFracIce = (1.0 - fracliquid) * mLayerTheta * (ro_liq/917);
%
% % calculate the retention
%    retent = min(0.02,max(retent,0.75*mLayerVolFracLiq));
