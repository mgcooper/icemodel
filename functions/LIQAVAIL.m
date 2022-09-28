function [  h_resid,                                                    ...
            h_avail,                                                    ...
            h_drain,                                                    ...
            h_ice_j,                                                    ...
            h_liq_j,                                                    ...
            h_air_j ]   =   LIQAVAIL(h_ice_j,h_liq_j,h_air_j,T_old_j,  ...
                            Tfp,theta_resid,h_drain,liqflag,h_total)
%--------------------------------------------------------------------------
% Update theta_resid 
    if T_old_j < Tfp; theta_resid = 0.0; end

% Compute liqresid and liqavail in units of liquid water. liqavail is
% liquid water available to drain or freeze. liqresid is unavailable.
    h_resid         =   theta_resid*h_ice_j;
    h_avail         =   max(h_liq_j-h_resid,0.0);

% % these would be in an 'else' statement below but here is faster
  % availCap        =   porosity-h_resid;  % available capacity      [m]
  % relsat          =   h_avail./availCap  % relative saturation     [-]
    
% Set the liquid drainage decision flag, with layer 1 always unsaturated.
    if liqflag
        h_drain     =   h_drain + h_avail;
        h_liq_j     =   max(h_liq_j-h_avail,h_resid);
        h_air_j     =   h_total-h_liq_j-h_ice_j;
       %h_avail     =   0.0;
      % availCap    =   h_resid;           % available capacity      [m]
      % relsat      =   h_avail./availCap  % relative saturation     [-]
    end 

% NOTE: it's essential that liqflag be sent back as false if we're not in
% the first layer, otherwise all water except liqresid is drained in MELT
       
% with porosity defined as h_tot-h_ice, colbeck's water saturation is:
% Sw = liqsat = h_liq_j/porosity,
% with porosity defined as 1-frac_ice, then Sw = liqsat = frac_liq/porosity

% also note that my theta_resid is NOT the usual theta_resid, it is
% irreducible liquid water content as defined in the glaciological
% literature as volume water / volume ice. this is why liqresid =
% theta_resid * h_ice

% % these are correct:
% theta_resid = 0.02;                 % residual water [m3 water/ m3 ice] [-]
% liqresid    = theta_resid.*h_ice;   % residual water volume             [m]
% liqavail    = h_liq - liqresid;     % water available to drain/freeze   [m]
% porosity    = h_tot - h_ice;        % pore volume                       [m]
% availCap    = porosity - liqresid;  % available capacity                [m]
% relsat      = liqavail./availCap    % relative saturation               [-]

% ALSO see LIQFLUX

% %   the following relationship holds and could be useful (eq. 6)
% %   liquid saturation, liquid volume per void volume % [m3/m3]
%    liqsat       =   0.1;
%    liqfrac      =   liqsat*porosity_snow; 
%   
% %   liqsat and liqfrac expressions above are from Jordan. The liqresid
% %   etc. below are from Colbeck (but Jordan might have them too)
% %   this would be needed to implement liquid fluxes
%     liqresid    =   0.02;
%     liqavail    =   porosity_snow - liqresid; 
%     if liqavail > 0.0
%         relsat  =   (h_liq - liqresid) / liqavail;
%     else
%         relsat  =   1.0;
%     end
