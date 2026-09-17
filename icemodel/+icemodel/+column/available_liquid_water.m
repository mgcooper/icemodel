function [h_resid, h_avail, h_drain, h_ice, h_liq, h_air] = available_liquid_water( ...
      h_ice, h_liq, h_air, T_ice_old, Tf, theta_resid, h_drain, liqflag, ...
      h_total)
   %AVAILABLE_LIQUID_WATER Compute available liquid water in one control volume.
   %
   % Note: The model does not use this function. Call it with scalar values,
   % in a loop or for the top layer. It can also work with vectors.
   %
   %#codegen

   % Update theta_resid
   if T_ice_old < Tf; theta_resid = 0.0; end

   % Compute liqresid and liqavail in units of liquid water. liqavail is
   % liquid water available to drain or freeze. liqresid is unavailable.
   h_resid = theta_resid * h_ice;
   h_avail = max(h_liq - h_resid, 0.0);

   % % these would be in an 'else' statement below but here is faster
   % availCap = porosity - h_resid;    % available capacity      [m]
   % relsat = h_avail ./ availCap      % relative saturation     [-]

   % Set the liquid drainage decision flag, with layer 1 always unsaturated.
   if liqflag
      h_drain = h_drain + h_avail;
      h_liq = max(h_liq - h_avail, h_resid);
      h_air = h_total - h_liq - h_ice;
      % h_avail = 0.0;
      % availCap = h_resid;         % available capacity      [m]
      % relsat = h_avail./availCap  % relative saturation     [-]
   end
end

% NOTE: liqflag must return false for every layer below the first layer. If it
% returns true, MELT drains all water except liqresid.

% If porosity is defined as h_tot - h_ice, Colbeck's water saturation is:
% Sw = liqsat = h_liq_j / porosity,
% If porosity is defined as 1 - f_ice, then:
% Sw = liqsat = f_liq / porosity

% theta_resid here is NOT the usual theta_resid. It is the irreducible liquid
% water content of the glaciological literature: volume water / volume ice.
% Thus:
% liqresid = theta_resid * h_ice

% Confirmed correct:
%
% theta_resid = 0.02;               % residual water [m3 water/ m3 ice] [-]
% liqresid = theta_resid * h_ice;   % residual water volume             [m]
% liqavail = h_liq - liqresid;      % water available to drain/freeze   [m]
% porosity = h_tot - h_ice;         % pore volume                       [m]
% availCap = porosity - liqresid;   % available capacity                [m]
% relsat = liqavail ./ availCap     % relative saturation               [-]

% SNTHERM, Eq. 6, "liquid saturation", s, liquid volume per void volume [m3/m3]:
% liqsat = 0.1; % made up number
% liqfrac = liqsat * porosity;
%
% Following Colbeck (might be same as SNTHERM), and using volumes not fractions:
% liqresid = 0.02;
% liqavail = porosity - liqresid;
% if liqavail > 0.0
%    relsat = (h_liq - liqresid) / liqavail;
% else
%    relsat = 1.0;
% end
