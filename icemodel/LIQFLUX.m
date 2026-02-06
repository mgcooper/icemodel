function [q, dq_dtheta] = LIQFLUX(f_liq, f_ice, liqresid, grainsz)
   %LIQFLUX Compute the liquid water flux between snowpack layers.
   %
   % Parameters:
   % -----------
   % f_liq : vector
   %     Fraction of liquid water in each layer (volumetric).
   %
   % f_ice : vector
   %     Fraction of ice in each layer (volumetric).
   %
   % dt : scalar
   %     Time step (in seconds).
   %
   % liqresid : scalar (optional)
   %     Residual water volume per pore volume. Default is 0.02.
   %
   % grainsz : scalar (optional)
   %     Grain size in meters for Shimizu 1970 method. If supplied, this method
   %     will be used instead of Colbeck 1972.
   %
   % Returns:
   % --------
   % q : vector
   %     Liquid water flux (m^3/m^2/s) between the layers for the time step.
   %
   % Notes:
   % ------
   % This function can use either of two methods to compute the hydraulic
   % conductivity:
   % 1. Colbeck 1972
   % 2. Shimizu 1970, which requires grain size as an input.
   %
   % The function calculates the available capacity and relative saturation of
   % liquid water in the snowpack. The flux is then computed based on the
   % saturated hydraulic conductivity and relative saturation.
   %
   % The following analogues with Colbeck 1972 are used in this function:
   % liqresid  - residual water volume / pore volume [-] (Swi)
   % f_por     - pore fraction [-]                       (phi)
   % f_res     - residual water fraction [-]
   % liqsat    - water saturation [-]                    (Sw)
   % availcap  - available capacity [-]
   % relsat    - relative saturation [-]                 (Sstar)
   %
   % f_ice, f_liq, and f_air are volumetric fractions. f_wat is the
   % water-equivalent volume fraction, i.e., if the ice melted, the combined
   % liquid plus melted ice fraction.
   %
   % See also:
   %
   %#codegen

   % Note, in general, 2-7% of the pore space must be filled with water before
   % any can infiltrate, so when debugging, check that

   % Parse inputs
   if nargin < 3
      liqresid = 0.07;
   end

   if nargin < 4
      method = 1;
   else
      method = 2;
   end

   % Initialize
   N = numel(f_liq);
   q = zeros(N, 1);
   dq_dtheta = zeros(N, 1);

   % Compute pore and residual water fractions
   f_por = 1.0 - f_ice;
   f_res = liqresid * f_por;

   % Compute saturated hydraulic conductivity (Colbeck 1972)
   if method == 1
      k_snow = 3.4975e-7 * exp(15.9 * f_por);

      % If updating dynamic viscosity as a function of T:
      % n = icemodel.dynamicViscosityWater(T);
      % k_snow = (6.1313e-10 ./ n) .* exp(15.9 * f_por);

   else
      % If grain size is known, use Shimizu 1970
      % grainsz = 0.01;
      % n = icemodel.dynamicViscosityWater(T);
      f_wat = f_liq + 0.917 * f_ice;
      k_snow = 0.077 * 1000 * 9.81 / 0.001753 * grainsz ^ 2 * exp(-7.8 * f_wat);
   end

   % Compute available capacity and relative saturation
   availCap = max(0.0, f_por - f_res);
   relSat = (f_liq - f_res) ./ availCap;
   relSat(availCap <= 0) = 1;

   % Note: if relSat is used elsewhere, set negative values to 0 or NaN:
   % relSat(relSat < 0) = NaN;

   % Indices where flow can occur
   iflux = availCap > 0 & f_liq > f_res;

   % Compute the liquid water flux [m/s]
   if method == 1
      m = 2.0;
      q(iflux) = k_snow(iflux) .* relSat(iflux) .^ m; % * dt
   else
      % if grain size is known, use the van Genuchten equation
      % grainsz = 1e-2;
      n = (15.68 * exp(-0.46 * 2 * grainsz)) + 1;
      m = 1 - (1 / n);
      q(iflux) = k_snow(iflux) .* relSat(iflux) .^ 0.5 ...
         .* ((1 - ((1 - relSat(iflux) .^ (1/m)) .^ m)) .^ 2); % * dt
   end

   % Compute dq_dtheta
   if method == 1
      relSatDeriv = 1 ./ availCap;
      dq_dtheta(iflux) = m .* k_snow(iflux) ...
         .* relSat(iflux) .^ (m-1) .* relSatDeriv(iflux); % * dt
   else
      % Here, add dq_dtheta for Shimizu 1970 if needed
   end

   % Calculate the min/max retention (field capacity?)
   % retention = min(0.02, max(retent, 0.75 * f_liq));
end

% For reference, in terms of CV thicknesses:
% Swi = liqresid  = f_res/f_por  = h_res/h_por = h_res/(h_tot-h_ice)
% Sw  = liqsat    = f_liq/f_por  = h_liq/h_por = h_liq/(h_tot-h_ice)
% phi = porosity  = f_por        = h_por/h_tot = (h_tot-h_ice)/h_tot
%
% Sstar = relSat  = (f_liq-f_res)/availCap = (h_liq-h_res)/availCap
