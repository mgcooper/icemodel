function [T, f_ice, f_liq, f_wat, dFdT] = liquid_fraction_function(T, ...
      f_ice, f_liq, f_wat)
   %LIQUID_FRACTION_FUNCTION Project state onto the liquid-fraction function.
   %
   % Inputs:
   %  T - control volume temperature [K]
   %  f_ice - fraction of frozen water (-) (volumetric ice fraction)
   %  f_liq - fraction of unfrozen water (-) (volumetric liquid water fraction)
   % Mass fraction of liquid water:
   %  f_ell = 1.0 ./ (1.0 + (fcp * T_dep) .^ 2.0);
   %  f_liq = f_wat ./ (1.0 + (fcp * T_dep) .^ 2.0);
   %  f_ice = (1.0 - f_ell) .* f_wat * (ro_liq / ro_ice);
   %
   % Increases in f_liq are defined as positive, so the change in f_liq is the
   % new f_liq minus the old f_liq, and the new f_ice is the old f_ice minus
   % df_liq times the ratio of liquid water density to ice density (volume
   % expansion).
   %
   % See also: icemodel.column.liquid_fraction_derivative
   %
   %#codegen

   persistent Tf ro_ice ro_liq fcp
   if isempty(Tf)
      [Tf, ro_ice, ro_liq] = icemodel.physicalConstant('Tf', 'ro_ice', 'ro_liq');
      fcp = icemodel.parameterLookup('fcp');
   end

   % Compute volumetric fraction of liquid water eq 67, Jordan
   T_dep = Tf - min(T, Tf);

   % Ensure f_wat does not exceed maximum capacity by more than eps
   if nargin < 4 || isempty(f_wat)
      f_wat = icemodel.column.water_fraction(f_ice, f_liq);
   end
   f_wat = min(f_wat, ro_ice / ro_liq); % f_wat_old
   f_liq = f_wat ./ (1.0 + (fcp * T_dep) .^ 2.0);  % f_liq_new (eq 67, Jordan)
   f_ice = (f_wat - f_liq) * ro_liq / ro_ice; % f_ice_new

   if nargout > 3
      % Compute temperature by inverting the fraction of liquid water function
      T = Tf - sqrt(f_wat ./ f_liq - 1.0) / fcp;
   end

   if nargout > 4
      % Differentiate the liquid-fraction function w.r.t temperature.
      dFdT = icemodel.column.liquid_fraction_derivative(T, [], [], f_wat);
   end
end

% In terms of f_ell = g_liq / g_wat = 1 / (1 + (fcp * T_dep) ^ 2 ):
%
% df_ell_dT = 2 * T_dep * fcp ^ 2 ./ (1 + (fcp * T_dep) .^ 2) .^ 2;
% T = Tf - ((1.0 ./ f_ell - 1.0) ./ fcp ^ 2.0) .^ 0.50;
% T = Tf - sqrt((1 ./ fliq - 1)) ./ fcp
