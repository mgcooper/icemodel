function [f_liq, f_ice, f_wat, T, dFdT] = MELTCURVE(T, f_liq, f_ice, ...
      ro_wie, ro_iwe, Tf, fcp) %#codegen
   %MELTCURVE Melt fraction curve
   %
   % The change in frac_liq is the new frac_liq minus the old frac_liq.
   % Increases in frac_liq are defined as positive, so the new frac_ice is
   % the old frac_ice minus dfrac_liq scaled by volume expansion.
   %
   % DEFINITIONS:
   %  fcp = freezing-curve parameter
   %  T_dep = temperature depression (-)
   %  fliq = fraction of unfrozen water (-) (redefined as liquid water fraction)
   %
   % Mass fraction of liquid water
   %  fliquid = 1.0 ./ (1.0 + (fcp * T_dep) .^ 2.0);
   %  f_liq = f_wat ./ (1.0 + (fcp * T_dep) .^ 2.0);
   %  f_ice = (1.0 - fliquid) .* f_wat * (ro_liq / ro_ice);
   % 
   % See also:

   % Volumetric fraction of liquid water eq 67, Jordan
   T_dep = Tf - min(T, Tf);
   f_wat = f_liq + ro_iwe * f_ice;                % f_wat_old
   f_liq = f_wat ./ (1.0 + (fcp * T_dep) .^ 2.0);  % f_liq_new (eq 67, Jordan)
   f_ice = (f_wat - f_liq) * ro_wie;              % f_ice_new

   if nargout > 3
      % Compute temperature by inverting the fraction of liquid water function
      T = Tf - sqrt(f_wat ./ f_liq - 1.0) / fcp;
   end
   
   if nargout > 4
      % Differentiate the freezing curve w.r.t temperature
      dFdT = 2.0 * fcp ^ 2.0 * T_dep .* f_wat ...
         ./ (1.0 + fcp ^ 2.0 * T_dep .^ 2.0) .^ 2.0;
   end
   
   % Calculate liquid retention (need to pass prior value into the function)
   % retent = min(0.02, max(retent, 0.75 * f_liq));
end

% Using ro_ice and ro_liq instead of ro_iwe:
%     T_dep = Tf - min(T, Tf);
%     f_wat = f_liq + f_ice .* ro_ice ./ ro_liq;
%     f_liq = f_wat ./ (1.0 + (fcp * T_dep) .^ 2.0);
%     f_ice = (f_wat - f_liq) .* ro_liq ./ ro_ice;
% 
% ... rest is identical
%
% In terms of fliq = g_liq / g_wat = 1 / (1 + (fcp*T_dep) ^ 2 ):
% 
%     df_liq_dT = 2 .* T_dep .* fcp .^ 2 ./ (1 + (fcp .* T_dep) .^ 2) .^ 2;
%     T_liq = Tf - ((1.0 ./ fliq - 1.0) ./ fcp .^ 2.0) .^ (0.50);
%     T_liq = Tf - sqrt((1 ./ fliq - 1)) ./ fcp
