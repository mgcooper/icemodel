function [f_liq,f_ice,f_wat,T,dFdT] = MELTCURVE(T,f_liq,f_ice,ro_wie, ...
      ro_iwe,Tf,fcp) %#codegen
   %MELTCURVE Melt fraction curve
   %
   % The change in frac_liq is the new frac_liq minus the old frac_liq.
   % Increases in frac_liq are defined as positive, so the new frac_ice is
   % the old frac_ice minus dfrac_liq scaled by volume expansion.
   %
   % DEFINITIONS:
   %  fcp = freezing-curve parameter
   %  Tdep = temperature depression (-)
   %  fliq = fraction of unfrozen water (-) (redefined as liquid water fraction)
   %
   %
   % See also:

   Tdep = Tf-min(T,Tf);
   f_wat = f_liq+f_ice.*ro_iwe;                 % frac_wat_old
   f_liq = f_wat./(1.0+fcp^2.*Tdep.*Tdep);      % eq 67, Jordan
   f_ice = (f_wat-f_liq).*ro_wie;

   % differentiate the freezing curve w.r.t temperature
   dFdT = 2*fcp^2.*f_wat.*Tdep./((1+fcp^2.*Tdep.*Tdep).^2);

   % compute temperature by inverting the fraction of liquid water function
   T = Tf-sqrt((f_wat./f_liq-1))./fcp;
end

% in terms of fliq = gam_liq/gam_wat:
%df_liq_dT = 2.*Tdep.*fcp.^2./(1+(fcp.*Tdep).^2).^2;
%T_liq = Tf-((1.0./fliq-1.0)./fcp.^2.0).^(0.50);
%T_liq = Tf-sqrt((1./fliq-1))./fcp


% EXPLICIT VERSION:
%     Tdep = Tf-min(T,Tf);
%     f_wat = f_liq+f_ice.*ro_ice./ro_liq;                  % frac_wat_old
%     f_liq = f_wat./(1.0+(fcp.*Tdep).^2.0);              % eq 67, Jordan
%     f_ice = (f_wat - f_liq) .* ro_liq./ro_ice;
%
% % differentiate the freezing curve w.r.t temperature
%     dFdT = (2*fcp^2).*f_wat.*Tdep./((1+(fcp.*Tdep).^2).^2);
%
% % compute temperature by inverting the fraction of liquid water function
%     T_liq = Tf-sqrt((f_wat./f_liq-1))./fcp;
%
% % in terms of fliq = gam_liq/gam_wat:
%    %df_liq_dT = 2.*Tdep.*fcp.^2./(1+(fcp.*Tdep).^2).^2;
%    %T_liq = Tf-((1.0./fliq-1.0)./fcp.^2.0).^(0.50);
%    %T_liq = Tf-sqrt((1./fliq-1))./fcp
