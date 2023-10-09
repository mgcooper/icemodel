function [dFdT,f_wat] = FREEZECURVE(T,f_liq,f_ice,ro_iwe,Tf,fcp) %#codegen

    Tdep    = Tf-min(T,Tf);
%   dFdT    = 2.0.*Tdep.*fcp.^2.0./(1.0+(fcp.*Tdep).^2.0).^2.0;
    
% this is the version that uses frac_liq instead of fliq
    f_wat   = f_liq+f_ice.*ro_iwe;           % frac_wat_old
    dFdT    = 2*fcp^2.*f_wat.*Tdep./(1+(fcp^2.*Tdep.^2)).^2;