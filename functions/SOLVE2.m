%--------------------------------------------------------------------------
%   NEWTON-RHAPSON SOLVER for TSFC
%--------------------------------------------------------------------------

function [Tsfc,Tflag] = SOLVE(Tair,ea,AAA,CCC,DDD,EEE,FFF,B1,B2,Tf,xTsfc)
%--------------------------------------------------------------------------

    tol                 =   1.0e-4;
    maxiter             =   20;
    old                 =   Tair;
    
 if xTsfc >= Tf
    % Over water.
    A                   =   6.1121 * 100.0;
    B                   =   17.502;
    C                   =   240.97;
 else
    % Over ice.
    A                   =   6.1115 * 100.0;
    B                   =   22.452;
    C                   =   272.55;
 end
 
    for i=1:maxiter
        
% added dif1, dif2, and dif3 for clarity
        dif1            =   Tair - old;
        dif2            =   old - Tair;
        dif3            =   old - Tf;
        
% This accounts for an increase in turbulent fluxes under unstable conditions.
        other1          =   AAA * dif1;
        es0             =   A * exp((B*dif3)/(C+dif3));
        other2          =   FFF*CCC*(ea-es0);
        
        dother1         =   -AAA;
        dother2         =   -FFF*CCC*es0*B*C/((C+dif3)^2);
        
        if (old>Tair)       % Unstable case.
            B3          =   1.0 + B2*sqrt(dif2);
            stability   =   1.0 + B1*(dif2)/B3;
            dstability  =   B1/B3 - (B1*B2*(dif2))/(2.0*B3*B3*sqrt(dif2));
            fprime1     =   -4.0*DDD*old^3;
            fprime2     =   stability * dother1 + other1 * dstability;
            fprime3     =   stability * dother2 + other2 * dstability;
            fprime4     =   -0.0;
        elseif (old<Tair)   % Stable case.
            B8          =   B1 / 2.0;
            stability   =   1.0 / ((1.0 + B8 * (dif1))^2);
            dstability  =   2.0 * B8 / ((1.0 + B8 * (dif1))^3);
            fprime1     =   -4.0*DDD*old^3;
            fprime2     =   stability * dother1 + other1 * dstability;
            fprime3     =   stability * dother2 + other2 * dstability;
            fprime4     =   -0.0;
        else                % Neutrally stable case.
            stability   =   1.0;
            fprime1     =   -4.0*DDD*old^3;
            fprime2     =   dother1;
            fprime3     =   dother2;
            fprime4     =   -0.0;
        end
        
        % added f1-5 for clarity
        f1              =   EEE - DDD*old^4;                % net radiation + conduction (Qc is in EEE)
        f2              =   AAA*(dif1)*stability;           % sensible flux
        f3              =   FFF*CCC*(ea-es0)*stability;     % latent flux
        f4              =   0.0;                            % holdover for the case where Qc is taken out from EEE
        
        funct           =   f1 + f2 + f3 + f4;
        fprime          =   fprime1 + fprime2 + fprime3 + fprime4;
        
        Tsfc            =   old - funct/fprime;
        
        % if xnew converges pass it to SFCTEMP
        if (abs(Tsfc - old)<tol) 
            Tflag   =   0;
            return
        end
        
        old             =   Tsfc;
        Tflag           =   1;
        Tsfc            =   Tair;
    end