function [Terr,asflag] = TEMPERROR(T,T_old,ro_sno,cp_sno,f_liq,f_air, ...
      f_liq_old,drovdT,xLs,xLf,dz,dt_old,dt_new)
   %TEMPERROR Compute the linearization error 

   tol = 0.05; % K
   
   dT = T-T_old;
   dLiqMass = ro_liq.*dz.*(f_liq-f_liq_old);                % kg/m2
   dLiqHeat = dt_old.*-dLiqMass.*xLf./dt_new;               % J/m2
   dSpcHeat = dz.*(ro_sno.*cp_sno + xLs.*f_air.*drovdT);    % J/m2/K

   Terr = dLiqHeat./dSpcHeat - dT;                          % K

   if abs(Terr)>tol
      asflag = true;
   else
      asflag = false;
   end

   % dLiqMass = ro_liq.*dz.*(f_liq-f_liq_old);                % kg/m2
   % dLiqNrg = -dLiqMass.*xLf./dt_new;                       % W/m2
   % dLiqHeat = dt_old.*dLiqNrg;                              % J/m2
end
