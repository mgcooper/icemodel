function [Terr, asflag] = TEMPERROR(T, T_old, ro_sno, cp_sno, f_liq, ...
      f_air, xf_liq, drovdT, Ls, Lf, dz, dt_old, dt_new)
   %TEMPERROR Compute the linearization error 

   tol = 0.05; % K
   % % % % % % % % % % % % % % % % % % % % % % % % % % 
   
   deltaT = T - T_old;
   dSrcHeat = dt_new * (Qc + Sc * dz);
   dLiqHeat = ro_liq * Lf * dz * (f_liq - xf_liq);             % J/m2
   dSpcHeat = dz * (ro_sno * cp_sno + Ls * f_air * drovdT);    % J/m2/K
   
   Terr = (dSrcHeat - dLiqHeat) / dSpcHeat - deltaT;
   
   
   % % % % % % % % % % % % % % % % % % % % % % % % % % 
   deltaT = T - T_old;
   dLiqMass = ro_liq * dz .* (f_liq - xf_liq);                    % kg/m2
   dLiqHeat = dt_old * Lf * -dLiqMass ./ dt_new;                  % J/m2
   dSpcHeat = dz .* (ro_sno .* cp_sno + Ls .* f_air .* drovdT);   % J/m2/K

   Terr = dLiqHeat ./ dSpcHeat - deltaT;                          % K

   if abs(Terr) > tol
      asflag = true;
   else
      asflag = false;
   end

   % dLiqMass = ro_liq.*dz.*(f_liq-f_liq_old);                % kg/m2
   % dLiqNrg = -dLiqMass.*xLf./dt_new;                       % W/m2
   % dLiqHeat = dt_old.*dLiqNrg;                              % J/m2
end
