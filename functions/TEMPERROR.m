function [Terr, ok] = TEMPERROR(T, T_old, ro_sno, cp_sno, f_liq, ...
      f_air, xf_liq, drovdT, Ls, Lf, dz, dt_old, dt_new)
   
   % THIS NEEDS TO BE RECONCILED WITH THE CORRECT, FINAL VERSION AT THE BOTTOM
   % OF ICEENABL, WHICH WILL BE KEPT AS THE TOP LAYER, AND THIS CAN BE CALLED TO
   % GET ALL LAYERS
   
   %TEMPERROR Compute the linearization error 
   % 
   % Governing equation for the top layer:
   % 
   % dH / dT * dT / dt = 1 / dz * k * dT / dz + S
   % (dH_T / dT + dH_vap / dT + dH_liq / dT) * dT / dt = 1/dz * k * dT / dz + S
   % (dH_T / dT + dH_vap / dT) * dT + dH_liq = dt/dz * k * dT / dz + S * dt
   % (dH_T_dT + dH_vap_dT) * dT + dH_liq = dt/dz * k_dT_dz + S * dt
   % 
   % Switch to variables:
   % 
   % (dH_T_dT + dH_vap_dT) * dT + dH_liq = dt/dz * (k_dT_dz + S * dz)
   % (dH_T_dT + dH_vap_dT) * dT = dt/dz * (k_dT_dz + S * dz) - dH_liq 
   % dT = ( dt/dz * (k_dT_dz + S * dz) - dH_liq ) / (dH_T_dT + dH_vap_dT) 
   %
   % Thus:
   %
   % err = ( dt/dz * (k_dT_dz + S * dz) - dH_liq ) / (dH_T_dT + dH_vap_dT) - dT
   %
   % Governing equation for the surface:
   %
   % SEB = Fc + Fp * Ts + Qc - Qm 
   %
   % If Ts < Tf, there is no energy available for melting, Qm = 0, otherwise Qm
   % computed from the residual: Qm = Fc + Fp * Ts.
   
   % Note: Jordan uses the actual diagnosed top fluxes, so I need to determine
   % if Tf or Ts should be used in Fc + Fp * Ts.
   
   tol = 0.05; % K
   
   % Update necessary variables
   [~, drovdT] = VAPORHEAT(T, f_liq, f_ice, Tf, Rv, Ls);
   dHdT = cv_ice * f_ice + cv_liq * f_liq;
      
   deltaT = T(1) - T_old(1);                                      % K
   dH_liq = ro_liq * Lf * (f_liq(1) - f_liq_o(1));                % J m-3
   dH_T_dT = dHdT(1);                                             % J m-3 K-1
   dH_vap_dT = Ls * drovdT(1) * (1 - f_ice(1) - f_liq(1));        % J m-3 K-1
   k_dT_dz = a2 * (T(2) - T(1)) + Fc + Fp * ...
      (Fc + a1 * T(1)) / (a1 - Fp);                               % W m-2
   % k_dT_dz = a2 * (T(2) - T(1)) + Fc + Fp * Ts;
   dq_dz = Sc(1) * dz(1);                                         % W m-2
   dt_dz = dt/dz(1);                                              % s m-1
   
   % Top layer energy balance linearization error
   Terr = (dt_dz * (k_dT_dz + dq_dz) - dH_liq ) / (dH_T_dT + dH_vap_dT) - deltaT;
   
   % Surface energy balance linearization error
   % err0 = Fc + Fp * Ts + a1 * (T(1) - Ts);
   
   % % % % % % % % % % % % % % % % % % % % % % % % % % 
   % This was an updated version of the original but it's wrong e.g. 
   % Fc + Fp * Ts is missing, 
   
   deltaT = T - T_old;
   dSrcHeat = dt_new * (a2 * (T2 - T1) + Sc * dz);
   dLiqHeat = ro_liq * Lf * dz * (f_liq - xf_liq);             % J/m2
   dSpcHeat = dz * (ro_sno * cp_sno + Ls * f_air * drovdT);    % J/m2/K
   
   Terr = (dSrcHeat - dLiqHeat) / dSpcHeat - deltaT;
   
   ok = abs(Terr) < tol;
   
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% These notes were from testing in ICEENBAL, the key thing I learned was it's
% necessary to use the governing equation not the discretized equation and use
% T(1) - T_old(1) and f_liq(1) - f_liq_old(1) not the final values on the last
% iteration and not the total enthalpy 
%    
% If this is used, the error is small
%    deltaV = Ls * drovdT(1) * (1 - f_ice(1) - f_liq(1));
%    deltaH = dt / dz(1) * (a2 * (T(2) - T(1)) - a1 * (T(1) - Ts) + Sc(1) * dz(1));
%    deltaL = ro_liq * Lf * (f_liq(1) - f_liq_o(1)) * dz(1);
%    deltaC = dz(1) * (ro_sno(1) * cp_sno(1) + deltaV);
%    
%    err = (deltaH - deltaL) / deltaC - deltaT;
%    
% %    err = Fc + Fp * T_old(1) + a1 * (T(1) - T_old(1))
% %    
% %    (H(1) - H_old(1)) * dz(1)
% %    
% %    (H(1) - H_iter(1)) * dz(1)

   
%    % % % % % % % % % % % % % % % % % % % % % % % % % % 
%    % Below here is the original 
   
%    deltaT = T - T_old;
%    dLiqMass = ro_liq * dz .* (f_liq - xf_liq);                    % kg/m2
%    dLiqHeat = dt_old * Lf * -dLiqMass ./ dt_new;                  % J/m2
%    dSpcHeat = dz .* (ro_sno .* cp_sno + Ls .* f_air .* drovdT);   % J/m2/K
% 
%    Terr = dLiqHeat ./ dSpcHeat - deltaT;                          % K
% 
%    if abs(Terr) > tol
%       asflag = true;
%    else
%       asflag = false;
%    end
% 
%    % dLiqMass = ro_liq * dz .* (f_liq - f_liq_old);               % kg/m2
%    % dLiqNrg = -dLiqMass .* xLf / dt_new;                         % W/m2
%    % dLiqHeat = dt_old * dLiqNrg;                                 % J/m2
end
