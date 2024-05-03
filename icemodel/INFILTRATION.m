function [f_liq, f_ice, T, OK] = INFILTRATION(f_liq, f_ice, T, ro_ice, ...
      ro_liq, Tf, TL, fcp, dz, dt)
   %INFILTRATION
   %
   %#codegen

   OK = true;
   if all(T(2:end) < TL)
      OK = true;
      return
   end

   % Calculate fluxes between layers [m/s]
   q = LIQFLUX(f_liq(1:end-1), f_ice(1:end-1));

   % If this is inside a sub iteration, scale non-zero surface flux by dt_new
   % bc_N = q_sfc / subiter;

   % Boundary conditions
   bc_N = 0; % Infiltration rate from precipitation or surface water
   bc_S = 0; % No-flux condition at the bottom

   % Append bc's to make q the same size as f_liq
   q = [bc_N; q; bc_S];

   % No flux into sub-freezing layers
   q([T < TL ; true]) = 0.0;

   % Compute net liquid water flux for each layer d(f_liq) = d/dz q(z) * dt [-]
   df_liq = (q(1:end-1) - q(2:end)) ./ dz * dt;

   % Limit drainage and infill to available water / capacity
   f_por = 1.0 - f_ice;
   f_res = 0.07 .* f_por;
   for n = 1:numel(df_liq)
      if df_liq(n) == 0
         continue
      elseif df_liq(n) < 0
         % Calculate the maximum allowable drainage from the layer
         max_drainage = f_liq(n) - f_res(n);
         df_liq(n) = -min(abs(df_liq(n)), max_drainage);
      else
         % Don't allow more water to infiltrate than can be stored
         max_infill = ro_ice / ro_liq * (1.0 - f_ice(n)) - f_liq(n);
         df_liq(n) = min(df_liq(n), max_infill);
      end
   end

   % Rebalance downward fluxes where limited by max_infill in loop above
   for n = 1:(numel(df_liq) - 1)
      if df_liq(n) == 0
         continue
      elseif df_liq(n) < 0
         % Ensure drainage does not exceed available capacity in layer below
         % max_infill = 1.0 - f_ice(n+1) - f_liq(n+1);
         max_infill = 0.917 * (1.0 - f_ice(n+1)) - f_liq(n+1);
         if abs(df_liq(n)) > max_infill
            df_liq(n) = -max_infill; % extra_liq = abs(df_liq(n)) - max_infill;
         end
      end
   end

   % First condition should never trigger
   if any(-df_liq > f_liq)
      OK = false;
      return
   else
      % Update f_liq using explicit scheme
      f_liq = f_liq + df_liq;

      % for debugging
      xT = T; xf_liq = f_liq; xf_ice = f_ice;

      % Update the water fraction and temperature
      f_wat = f_liq + f_ice * ro_ice / ro_liq;
      T = Tf - sqrt((f_wat ./ f_liq - 1.0)) ./ fcp;

      % f_liq = f_wat ./ (1.0 + (fcp * (Tf - min(T, Tf))) .^ 2.0); % f_liq_new
      f_ice = (f_wat - f_liq) * ro_liq / ro_ice; % f_ice_new
   end

   if any(f_ice <= 1e-8)
      OK = false;
   end

   if any(T > Tf) || any(f_liq + f_ice - 1.0 > 1e-12)
      OK = false;
   end
end

% max(f_liq + f_ice - 1.0)
% balance = f_liq + f_ice - 1.0;
% find(balance > 10*eps)

% dT = T - xT;
% df_liq = f_liq - xf_liq;
% df_ice = f_ice - xf_ice;
%
% df_liq_m = 2102 .* dT ./ 334000 .* dz;

% After the BCs are appended to q:

% ----- q(1) = flux across the top (surface) interface into the first CV
%
%
%
% ----- q(2) = flux across the second interface into the second CV
%
%
%
% ----- q(3) = flux across the third interface into the third CV
%   .
%   .
%   .
% ----- q(N) = flux across the second from bottom interface into the bottom CV
%
%
%
% ----- q(N+1) = flux across the bottom interface = 0.0

% % Stability check 1: available water
% dq = (q(1:end-1) - q(2:end));
% df_liq = dq .* dt ./ dz;
% notOK = -df_liq > f_liq;
%
% % Stability check 2: storage capacity
%
%

% dq  = q(1:end-1) - q(2:end)
% -dq = q(2:end)   - q(1:end-1)
% dq < 0 means water drains out of CV
% dq > 0 means water enters CV
% dq < 0 && -df_liq > f_liq means more drains than exists = unstable
% dq > 0 && dq + f_liq + f_ice > 1 means more drains than can be stored

% Note, better than checking < TL, check if infiltrating water satisfies the
% cold content:
% f_cc = 0.0058 .* f_ice .* (Tf - T);

% if any(df_liq > 0.0 & T < TL)
%    % f_cc = 0.0058 .* f_ice .* (Tf - T);
%    % df_liq(df_liq > 0.0 & T < TL)
%    df_liq(T < TL) = 0.0;
% end

