function S = STABLEFN(Ta, Ts, wspd, scoef)
   %STABLEFN Compute the stability function
   %
   % Following Appendix 1 of Liston et al. 1999:
   %    S = 1.0;
   %    Ri = 9.81 * z_obs / (Ta * wspd ^ 2) * (Ta - Ts);
   %    De = wcoef * wspd;
   %    gamma = 5.3 * 9.4 * De / wspd * sqrt((z_obs / z_0));
   %    if Ri<0
   %       S = 1 - 9.4 * Ri / (1 + gamma * sqrt(abs(Ri)));
   %    elseif Ri > 0
   %       S = 1 / (1 + 4.7 * Ri) ^2;
   %    end
   %
   % For unstable, the absolute value of the Richardson number is implicit
   % because Ts>Ta implies Ts-Ta>0, see Appendix of Liston et al 1999
   %
   % Units: [1]
   % 
   % See also:

   S = 1.0; % Neutrally stable case.

   if (Ts > Ta) % Unstable case.

      S = 1 + scoef(2) / (Ta * wspd ^ 2) * (Ts - Ta) ...
         / (1 + scoef(3) / (sqrt(Ta) * wspd) * sqrt(Ts - Ta));

   elseif (Ts < Ta) % Stable case.

      S = 1 / ((1 + scoef(2) / (2 * wspd ^ 2) * (Ta - Ts) / Ta) ^ 2);
      
      % S = 1 / ((1 + scoef(2) / (Ta * wspd ^ 2) / 2 * (Ta - Ts)) ^ 2);
   end

end

% if i passed in z_obs, I could compute Ri, and then unstable S:
% S = 1-scoef(2)/(Ta*wspd^2)*(Ta-Ts)/(1+scoef(1)*sqrt(9.81*2/(Ta*wspd^2)*(Ts-Ta)));

% had these whiteboard notes might need (but check scoef indexing if used)
% B1 = 9.4*gravity*zobs/(Ta*wspd^2) = scoef(2)/(Ta*wspd^2)
% B2 = 5.3*9.4*wcoef*sqrt(zobs/z0)*sqrt(gravity*zobs/(Ta*wspd^2));

% % for clarity:
% a1 = 5.3 * 9.4;                                    % [-]
% z1 = z_obs / z_0;                                  % [-]
% C1 = a1 * (kappa / log(z1)) ^ 2 * sqrt(z1);        % [-]
% C2 = gravity * z_obs / (Ta * wspd ^ 2);            % [K-1]
% B1 = 9.4 * C2;                                     % [K-1]
% B2 = C1 * sqrt(C2);                                % [K-1]
% 
% if (Ts > Ta)                                       % Unstable case.
%    B3 = 1.0 + B2 * sqrt(Ts - Ta);                  % [-]
%    S = 1.0 + B1 * (Ts - Ta) / B3;                  % [-]
% elseif (Ts < Ta)                                   % Stable case.
%    B8 = B1 / 2.0;                                  % [-]
%    S = 1.0 / (1.0 + B8 * (Ta - Ts)) ^ 2;           % [-]
% else                                               % Neutrally stable case.
%    S = 1.0;                                        % [-]
% end
