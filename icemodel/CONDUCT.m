function Qc = CONDUCT(k_eff, T, dz, Ts, level)
   %CONDUCT Compute conductive heat flux, positive into the surface
   %
   %  Qc = CONDUCT(k_eff, T, dz, Ts)
   %
   %  Units: [W m-2] = [W m-1 K-1] * [K] * [m-1]
   %
   % See also: ENBALANCE, SEBSOLVE, SFCTEMP, MFENERGY
   %
   %#codegen

   if nargin < 5
      level = 1;
   end

   if level == 1
      % Conduction from layer 1 into the surface (across a 1/2 cv)
      Qc = k_eff(1) * (T(1) - Ts) / (dz(1) / 2);
   elseif level == 2
      % Conduction from layer 2 into layer 1 (across a full cv)
      Qc = (k_eff(1) + k_eff(2)) / 2.0 * (T(2) - T(1)) / (dz(1) + dz(2)) / 2.0;
   end
end
