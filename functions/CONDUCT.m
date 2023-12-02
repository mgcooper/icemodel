function Qc = CONDUCT(k_eff, T, dz, Ts)
   %CONDUCT Compute conductive heat flux, positive into the surface
   %
   %  Qc = CONDUCT(k_eff, T, dz, Ts)
   %
   %  Units: [W m-2] = [W m-1 K-1] * [K] * [m-1]
   % 
   % See also:

   % Conduction is across a 1/2 cv between the top node and the surface node
   Qc = k_eff(1) * (T(1) - Ts) / (dz(1) / 2); 
end
