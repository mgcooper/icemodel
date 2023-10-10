function Qc = CONDUCT(k_eff,T,dz,xTsfc)
   %CONDUCT Compute conductive heat flux, positive into the surface
   %
   %  Qc = CONDUCT(k_eff,T,dz,xTsfc)
   %
   % See also:

   % Conduction is across a 1/2 cv between the top node and the surface node
   Qc = -k_eff(1)*(xTsfc-T(1))/(dz(1)/2); % [W m-2] = [W m-1 K-1] * [K] * [m-1]

   % Qc = - (k_eff(1)+k_eff(2))/2.0*(T(1)-T(2))/((dz(1)/2+dz(2))/2.0);
   % k1 = mean(k_eff(1:3));
   % Qc = -k1*(3*T_old(1)-4*T_old(2)+T_old(3))/(2*dz(1));
   % dT/dz = (T(3) - 4*T(2) + 3*T(1)) / 2*dz
end
