function stability = stableSTABLEFN(Tair,Tsfc,wspd,gravity,xkappa,opts)

   % added opts
   z_obs = opts.z_obs;
   z_0 = opts.z_0;
   % added a1,z1 for clarity
   a1 = 5.3*9.4;                                    % [-]
   z1 = z_obs/z_0;                                  % [-]
   C1 = a1*(xkappa/(log(z1)))^2*sqrt(z1);           % [-]
   C2 = gravity*z_obs/(Tair*wspd^2);                % [K-1]
   B1 = 9.4*C2;                                     % [K-1]
   B2 = C1*sqrt(C2);                                % [K-1]
   B8 = B1/2.0;                                     % [-]
   stability = 1.0/((1.0 + B8*(Tair-Tsfc))^2);             % [-]
end
