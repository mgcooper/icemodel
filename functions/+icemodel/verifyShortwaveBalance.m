function verifyShortwaveBalance(Qsi, Sc, albedo, chi, dz, method, Qseb)

   if nargin < 7 && method == 3
      error('Qseb must be supplied with method 3')
   end
   
   % Compute absorbed radiation from forcing data
   totalNetShortwave = Qsi * (1 - albedo);
   topNodeNetShortwave = chi * Qsi * (1 - albedo);
   interiorNetShortwave = (1 - chi) * Qsi * (1 - albedo);

   % Compute absorbed radiation from the two-stream source term solution
   totalSourceTerm = sum(Sc .* dz);
   topNodeSourceTerm = Sc(1) * dz(1);
   interiorSourceTerm = sum(Sc(2:end) .* dz(2:end));
   
   % Confirm the total, top node, and interior nodes

   switch method
      
      case 1 % the standard method 
         
         assertWithRelTol(totalNetShortwave, totalSourceTerm)
         assertWithRelTol(topNodeNetShortwave, topNodeSourceTerm)
         assertWithRelTol(interiorNetShortwave, interiorSourceTerm)

      case 2 % the Sc(1) = 0 method

         assertWithRelTol( totalNetShortwave ...
            - ( topNodeNetShortwave + totalSourceTerm ), 0)
         
         assertWithRelTol(interiorNetShortwave - topNodeSourceTerm, 0)
      
      case 3 % the Qseb method
         
         % The total absorbed radiation should equal the portion allocated to
         % the surface plus the sum of the subsurface absorption
         
         % Verify the total, surface (SEB), and interior
         assertWithRelTol(totalNetShortwave, totalSourceTerm + Qseb)
         assertWithRelTol(topNodeNetShortwave, Qseb)
         assertWithRelTol(interiorNetShortwave, totalSourceTerm)

         % so we can either pass Qseb or chi out of extcoefs   
         
         % For clarity, the total absorbed:
         Qsi * (1 - albedo) * chi + sum(Sc .* dz)
         % Compare that with the known total absorbed:
         Qsi * (1 - albedo)
         
      otherwise
         
   end
   