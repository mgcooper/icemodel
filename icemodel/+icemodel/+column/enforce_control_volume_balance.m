function [h_ice, h_liq, h_air, x_ice, x_liq] = enforce_control_volume_balance( ...
      h_ice, h_liq, h_res, h_tot)
   %ENFORCE_CONTROL_VOLUME_BALANCE Enforce the total-volume constraint.
   %#codegen
   % The model does not use this function. It is written for the
   % thickness-based formulation. It can also work for the volumetric
   % fraction-based forms, but no test covers that case.

   % Check if ice+liq exceeds available pore space
   x_ice = max(0.0, h_ice + h_res - h_tot);
   x_liq = max(0.0, h_ice + h_liq - h_tot);

   if x_ice > 0
      % ice + residual water exceeds pore space.

      if h_ice <= h_tot
         % The cv can acommodate the ice, and some or no residual water

         x_liq = h_res - (h_tot - h_ice); % drain excess
         x_ice = 0.0;
         h_res = h_tot - h_ice;           % reduce h_resid

      else
         % The cv cannot acommodate the ice even without residual water
         x_liq = h_liq;
         x_ice = h_ice - h_tot;           % this is the error
         h_res = 0.0;                     % reduce h_resid entirely
      end

      h_ice = h_tot - h_res;              % h_resid can be zero or +ive
      h_liq = h_res;                      % h_liq can be zero or +ive
      h_air = 0.0;

   elseif x_liq > 0
      % The cv can acommodate the ice + some but not all free water

      h_liq = h_liq - x_liq;              % drain free water as needed
      h_air = 0.0;

   else
      % The cv can acommodate all the ice + free water
      h_air = h_tot - h_ice - h_liq;      % reduce the air
   end

   % A single h_air = h_tot - h_ice - h_liq at the end also works instead of
   % h_air = 0.0 in the if/else statements above. That form usually gives
   % rounding error (~1e-18), so keep the h_air = 0.0 statements above.
end
