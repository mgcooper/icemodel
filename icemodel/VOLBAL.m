function [h_ice, h_liq, h_air, x_ice, x_liq] = VOLBAL(h_ice, h_liq, h_res, h_tot)

   % Note: this is not used in the model, it was designed for the
   % thickness-based formulation. It might work as-is for the
   % volumetric fraction-based forms, but needs to be tested.

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

   % Note: It is sufficient to do h_air = h_tot - h_ice - h_liq at the end,
   % rather than h_air = 0.0 in the if/else statements above, the former
   % usually evaluates to rounding error (~1e-18), so keep the h_air = 0
   % statements in the if/else block above.
end
