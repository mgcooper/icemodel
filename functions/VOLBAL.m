function [  h_ice_j,                                                    ...
            h_liq_j,                                                    ...
            h_air_j,                                                    ...
            xice,                                                       ...
            xliq    ]   =   VOLBAL(h_ice_j,h_liq_j,h_air_j,h_res,h_tot)

% Check if ice+liq exceeds available pore space
xice = max(0.0,h_ice_j+h_res-h_tot);
xliq = max(0.0,h_ice_j+h_liq_j-h_tot);

if xice > 0 % ice + residual water exceeds pore space.

   % This means the cv can acommodate the ice, and some or no residual water
   if h_ice_j <= h_tot
      
      xliq = h_res - (h_tot - h_ice_j);   % drain excess
      xice = 0.0;
      h_res = h_tot - h_ice_j;            % reduce h_resid
      
   else % the cv cannot acommodate the ice even without residual water
      
      xliq = h_liq_j;
      xice = h_ice_j - h_tot;             % this is the error
      h_res = 0.0;                        % reduce h_resid entirely
   end
   
   h_ice_j = h_tot-h_res;                 % h_resid can be zero or +ive
   h_liq_j = h_res;                       % h_liq can be zero or +ive
   h_air_j = 0.0;
   
elseif xliq > 0 % the cv can acommodate the ice + some but not all free water

   h_liq_j = h_liq_j-xliq;                % drain free water as needed
   h_air_j = 0.0;
   
else % the cv can acommodate all the ice + free water
   
   h_air_j = h_tot - h_ice_j - h_liq_j;   % reduce the air
end

% Note: athough it is sufficient to do h_air = h_tot - h_ice - h_liq at
% the end, rather than h_air = 0.0 in the if/else statements above, the
% former usually evaluates to rounding error (~1e-18), so keep the h_air=0
% statements in the if/else block above.