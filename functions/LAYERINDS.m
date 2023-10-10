function [j1,j2] = LAYERINDS(ji,f_ice) %#codegen
%LAYERINDS compute layer indices to merge

j1 = ji; % first layer to merge

% Combine the top layer with the bottom layer
if j1 == 1
   j2 = j1 + 1;

   % Combine any other layer with its thinnest neighbor that is > 0
else
   j2a = j1 - 1; % layer above
   j2b = j1 + 1; % layer below

   % If the layer is zero, choose the thinnest neighbor that is > 0
   % If the layer is not zero, choose the neighbor that is zero.

   if f_ice(j1) == 0

      if f_ice(j2a) < f_ice(j2b)
         j2 = j2a;
      elseif f_ice(j2b) < f_ice(j2a) && f_ice(j2b) > 0.0
         j2 = j2b;
      elseif f_ice(j2b) == 0.0 && f_ice(j2a) > 0.0
         j2 = j2a;
      elseif f_ice(j2a) == 0.0 && f_ice(j2b) > 0.0
         j2 = j2b;
      else % they are both zero
         j2 = j2a;
      end

   else

      if f_ice(j2a) == 0.0 && f_ice(j2b) > 0.0
         j2 = j2a;
      elseif f_ice(j2b) == 0.0 && f_ice(j2a) > 0.0
         j2 = j2b;
      else % they are both zero
         j2 = j2a;
      end

      % finally, check if j2b is thinner than j2a
      % if h_ice(j2b)<h_ice(j2a)
      %    j2 = j2b;
      % end
   end
end
