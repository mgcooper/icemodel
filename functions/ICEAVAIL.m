function [h_ice, h_liq, h_air, potmelt, ASflag] = ICEAVAIL(h_ice, h_liq, ...
      h_air, isubl, ro_ice, ro_liq, ro_sno, cp_sno, T_old, Tf, Lf, h_tot, ...
      f_ice, th_resid, lay2subl)
   %ICEAVAIL Compute ice available to sublimate

   % Reduce h_ice(1) by surf_subl in units of ice equivalent
   h_tmp = h_ice;
   h_sub = isubl * ro_liq / ro_ice;
   %h_res = opts.theta_resid*h_ice(2);

   % If the top layer density is <600 kg/m3, sublimate a small amount from the
   % second layer down so that the layer combination is smooth.

   %if lay2subl == true && f_ice(1) * ro_ice < 600 && h_ice(1) < h_ice(2)
   if lay2subl == true && ro_sno(1) < 600 && h_ice(1)<h_ice(2)
      %if lay2subl == true && f_ice(1) * ro_ice < 600 && h_liq(2) <= h_res
      h_ice(1) = h_ice(1) - 0.8*h_sub;
      h_ice(2) = h_ice(2) - 0.2*h_sub;
   else
      h_ice(1) = h_ice(1) - h_sub;
   end

   % Calculate liquid water equivalent of the ice in each layer
   iwe = h_ice .* ro_ice./ro_liq;
   potmelt = cp_sno.*(T_old-Tf)./Lf.*h_tot;
   ASflag = any(potmelt>iwe);

   % Repeat for freeze
   % ASFflag = any(abs(potmelt(potmelt<0))>h_liq(potmelt<0));

   % If melt energy exceeds available ice, reset h_ice and return
   if ASflag == true
      h_ice = h_tmp; return;
   else
      % Calculate new air volume
      h_air(1) = h_tot - h_liq(1) - h_ice(1);
      h_air(2) = h_tot - h_liq(2) - h_ice(2);
   end
end

% % If the top layer ice density is <200 kg/m3, sublimate the second layer
% idens = ro_ice*h_ice(1)/dy_p(1);
% if idens >= 200
%    h_ice(1) = h_ice(1) - h_sub;
% else
%    h_ice(2) = h_ice(2) - h_sub;
% end
%
% % if h_ice(1) gets below h_min, sublimate layer 2, so the layers combine
% h_min = opts.h_ice_min;
% if (h_ice(1) - h_sub) < h_min
%    h_ice(2) = h_ice(2) - h_sub + h_ice(1) - h_min;
%    h_ice(1) = h_min;
%
%    %       % the two lines above implement this:
%    %         del_h1 = h_ice(1)-h_min;
%    %         del_h2 = h_sub - del_h1;
%    %         h_ice(1) = h_min;
%    %         h_ice(2) = h_ice(2) - del_h2;
% else
%    h_ice(1) = h_ice(1) - h_sub;
% end
