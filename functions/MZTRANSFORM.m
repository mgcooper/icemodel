function [T, f_liq, f_ice, OK] = MZTRANSFORM(T, T_old, f_liq, f_wat, ...
      ro_liq, ro_wie, Tf, TL, TH, fcp, flmin, flmax, iM, OK)
   %MZTRANSFORM Melt zone transformation
   %
   % This implements two numerical checks:
   %  1) if a node within the melt zone exited the meltzone by more than 5% 
   %     in either direction, and 
   %  2) if any nodes skipped the melt zone, meaning a node below the lower 
   %     melt zone temperature TL ended the step above the upper melt zone
   %     temperature TH (i.e., completely melted in one step) 
   % 
   % Notes: 
   %  - Told is needed for the second check.
   %  - flmin/max are f_wat at T=TL/TH, not min/max possible f_wat.
   %  - iM is the melt zone nodes at the start of the step, f_wat is the 
   %    water fraction at the start of the step.
   %
   % See also:
   
   % These must be true:
   % 
   % f_ice + f_liq * ro_liq / ro_ice <= 1
   % f_ice * ro_ice / ro_liq + f_liq <= ro_ice / ro_liq
   % f_wat <= ro_ice / ro_liq
   % 
   % If all ice melts:
   % 
   % f_wat = f_liq(new) 
   %       = f_liq(old) + df_liq
   % 
   % Thus the requirement is:
   % 
   % f_liq(old) + df_liq <= ro_ice / ro_liq
   % df_liq <= ro_ice / ro_liq - f_liq(old)
   %
   % But the implementation should probably use a small f_ice threshold to
   % disallow complete melting, and/or call combine layers.
   
   if any( (T(iM) / ro_liq) >= (1 / ro_wie - f_liq(iM)) ) 
      OK = false;
   end

   % Update liquid fraction of melting layers (f_liq = f_liq_o + P/ro)
   f_liq(iM) = f_liq(iM) + T(iM) / ro_liq; % line 79 of ftemp.f
   
   % This is what the transformation above does:
   %
   % f_liq = P / ro_liq + f_liq_o
   %       = ro_liq * (f_liq - f_liq_o) / ro_liq + f_liq_o
   %       =          (f_liq - f_liq_o)          + f_liq_o
   % f_liq = f_liq 
   
   % Don't allow f_liq to exceed f_wat (see line 113 of ftemp.f)
   if any((f_liq(iM) - f_wat(iM)) ./ f_wat(iM) > 0.01)
      OK = false;
   else
      f_liq(iM) = min(f_liq(iM), f_wat(iM) - eps);
      % This might need to be:
      % f_liq(iM) = min(f_liq(iM), f_wat(iM) - f_ice(iM) / ro_wie);
   end

   % Update frac ice
   f_ice = (f_wat - f_liq) * ro_wie;
   
   % Check for nodes that were within the melt zone but now have a liquid
   % fraction below the lower limit (phase boundary overshoot).
   if any(f_liq(iM) < 0.95 * flmin(iM)) || any(f_liq(iM) > 1.05 * flmax(iM))
      OK = false;
      return
   end
   
   % Transform melt-zone frac liq to T (NOTE g_wat / g_liq = f_wat / f_liq)
   T(iM) = Tf - sqrt(f_wat(iM) ./ f_liq(iM) - 1.0) / fcp; % Eq. 133a
   
   % Note: above uses the new f_liq, which is correct. Don't use this update:
   % T(iM) = T_old(iM) + T(iM) ./ (ro_liq * dFdT(iM));

   % Check for nodes that were outside the melt zone and fully overshot it.
   if any(T_old(~iM) < TL & T(~iM) >= TH) ...
         || any(T_old(~iM) > TH & T(~iM) < TL) || any(iscomplex(T))
      OK  = false;
   end
end
