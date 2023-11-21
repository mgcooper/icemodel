function [T,f_liq,f_ice,OK] = MZTRANSFORM(T,T_old,f_liq,f_wat,ro_liq,   ...
      ro_wie,Tf,TL,TH,fcp,flmin,flmax,iM,OK)
   %MZTRANSFORM melt zone transformation
   %
   % This implements two numerical checks:
   % 1) if a node within the melt zone exited the meltzone by more than 5% in
   % either direction, and 2) if any nodes skipped the melt zone, meaning a
   % node below the lower melt zone temperature TL ended the step above the
   % upper melt zone temperature TH (i.e., completely melted in one step)
   % note: Told is needed for the second check
   % iM is the melt zone nodes at the start of the step, f_wat is the water
   % fraction at the start of the step
   %
   % note: flmin/max are f_wat at T=TL/TH, NOT min/max possible f_wat

   % Update liquid fraction of melting layers (f_liq = f_liq_o + P/ro)
   f_liq(iM) = f_liq(iM) + T(iM) / ro_liq; % line 79 of ftemp.f

   % Update frac ice
   f_ice = (f_wat - f_liq) * ro_wie;
   
   % Transform melt-zone frac liq to T (NOTE g_wat/g_liq=f_wat/f_liq)
   T(iM) = Tf - sqrt(f_wat(iM) ./ f_liq(iM) - 1.0) / fcp;
   
   % Note: above uses the new f_liq, which is correct. Don't use this update:
   % T(iM) = T_old(iM) + T(iM) ./ (ro_liq * dFdT(iM));
   
   % Check for nodes that were within the melt zone but now have a liquid
   % fraction below the lower limit (phase boundary overshoot).
   if any(f_liq(iM) < 0.95 * flmin(iM)) || any(f_liq(iM) > 1.05 * flmax(iM))
      OK = false;
   end

   % Check for nodes that were outside the melt zone and fully overshot it.
   if any(T_old(~iM) < TL & T(~iM) >= TH) || any(T_old(~iM) > TH & T(~iM) < TL)
      OK  = false;
   end
end
