function [T,f_liq,f_ice,OK] = MZTRANSFORM(T,T_old,f_liq,f_wat,ro_liq,   ...
      ro_wie,Tf,TL,TH,fcp,flmin,flmax,iM,OK)
   %MZTRANSFORM melt zone transformation
   %
   % this implements two numerical checks:
   % 1) if a node within the melt zone exited the meltzone by more than 5% in
   % either direction, and 2) if any nodes skipped the melt zone, meaning a
   % node below the lower melt zone temperature TL ended the step above the
   % upper melt zone temperature TH (i.e., completely melted in one step)
   % note: Told is needed for the second check
   % iM is the melt zone nodes at the start of the step, f_wat is the water
   % fraction at the start of the step
   % note: flmin/max are f_wat at T=TL/TH, NOT min/max possible f_wat

   % transform melt-zone layers from T to liq density
   f_liq(iM) = T(iM)./ro_liq + f_liq(iM); % line 79 of ftemp.f

   % transform melt-zone layers frac liq to T (NOTE g_wat/g_liq=f_wat/f_liq)
   T(iM) = Tf-sqrt(f_wat(iM)./f_liq(iM) - 1.0)./fcp;

   % Check if phase boundary is overshot. This is essential because melt-zone
   % switches will not give a valid solution outside the melt zone. Also,
   % since this is only called on nodes that were within the melt zone, it's
   % sufficient to check if the boundary is overshot in either direction.
   if any(f_liq(iM) < 0.95*flmin(iM)) || any(f_liq(iM) > 1.05*flmax(iM))
      OK = false;
   end

   % Check if any frozen nodes totally skipped the melt-zone or vise-versa
   if any(T_old(~iM)<TL & T(~iM)>TH) || any(T_old(~iM)>TH & T(~iM)<TL)
      OK  = false;
   end

   % update frac ice
   f_ice = (f_wat-f_liq).*ro_wie;
end
