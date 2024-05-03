function [T, f_ice, f_liq, OK] = MZTRANSFORM(T, T_old, f_ice, f_liq, f_wat, ...
      ro_ice, ro_liq, Tf, TL, TH, fcp, f_liq_min, f_liq_max, iM, OK)
   %MZTRANSFORM Liquid fraction change to temperature-enthalpy transformation
   %
   % This function uses the change in liquid fraction returned by the numerical
   % solution of the enthalpy equation to update the temperature of nodes
   % undergoing phase change, using the temperature-enthalpy relationship.
   % This ensures the change in enthalpy due to solid-liquid phase change is
   % accounted for in terms of both the latent heat and the specific heat of the
   % control volume, i.e., that the temperature, liquid fraction, and enthalpy
   % are consistent. This is sometimes referred to as a "corrector" step. Note
   % that for nodes undergoing phase change, the input to this function, T, is
   % the change in liquid fraction due to solid-liquid phase change multiplied
   % by the density of liquid water. For nodes which are not undergoing phase
   % change, T is the temperature.
   %
   % This function also implements two numerical checks:
   %  1) if a node within the melt zone exited the meltzone by more than 5%
   %     in either direction, and
   %  2) if any nodes skipped the melt zone, meaning a node below the lower
   %     melt zone temperature TL ended the step above the upper melt zone
   %     temperature TH (i.e., completely melted in one step)
   % If either are true, OK is set false, control returns to the main program,
   % the time step is shortened, and the solution is repeated until a valid
   % solution is obtained.
   %
   % Notes:
   %  - Told is needed for the second check.
   %  - f_liq_min/max are f_liq at T=TL/TH, given the current f_wat, not the
   %    min/max possible f_liq (but f_liq_max = f_wat for an ice model)
   %  - iM is the melt zone nodes at each inner iteration, f_wat is the
   %    water fraction at the start of the outer iterations.
   %
   % See also:
   %
   %#codegen

   % Check if the change in liquid water exceeds available ice
   if any( (T(iM) / ro_liq) >= (ro_ice / ro_liq - f_liq(iM)) )
      OK = false;
      return
   end

   % Update the liquid fraction of melting layers (f_liq = f_liq_o + P/ro)
   f_liq(iM) = f_liq(iM) + T(iM) / ro_liq; % line 79 of ftemp.f

   % This is what the transformation above does:
   %
   % f_liq = P / ro_liq + f_liq_o
   %       = ro_liq * (f_liq - f_liq_o) / ro_liq + f_liq_o
   %       =          (f_liq - f_liq_o)          + f_liq_o
   % f_liq = f_liq

   % Don't allow f_liq to exceed f_wat by more than 1% (line 113 of ftemp.f)
   if any((f_liq(iM) - f_wat(iM)) ./ f_wat(iM) > 0.01)
      OK = false;
      return
   end
   f_liq(iM) = min(f_liq(iM), f_wat(iM) - eps);

   % Update the ice fraction
   f_ice = (f_wat - f_liq) * ro_liq / ro_ice;

   % Check for nodes that were within the melt zone but now have a liquid
   % fraction below the lower limit (lower phase boundary overshoot).
   if any(f_liq(iM) < 0.95*f_liq_min(iM)) || any(f_liq(iM) > 1.05*f_liq_max(iM))
      OK = false;
      return
   end

   % Update the temperature (NOTE g_wat / g_liq = f_wat / f_liq)
   T(iM) = Tf - sqrt(f_wat(iM) ./ f_liq(iM) - 1.0) / fcp; % Eq. 133a

   % Note: above uses the new f_liq, which is correct. Don't use this update:
   % T(iM) = T_old(iM) + T(iM) ./ (ro_liq * dFdT(iM));

   % Check for nodes that were outside the melt zone (T<TL) and fully melted in
   % one step (upper phase boundary overshoot).
   if any(T_old(~iM) < TL & T(~iM) >= TH) ...
         || any(T_old(~iM) > TH & T(~iM) < TL) || any(iscomplex(T))
      OK  = false;
      return
   end

   % Update f_liq/ice/T to be consistent with the enthalpy-temperature
   % relationship (corrector step)
   [T, f_ice, f_liq] = MELTCURVE(T, f_ice, f_liq, ro_ice, ro_liq, fcp, Tf);
end
