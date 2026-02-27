function [f_ice, f_liq, d_con, xd_sbl] = ICESUBL(f_ice, f_liq, d_con, ...
      ro_ice, ro_liq, Ls, Lv, d_pevp, f_ice_min, f_liq_resid)
   %ICESUBL Apply potential vapor-driven fraction change at the top layer
   %
   % f_ice = fraction of ice by volume in each control volume
   % f_liq = fraction of liquid water by volume in each control volume
   % d_pevp = potential vapor-driven change in top-layer liquid fraction
   % d_con = condensation which exceeds control volume available porosity
   % xd_sbl = sublimation which exceeds control volume fractional ice content
   %
   %#codegen

   % Compute top-layer liquid water and ice fraction
   f_liq_top = f_liq(1);
   f_ice_top = f_ice(1);

   % Compute top-layer residual water fraction
   f_res = f_liq_resid * (1.0 - f_ice_top);

   % Clarify the CV budget. Note for glacier ice or any medium w/"closed" pores,
   % "availableCapacity" is a misnomer, it should be "potentialCapacity", and
   % "availableCapacity" is actually the airFraction, or "availablePorosity",
   % and should be adjusted for f_bub.
   %
   % f_por = 1 - f_ice;             % porosityFraction
   % f_res = f_liq_resid * f_por;   % residualWaterFraction
   % f_air = 1 - f_ice - f_liq;     % airFraction or availablePorosity
   % f_ava = f_por - f_res;         % availableCapacity or potentialCapacity
   % f_sat = f_liq - f_res;         % availableWater
   % S_rel = f_sat / f_ava;         % relativeSaturation

   % Initialize variables
   xd_evp = 0;
   xd_sbl = 0;

   % Determine whether evaporation or sublimation occurs. Evaporation occurs if
   % f_liq is > 0.02. Otherwise sublimation (or condensation) occurs.
   if f_liq_top > 0.02 || d_pevp > 0

      if d_pevp < 0 % evaporation

         if f_liq_top < f_res % f_liq_top - f_res < 0
            % only residual water exists, send d_pevp to SUBL
            xd_evp = d_pevp;

            %fprintf('metstep = %d, f_liq(1) < f_res\n',metstep)

         elseif abs(d_pevp) <= (f_liq_top - f_res) % availWater >= evap
            % some water evaporates
            % d_aevp = d_pevp;
            f_liq(1) = f_liq(1) - abs(d_pevp); % d_aevp = d_pevp, xd_evp = 0

         else
            % all available water evaporates, residual water remains
            d_aevp = -(f_liq_top - f_res);
            f_liq(1) = f_res;
            xd_evp = d_pevp - d_aevp; % remaining energy goes to sublimation
         end

         % Send energy demand not satisfied by evaporation to SUBL
         [f_ice, xd_sbl] = SUBL(xd_evp, f_ice, f_ice_min, ro_ice, ro_liq, Lv, Ls);

      elseif d_pevp > 0 % condensation

         d_aevp = ro_ice / ro_liq * (1.0 - f_ice_top) - f_liq_top;
         if d_pevp <= d_aevp
            % all condensation stored
            f_liq(1) = f_liq(1) + d_pevp;  % d_aevp = d_pevp, xd_evp = 0

         else
            % some condensation stored, some converts to runoff
            f_liq(1) = f_liq(1) + d_aevp;
            d_con(1) = d_con(1) + (d_pevp - d_aevp); % xd_evp = d_pevp - d_aevp
            % fprintf('condensation exceeds porosity: %.6f\n', (d_pevp - d_aevp))
         end
      end

   else % sublimation
      [f_ice, xd_sbl] = SUBL(d_pevp, f_ice, f_ice_min, ro_ice, ro_liq, Lv, Ls);
   end
end

%% Sublimation
function [f_ice, xd_sbl] = SUBL(d_pevp, f_ice, f_ice_min, ro_ice, ro_liq, Lv, Ls)

   xd_sbl = 0;
   if d_pevp >= 0
      return
   end

   % Convert potential evaporation to potential sublimation:
   % d_pevp = Qe / (Lv * ro_liq) * dt_new / dz_therm;
   % d_psbl  = Qe / (Ls * ro_ice) * dt_new / dz_therm;
   % Therefore:
   % d_psbl = d_pevp * (Lv * ro_liq) / (Ls * ro_ice)

   % Note: in icemodel, roL is set to roLv or roLs depending on liqflag, so Qe
   % is already computed wrt to them. That way evap/subl are computed using the
   % same formula: e = Qe / roL * dt / dz.
   %
   % The conversion here conserves heat when the surface latent heat flux, Qe,
   % cannot be satisfied by evaporation alone (all available water evaporates),
   % and the remainder is allocated to sublimation of ice using the conversion:
   % d_sbl = xd_evp * roLv / roLs

   % Convert potential evap to potential subl in ice frac-equivalent thickness
   d_psbl = d_pevp * (Lv * ro_liq) / (Ls * ro_ice);
   f_ice_top = f_ice(1);

   % Note: layer combination is based on f_ice, so requiring f_ice_top < 0
   % should suffice (rather than <f_ice_min). If f_ice(1) + d_psbl < 0, it will
   % error, otherwise the layers will combine if f_ice(1) + d_psbl < f_min.

   if f_ice_top < f_ice_min

      if abs(d_psbl) < f_ice_top
         % some ice sublimates, and the top layer will be combined in ICEMF
         f_ice(1) = f_ice(1) - abs(d_psbl);

      elseif abs(d_psbl) >= f_ice_top
         % all ice sublimates (d_asbl = f_ice(1))
         f_ice(1) = 0.0;

         % sublimation which cannot be satisfied
         xd_sbl = d_psbl + f_ice_top;
      end

   elseif abs(d_psbl) < (f_ice_top - f_ice_min)
      % some ice sublimates (d_asbl = d_psbl)
      f_ice(1) = f_ice(1) - abs(d_psbl);

   elseif abs(d_psbl) >= (f_ice_top - f_ice_min)
      % all ice sublimates (d_asbl = -(fi - f_ice_min))
      f_ice(1) = f_ice_min; % keep minimum ice thickness

      % sublimation which cannot be satisfied
      xd_sbl = d_psbl + f_ice_top - f_ice_min;
   end
end
