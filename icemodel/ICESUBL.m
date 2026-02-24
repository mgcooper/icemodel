function [f_ice, f_liq, d_con, xsubl] = ICESUBL(f_ice, f_liq, d_con, ...
      ro_ice, ro_liq, Ls, Lv, dz_therm, dt_new, Qe, f_ice_min, f_liq_resid)
   %SUBL Compute vapor mass flux for the top layer control volume
   %
   % f_ice = fraction of ice by volume in each control volume
   % f_liq = fraction of liquid water by volume in each control volume
   % d_con = condensation which exceeds control volume available porosity
   % xsubl = sublimation which exceeds control volume fractional ice conent
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

   % Compute potential evaporation
   pevap = Qe / (Lv * ro_liq) * dt_new / dz_therm;

   % Initialize variables
   xevap = 0;
   xsubl = 0;

   % Determine whether evaporation or sublimation occurs. Evaporation occurs if
   % f_liq is > 0.02. Otherwise sublimation (or condensation) occurs.
   if f_liq_top > 0.02 || pevap > 0

      if pevap < 0 % evaporation

         if f_liq_top < f_res % f_liq_top - f_res < 0
            % only residual water exists, send pevap to SUBL
            xevap = pevap;

            %fprintf('metiter = %d, f_liq(1) < f_res\n',metiter)

         elseif abs(pevap) <= (f_liq_top - f_res) % availWater >= evap
            % some water evaporates
            % d_evap = pevap;
            f_liq(1) = f_liq(1) - abs(pevap); % aevap = pevap, xevap = 0

         else
            % all available water evaporates, residual water remains
            aevap = -(f_liq_top - f_res);
            f_liq(1) = f_res;
            xevap = pevap - aevap; % remaining energy goes to sublimation
         end

         % Send energy demand not satisfied by evaporation to SUBL
         [f_ice, xsubl] = SUBL(xevap, f_ice, f_ice_min, ro_ice, ro_liq, Lv, Ls);

      elseif pevap > 0 % condensation

         aevap = ro_ice / ro_liq * (1.0 - f_ice_top) - f_liq_top;
         if pevap <= aevap
            % all condensation stored
            f_liq(1) = f_liq(1) + pevap;  % aevap = pevap, xevap = 0

         else
            % some condensation stored, some converts to runoff
            f_liq(1) = f_liq(1) + aevap;
            d_con(1) = d_con(1) + (pevap - aevap); % xevap = pevap - aevap
            % fprintf('condensation exceeds porosity: %.6f\n', (pevap - aevap))
         end
      end

   else % sublimation
      [f_ice, xsubl] = SUBL(pevap, f_ice, f_ice_min, ro_ice, ro_liq, Lv, Ls);
   end
end

%% Sublimation
function [f_ice, xsubl] = SUBL(pevap, f_ice, f_ice_min, ro_ice, ro_liq, Lv, Ls)

   xsubl = 0;
   if pevap >= 0
      return
   end

   % Convert potential evaporation to potential sublimation:
   % pevap = Qe / (Lv * ro_liq) * dt_new / dz_therm;
   % psubl = Qe / (Ls * ro_ice) * dt_new / dz_therm;
   % Therefore:
   % psubl = pevap * (Lv * ro_liq) / (Ls * ro_ice)

   % Note: in icemodel, roL is set to roLv or roLs depending on liqflag, so Qe
   % is already computed wrt to them. That way evap/subl are computed using the
   % same formula: e = Qe / roL * dt / dz.
   %
   % The conversion here conserves heat when the surface latent heat flux, Qe,
   % cannot be satisfied by evaporation alone (all available water evaporates),
   % and the remainder is allocated to sublimation of ice using the conversion:
   % s = xe * roLv / roLs

   % Convert potential evap to potential subl in ice frac-equivalent thickness
   psubl = pevap * (Lv * ro_liq) / (Ls * ro_ice);
   f_ice_top = f_ice(1);

   % Note: layer combination is based on f_ice, so requiring f_ice_top < 0
   % should suffice (rather than <f_ice_min). If f_ice(1)+psubl<0, it will
   % error, otherwise the layers will combine if f_ice(1)+psubl<f_min.

   if f_ice_top < f_ice_min

      if abs(psubl) < f_ice_top
         % some ice sublimates, and the top layer will be combined in ICEMF
         f_ice(1) = f_ice(1) - abs(psubl);

      elseif abs(psubl) >= f_ice_top
         % all ice sublimates (asubl = f_ice(1) )
         f_ice(1) = 0.0;

         % sublimation which cannot be satisfied
         xsubl = psubl + f_ice_top;
      end

   elseif abs(psubl) < (f_ice_top - f_ice_min)
      % some ice sublimates (asubl = psubl)
      f_ice(1) = f_ice(1) - abs(psubl);

   elseif abs(psubl) >= (f_ice_top - f_ice_min)
      % all ice sublimates (asubl = -(fi - f_ice_min) )
      f_ice(1) = f_ice_min; % keep minimum ice thickness

      % sublimation which cannot be satisfied
      xsubl = psubl + f_ice_top - f_ice_min;
   end
end
