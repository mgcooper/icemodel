function [f_ice, f_liq, d_rof, d_sbl_err] = apply_surface_vapor_mass_change( ...
      f_ice, f_liq, d_rof, d_pevp, f_ice_min, f_res_por)
   %APPLY_SURFACE_VAPOR_MASS_CHANGE Apply the surface vapor-mass increment.
   %
   % f_ice     = fraction of ice by volume in each control volume
   % f_liq     = fraction of liquid water by volume in each control volume
   % d_rof     = condensation which exceeds control volume available porosity
   % d_pevp    = potential vapor-driven change in top-layer liquid fraction
   % d_sbl_err = vapor-driven ice change which exceeds control-volume limits
   % f_ice_min = minimum retained surface ice fraction before remeshing
   % f_res_por = residual liquid-water fraction per pore volume [-]
   %
   %#codegen

   % Clarify the CV budget for a future refactor that accounts for f_bub. For
   % glacier ice or any medium w/closed pores, "availableCapacity" is a misnomer
   % if f_por is defined as 1-f_ice since it does not account for f_bub.
   %
   % Canonical definitions for snow:
   %
   % f_por = 1 - f_ice;          % total porosity
   % f_res = max(f_res_por*f_por, f_liq_min)  % volumetric residual floor
   % f_air = 1 - f_ice - f_liq;  % air fraction
   % f_ava = f_por - f_res;      % "availableCapacity" but should exclude f_bub
   % f_sat = f_liq - f_res;      % "availableWater"
   % S_rel = f_sat / f_ava;      % "relativeSaturation"

   persistent ro_ice ro_liq Ls Lv
   if isempty(ro_ice)
      [ro_ice, ro_liq, Ls, Lv] = icemodel.physicalConstant( ...
         'ro_ice', 'ro_liq', 'Ls', 'Lv');
   end

   % Option to check if f_liq_top ever falls below f_res.
   debug = false;

   % Compute the volumetric residual-water floor and wet-surface flag.
   % The floor is max(capillary, Jordan thermodynamic minimum).
   f_res = icemodel.column.residual_water_fraction(f_ice, f_liq, f_res_por);

   % Liquid-film flag.
   wetflag = f_liq > f_res;

   % Retain initial top-layer fractions.
   f_liq_top = f_liq;
   f_ice_top = f_ice;

   % Initialize potential deposition that cannot be satisfied by the cv budget.
   d_sbl_err = 0;

   % If a liquid film is present, partition vapor exchange through the liquid
   % reservoir first. Otherwise route it directly to the ice phase so dry/cold
   % deposition forms ice rather than spurious liquid water.
   if wetflag

      if d_pevp < 0 % evaporation

         if f_liq_top < f_res % f_liq_top - f_res < 0
            % only residual water exists, send d_pevp to sublimation

            if debug == true
               fprintf('metstep = %d, f_liq(1) < f_res\n', metstep)
            end

         elseif abs(d_pevp) <= (f_liq_top - f_res) % availWater >= evap
            % some water evaporates, d_pevp fully satisfied
            d_aevp = d_pevp;
            d_pevp = 0;

            f_liq = f_liq - abs(d_aevp);
         else
            % all available water evaporates, residual water remains
            d_aevp = -(f_liq_top - f_res);
            d_pevp = d_pevp - d_aevp; % send remaining d_pevp to sublimation

            f_liq = f_res;
         end

         % Send energy demand not satisfied by evaporation to sublimation
         [f_ice, d_sbl_err] = sublimation(d_pevp, f_ice, f_liq, f_ice_min, ...
            ro_ice, ro_liq, Lv, Ls);

      elseif d_pevp > 0 % condensation

         % Compute the capacity for condensation. This is the air fraction
         % converted to liquid-water-equivalent fraction assuming the air
         % fraction is ice, to account for the expansion that would occur if the
         % air fraction filled with water and later refroze, to prevent
         % the cv budget from exceeding 1.0.
         d_aevp_max = ro_ice / ro_liq * (1.0 - f_ice_top) - f_liq_top;

         if d_pevp <= d_aevp_max
            % all condensation stored
            d_aevp = d_pevp; % d_pevp = 0

            f_liq = f_liq + d_aevp;

         else
            % some condensation stored, some converts to runoff
            d_aevp = d_aevp_max;

            f_liq = f_liq + d_aevp;

            % Update d_pevp. This is "extra" condensation that converts to
            % runoff. In practice this never occurs hence debug is disabled.
            d_pevp = d_pevp - d_aevp_max;
            d_rof = d_rof + d_pevp;

            if debug == true
               fprintf( ...
                  'condensation exceeds porosity: %.6f\n', d_pevp)
            end
         end
      end

   else % dry/cold ice: sublimation or direct deposition to ice
      [f_ice, d_sbl_err] = sublimation(d_pevp, f_ice, f_liq, f_ice_min, ...
         ro_ice, ro_liq, Lv, Ls);
   end

   if debug == true && d_sbl_err > 0
      fprintf('unsatisfied sublimation: %.6f\n', d_sbl_err)
   end
end

%% Vapor exchange with the ice phase
function [f_ice, d_sbl_err] = sublimation(d_pevp, f_ice, f_liq, f_ice_min, ...
      ro_ice, ro_liq, Lv, Ls)
   %SUBLIMATION Convert potential evaporation to sublimation.
   %
   % Convert potential evaporation to potential sublimation:
   %
   % d_pevp = Qe / (Lv * ro_liq) * dt_new / dz_therm;
   % d_psbl = Qe / (Ls * ro_ice) * dt_new / dz_therm;
   %
   % Therefore:
   %
   % d_psbl = d_pevp * (Lv * ro_liq) / (Ls * ro_ice)
   %
   % Note: in icemodel, ro_air_Lv is set to ro_air * Lv or ro_air * Ls depending
   % on liqflag, so Qe is already computed wrt to them. That way evap/subl are
   % computed using the same formula: e = Qe / ro_air_Lv * dt / dz.
   %
   % The conversion here conserves heat when the surface latent heat flux, Qe,
   % cannot be satisfied by evaporation alone (all available water evaporates),
   % and the remainder is allocated to sublimation of ice by sending the excess
   % d_pevp from the evaporation branch to this function.

   % Unlike the excess condensation case which is sent to runoff, "excess"
   % deposition cannot be satisfied by the control volume budget and is assigned
   % to d_sbl_err. Initialize this value to 0.
   d_sbl_err = 0;

   % Early return if there's no energy for sublimation.
   if d_pevp == 0
      return
   end

   % Convert potential evap to potential subl in ice frac-equivalent thickness
   d_psbl = d_pevp * (Lv * ro_liq) / (Ls * ro_ice);

   % Retain initial values for the top layer f_ice/liq.
   f_liq_top = f_liq;
   f_ice_top = f_ice;
   f_air_top = 1.0 - f_liq_top - f_ice_top;

   % Budget deposition
   if d_psbl > 0

      if f_air_top <= 0
         % no pore/air space remains for new ice
         d_sbl_err = d_psbl;

      elseif d_psbl <= f_air_top
         % all deposition can be stored as new ice
         f_ice = f_ice + d_psbl;

      else
         % fill the remaining air space and return the unsatisfied remainder
         f_ice = f_ice + f_air_top;
         d_sbl_err = d_psbl - f_air_top;
      end
      return
   end

   % Note: layer combination is based on f_ice, so requiring f_ice_top < 0
   % should suffice (rather than <f_ice_min). If f_ice(1) + d_psbl < 0, it will
   % error, otherwise the layers will combine if f_ice(1) + d_psbl < f_min.

   % Budget sublimation
   if f_ice_top < f_ice_min

      if abs(d_psbl) < f_ice_top
         % some ice sublimates, and the top layer will be combined by the
         % follow-on merge_thin_layers step
         f_ice = f_ice - abs(d_psbl);

      elseif abs(d_psbl) >= f_ice_top
         % all ice sublimates (d_asbl = f_ice(1))
         f_ice = 0.0;

         % sublimation which cannot be satisfied
         d_sbl_err = d_psbl + f_ice_top;
      end

   elseif abs(d_psbl) < (f_ice_top - f_ice_min)
      % some ice sublimates (d_asbl = d_psbl)
      f_ice = f_ice - abs(d_psbl);

   elseif abs(d_psbl) >= (f_ice_top - f_ice_min)
      % all ice sublimates (d_asbl = -(fi - f_ice_min))
      f_ice = f_ice_min; % keep minimum ice thickness

      % sublimation which cannot be satisfied
      d_sbl_err = d_psbl + f_ice_top - f_ice_min;
   end
end
