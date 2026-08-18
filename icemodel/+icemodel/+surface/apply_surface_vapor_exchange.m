function [f_ice, f_liq, d_rof, d_sbl_err, d_applied] = ...
      apply_surface_vapor_exchange( ...
      f_ice, f_liq, d_rof, d_pevp, f_ice_min, f_res_por)
   %APPLY_SURFACE_VAPOR_EXCHANGE Apply the surface vapor exchange to the
   % top cell.
   %
   % f_ice     = fraction of ice by volume in each control volume
   % f_liq     = fraction of liquid water by volume in each control volume
   % d_rof     = condensation which exceeds control volume available porosity
   % d_pevp    = potential vapor-driven change in the top cell's liquid
   %             fraction, an energy demand the surface energy balance fixed
   % d_sbl_err = per-cell vapor-driven change which exceeds control-volume
   %             limits, in ice-fraction units; positive is rejected
   %             deposition, negative is unsatisfied sublimation
   % d_applied = the exchange the top cell realized, signed like d_pevp, as
   %             the liquid-water volume fraction of the mass that actually
   %             crossed the surface after every limit, including
   %             condensation overflow, which crossed before it ran off.
   %             Grain growth consumes this: only demand the limits
   %             rejected never crossed, so only that must not drive
   %             growth.
   % f_ice_min = minimum retained ice fraction before remeshing
   % f_res_por = residual liquid-water fraction per pore volume [-]
   %
   % Three limits: condensation into pore liquid capped by pore capacity,
   % deposition into ice capped by available air space, and sublimation
   % floored at f_ice_min. The limits clamp and record; they never reject the
   % substep. The fluxes are small, and the ok rejection machinery exists for
   % an inconsistent solve rather than for a limited source.
   %
   % Condensation the top cell cannot hold leaves as runoff, because surface
   % condensate physically runs off.
   %
   % This applies the surface exchange, which the surface energy balance
   % fixes as an ENERGY: the same Qe sublimates less mass than it evaporates,
   % so this function may spend part of one demand on liquid and the rest on
   % ice. icemodel.column.apply_vapor_transport applies the
   % interior transport instead, which Fick's law fixes as a MASS. Sending
   % interior transport through here would let a limited cell apply a
   % different mass than arrived.
   %
   %#codegen

   % Reference definitions for the control-volume budget. For glacier ice or
   % any medium w/closed pores, "availableCapacity" is a misnomer if f_por is
   % defined as 1-f_ice since it does not account for f_bub.
   %
   % Canonical definitions for snow:
   %
   % f_por = 1 - f_ice;          % total porosity
   % f_res = max(f_res_por*f_por, f_liq_min)  % volumetric residual floor
   % f_air = 1 - f_ice - f_liq;  % air fraction
   % f_ava = f_por - f_res;      % "availableCapacity" but should exclude f_bub
   % f_sat = f_liq - f_res;      % "availableWater"
   % S_rel = f_sat / f_ava;      % "relativeSaturation"
   %
   % Settled here, per DesignSpec decision 6. The model carries no bubble
   % fraction: f_ice and f_liq are the only phase state, so 1 - f_ice is the
   % only porosity the column can compute. The arithmetic below is therefore
   % correct for what the model represents, and the name is what overstates
   % it. Read f_por and the condensation capacity as open porosity in snow.
   % In bubbly glacier ice they are the total non-ice volume. Some of that
   % volume is closed to liquid there, so the capacity is an upper bound.
   %
   % Tracking f_bub would need a new state variable, a densification path
   % that produces it, and a rule for when bubbles close off. That belongs to
   % the snow-firn epic, not to this function. Until then the overstatement
   % is bounded and visible: the condensation limit binds only when a cell is
   % nearly saturated, and every amount it rejects is recorded rather than
   % dropped.

   % The surface exchange reaches the top cell. d_sbl_err carries one entry
   % per cell so the ledger integrates one shape whichever path filled it,
   % but only the first entry can be nonzero here.
   d_sbl_err = zeros(numel(f_ice), 1);

   [f_ice(1), f_liq(1), overflow, d_sbl_err(1), d_applied] = applyCell( ...
      f_ice(1), f_liq(1), d_pevp, f_ice_min, f_res_por);

   % Surface condensate the top cell cannot hold runs off.
   d_rof = d_rof + overflow;
end

function [f_ice, f_liq, d_rof, d_sbl_err, d_applied] = applyCell( ...
      f_ice, f_liq, d_pevp, f_ice_min, f_res_por)
   %APPLYCELL Apply one cell's vapor-mass increment under the three limits.
   %
   % D_ROF is the condensation this cell could not hold. The caller decides
   % where it goes, because that depends on whether the cell is the surface.
   %
   % D_APPLIED is the realized exchange as a signed liquid-water volume
   % fraction, formed from the state the limits actually produced: the
   % liquid change plus the ice change converted to the same mass basis.

   persistent ro_ice ro_liq
   if isempty(ro_ice)
      [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   end

   % Option to check if f_liq_top ever falls below f_res.
   debug = false;

   % This cell holds no overflow until the condensation branch finds some.
   d_rof = 0;

   % Liquid-film flag and the volumetric residual-water floor the evaporation
   % branch draws down to. The floor is max(capillary, Jordan thermodynamic
   % minimum). One function owns this decision. The coupled path has to pick
   % the same phase when it converts transported mass into the energy
   % demand this function consumes. Taking the floor from the same call keeps
   % it to one evaluation and makes it the floor the decision actually used.
   [wetflag, f_res] = icemodel.column.vapor_exchange_is_wet( ...
      f_ice, f_liq, f_res_por);

   % Retain initial top-layer fractions.
   f_liq_top = f_liq;
   f_ice_top = f_ice;

   % Initialize potential deposition that cannot be satisfied by the cv budget.
   d_sbl_err = 0;

   % If a liquid film is present, partition liquid vapor exchange first.
   % Otherwise route the exchange straight to the ice phase, so dry or cold
   % deposition forms ice instead of liquid water that is not there.
   if wetflag

      if d_pevp < 0 % evaporation

         if f_liq_top < f_res % f_liq_top - f_res < 0
            % only residual water exists, send d_pevp to sublimation

            if debug == true
               fprintf('f_liq(1) < f_res, d_pevp sent to sublimation\n')
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
         [f_ice, d_sbl_err] = sublimation(d_pevp, f_ice, f_liq, f_ice_min);

      elseif d_pevp > 0 % condensation

         % Compute the capacity for condensation.
         d_aevp_max = icemodel.column.max_liquid_fraction_change( ...
            f_ice_top, f_liq_top);

         if d_pevp <= d_aevp_max
            % all condensation stored
            d_aevp = d_pevp; % d_pevp = 0

            f_liq = f_liq + d_aevp;

         else
            % some condensation stored, some converts to runoff
            d_aevp = d_aevp_max;

            f_liq = f_liq + d_aevp;

            % Condensation beyond what this cell's pore space can hold. The
            % excess cannot be stored, so the caller routes it: to runoff
            % from the top cell (d_rof sends it to diagnose_column_runoff),
            % or to the unapplied accounting from an interior cell. The
            % excess is real water, so dropping it would break the budget.
            d_pevp = d_pevp - d_aevp_max;
            d_rof = d_rof + d_pevp;

            if debug == true
               fprintf( ...
                  'condensation exceeds porosity: %.6f\n', d_pevp)
            end
         end
      end

   else % dry/cold ice: sublimation or direct deposition to ice
      [f_ice, d_sbl_err] = sublimation(d_pevp, f_ice, f_liq, f_ice_min);
   end

   if debug == true && d_sbl_err > 0
      fprintf('rejected deposition: %.6f\n', d_sbl_err)
   end

   % The realized exchange, from the state the limits produced plus the
   % overflow. The liquid change is already on the liquid basis; the ice
   % change carries the same mass at a different density, so the density
   % ratio converts it. Overflow condensate crossed the surface before it
   % ran off, so it counts: only the amounts d_sbl_err rejected never
   % crossed.
   d_applied = (f_liq - f_liq_top) ...
      + (f_ice - f_ice_top) * ro_ice / ro_liq + d_rof;
end

%% Vapor exchange with the ice phase
function [f_ice, d_sbl_err] = sublimation(d_pevp, f_ice, f_liq, f_ice_min)
   %SUBLIMATION Convert potential evaporation to sublimation.
   %
   % icemodel.column.potential_sublimation holds the conversion and
   % derives it from Qe, Lv, Ls, ro_liq, and ro_ice.
   %
   % In icemodel, ro_air_Lv is set to ro_air * Lv or ro_air * Ls depending on
   % liqflag, so Qe is already computed wrt to them. That way evap/subl are
   % computed using the same formula: e = Qe / ro_air_Lv * dt / dz.
   %
   % This conversion conserves heat when evaporation alone cannot satisfy the
   % surface latent heat flux Qe, that is, when all available water evaporates.
   % The evaporation branch then sends the remaining d_pevp to this function,
   % which applies it to sublimation of ice.

   % The control volume sends excess condensation to runoff. It cannot do the
   % same with excess deposition, so this function reports that amount in
   % d_sbl_err. Initialize the value to 0.
   d_sbl_err = 0;

   % Early return if there's no energy for sublimation.
   if d_pevp == 0
      return
   end

   % Convert potential evap to potential subl in ice frac-equivalent thickness.
   % icemodel.column.merge_thin_layers predicts against this same conversion,
   % so both call one helper and cannot drift apart.
   d_psbl = icemodel.column.potential_sublimation(d_pevp);

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

   % Layer combination is based on f_ice, so requiring f_ice_top < 0 should
   % suffice (rather than <f_ice_min). If f_ice(1) + d_psbl < 0, it will
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
