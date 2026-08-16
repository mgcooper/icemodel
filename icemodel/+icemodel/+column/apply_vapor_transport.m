function [f_ice, f_liq, d_sbl_err] = apply_vapor_transport( ...
      f_ice, f_liq, d_vap, f_ice_min, f_res_por)
   %APPLY_VAPOR_TRANSPORT Apply interior vapor transport as mass.
   %
   %  [f_ice, f_liq, d_sbl_err] = ...
   %     icemodel.column.apply_vapor_transport( ...
   %     f_ice, f_liq, d_vap, f_ice_min, f_res_por)
   %
   % D_VAP is the mass each cell gained from its neighbours over one
   % substep, as a liquid-water volume fraction, from
   % icemodel.column.couple_vapor_transport, where sum(d_vap .* dz) is
   % zero.
   %
   % Mass, not energy. This is the difference from
   % icemodel.surface.apply_surface_vapor_exchange, which spends an energy
   % demand the surface energy balance fixed and lets the mass follow. Fick's
   % law fixes the mass instead, so this applies exactly that mass and lets
   % the energy follow. A cell therefore never spends part of its exchange at
   % one latent heat and the rest at another: it takes the mass it received.
   %
   % Each cell exchanges with the phase it holds, decided by
   % icemodel.column.vapor_exchange_is_wet, the same function the surface
   % applier uses.
   %
   % The same three limits apply, and they clamp and record rather than
   % reject: deposition is capped by the air space, sublimation is floored at
   % f_ice_min, and condensation into a wet cell is capped by
   % icemodel.column.max_liquid_fraction_change. Whatever a cell cannot take
   % is recorded in D_SBL_ERR, on the ice-phase energy basis, so no mass
   % leaves without a record.
   %
   % The closure identity does not balance when a clamp binds. The shortfall
   % enters the ledger's unapplied term, while
   % icemodel.column.accumulate_redistribution_budget records only the
   % storage that survived the clamp, so the identity carries a residual the
   % size of the shortfall. Bead icemodel-55x carries the fix.
   %
   % Inputs
   %   f_ice     - Ice fraction [-] (JJ x 1).
   %   f_liq     - Liquid-water fraction [-] (JJ x 1).
   %   d_vap    - Vapor gained per cell as a liquid-water volume fraction
   %              [-] (JJ x 1). Positive is gain. A dry cell converts it to
   %              the ice basis before applying it.
   %   f_ice_min - Minimum retained ice fraction [-].
   %   f_res_por - Residual liquid-water fraction per pore volume [-].
   %
   % Outputs
   %   f_ice, f_liq - Updated fractions.
   %   d_sbl_err    - Unapplied amount per cell, in ice-fraction units
   %                  [-] (JJ x 1).
   %
   % See also: icemodel.column.couple_vapor_transport,
   %  icemodel.surface.apply_surface_vapor_exchange,
   %  icemodel.column.vapor_exchange_is_wet
   %
   %#codegen

   persistent ro_ice ro_liq
   if isempty(ro_ice)
      [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   end

   n_cells = numel(f_ice);
   d_sbl_err = zeros(n_cells, 1);
   % One call gives both the decision and the floor it used. The column then
   % evaluates residual_water_fraction once per cell rather than twice, and
   % the floor cannot disagree with the decision.
   [wet, f_res] = icemodel.column.vapor_exchange_is_wet( ...
      f_ice, f_liq, f_res_por);

   for j = 1:n_cells
      if d_vap(j) == 0
         continue
      end

      if wet(j)
         % A liquid film takes the exchange. The mass is already a liquid
         % volume fraction, so it applies directly.
         if d_vap(j) > 0
            % Condensation, capped by the pore space this cell can hold.
            capacity = icemodel.column.max_liquid_fraction_change( ...
               f_ice(j), f_liq(j));
            applied = min(d_vap(j), max(capacity, 0));
            f_liq(j) = f_liq(j) + applied;
            d_sbl_err(j) = icemodel.column.potential_sublimation( ...
               d_vap(j) - applied);
         else
            % Evaporation, floored at the residual the film retains.
            available = max(f_liq(j) - f_res(j), 0);
            applied = min(abs(d_vap(j)), available);
            f_liq(j) = f_liq(j) - applied;
            d_sbl_err(j) = icemodel.column.potential_sublimation( ...
               d_vap(j) + applied);
         end
      else
         % No film, so the exchange is with ice. Convert the liquid-water
         % volume fraction to the ice volume fraction of the same mass.
         d_ice = d_vap(j) * ro_liq / ro_ice;
         if d_ice > 0
            % Deposition, capped by the air space.
            f_air = 1.0 - f_ice(j) - f_liq(j);
            applied = min(d_ice, max(f_air, 0));
            f_ice(j) = f_ice(j) + applied;
            d_sbl_err(j) = d_ice - applied;
         else
            % Sublimation, floored at the retained ice minimum.
            available = max(f_ice(j) - f_ice_min, 0);
            applied = min(abs(d_ice), available);
            f_ice(j) = f_ice(j) - applied;
            d_sbl_err(j) = d_ice + applied;
         end
      end
   end
end
