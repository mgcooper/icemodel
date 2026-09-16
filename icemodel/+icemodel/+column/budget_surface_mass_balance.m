function [T_ice, f_ice, f_liq, d_liq, d_evp, d_rof, d_vap_liq, d_vap_ice, ...
      d_vap, budget] = budget_surface_mass_balance(T_ice, f_ice, f_liq, ...
      xf_liq, d_pevp, d_liq, d_evp, d_rof, d_vap_liq, d_vap_ice, ...
      f_res_por, f_ice_min, budget, dz)
   %BUDGET_SURFACE_MASS_BALANCE Apply and budget the surface mass balance.
   %
   %  [T_ice, f_ice, f_liq, d_liq, d_evp, d_rof, d_vap_liq, d_vap_ice, ...
   %     d_vap, budget] = icemodel.column.budget_surface_mass_balance( ...
   %     T_ice, f_ice, f_liq, xf_liq, d_pevp, d_liq, d_evp, d_rof, ...
   %     d_vap_liq, d_vap_ice, f_res_por, f_ice_min, budget, dz)
   %
   % Record melt and freeze after the enthalpy solve, apply surface vapor
   % exchange, and add the applied vapor changes to the forcing-step budget.
   %
   % Sign conventions:
   %   d_liq > 0 melt, d_liq < 0 freeze.
   %   d_evp > 0 condensation, d_evp < 0 evaporation.
   %   d_rof > 0 condensation above the top-cell porosity (runs off).
   %
   % Inputs
   %   T_ice      - Column temperature state [K].
   %   f_ice      - Column ice fraction [-].
   %   f_liq      - Column liquid-water fraction [-].
   %   xf_liq     - Column liquid-water fraction for the previous substep [-].
   %   d_pevp     - Surface vapor energy demand on the Lv basis [-], from
   %                icemodel.surface.potential_surface_vapor_demand.
   %   d_liq      - Accumulated phase-change increments per cell [-].
   %   d_evp      - Accumulated liquid-basis vapor exchange increments per cell [-].
   %   d_rof      - Accumulated condensation overflow (runoff) [-].
   %   d_vap_liq  - Accumulated vapor-driven liquid increments per cell [-]
   %                (exchange plus transport; written out as ice2.df_vap_liq).
   %   d_vap_ice  - Accumulated vapor-driven ice increments per cell [-],
   %                ice fraction basis (written out as ice2.df_vap_ice).
   %   f_res_por  - Residual liquid-water fraction per pore volume [-].
   %   f_ice_min  - Minimum retained ice fraction before remeshing [-].
   %   budget     - Forcing-step budget (icemodel.column.initialize_budget_state).
   %   dz         - Control-volume thickness [m].
   %
   % Outputs
   %   T_ice      - Updated column temperature [K].
   %   f_ice      - Updated column solid ice phase fraction [-].
   %   f_liq      - Updated column liquid phase fraction [-].
   %   d_liq      - Updated melt/freeze phase-change accumulator [-].
   %   d_evp      - Updated liquid-basis vapor exchange accumulator [-].
   %   d_rof      - Updated condensation overflow (runoff) [-].
   %   d_vap      - Realized exchange this substep, liquid-water basis,
   %                signed like d_pevp [-]. Rejected demand never crossed
   %                the surface, so diagnostics consume this value.
   %   budget     - Budget with the exchange increments added.
   %
   % See also: icemodel.surface.apply_surface_vapor_exchange,
   %           icemodel.column.accumulate_vapor_exchange,
   %           icemodel.column.merge_thin_layers
   %
   %#codegen

   % Record solid/liquid phase change after the enthalpy solve, before vapor
   % exchange. The only process that changes f_liq between the last xf_liq
   % checkpoint and this assignment is melt/freeze phase change.
   d_liq = d_liq + f_liq - xf_liq;

   % Apply surface vapor exchange and return its applied phase increments.
   % apply_surface_vapor_exchange returns the applied per-cell increments
   % directly, unlike the state-differencing used for d_liq.
   [f_ice, f_liq, d_rof, d_vap_liq_x, d_vap_ice_x, d_vap] = ...
      icemodel.surface.apply_surface_vapor_exchange( ...
      f_ice, f_liq, d_rof, d_pevp, dz, f_ice_min, f_res_por);

   % Add evaporation/condensation to the liquid vapor diagnostic.
   % Sublimation/deposition modify f_ice and are tracked in the budget's
   % vapor_solid value.
   d_evp = d_evp + d_vap_liq_x;

   % Add surface exchange to the per-cell vapor diagnostics. These accumulators
   % track every vapor-driven storage change per cell; the surface exchange is
   % budgeted here, subsurface transport is budgeted in couple_vapor_step.
   d_vap_liq = d_vap_liq + d_vap_liq_x;
   d_vap_ice = d_vap_ice + d_vap_ice_x;

   % Add this substep's exchange to the forcing-step budget.
   budget = icemodel.column.accumulate_vapor_exchange( ...
      budget, d_pevp, d_vap_liq_x, d_vap_ice_x, d_rof, dz);

   % Below here (keep for reference, don't delete):
   % f_liq = ... some modification to f_liq
   % d_XXX = f_liq - xf_liq - d_evp
   %
   % Equivalently, re-assign xf_liq first:
   % xf_liq = f_liq;
   % f_liq = ... some modification to f_liq
   % d_XXX = f_liq - xf_liq
end
