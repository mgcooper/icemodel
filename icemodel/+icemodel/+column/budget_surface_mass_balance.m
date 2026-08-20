function [T, f_ice, f_liq, d_liq, d_evp, d_rof, d_applied, budget] = ...
      budget_surface_mass_balance(T, f_ice, f_liq, xf_liq, d_pevp, d_liq, ...
      d_evp, d_rof, f_res_por, f_ice_min, budget, dz)
   %BUDGET_SURFACE_MASS_BALANCE Budget surface mass-balance increments.
   %
   % budget_surface_mass_balance updates the cumulative liquid-water and
   % vapor-driven mass-change increments over the current full step. It uses
   % the phase state that the latest substep solve already updated, applies
   % the surface vapor exchange, and accumulates the surface vapor budget,
   % so the driver sees one call for the whole surface mass balance.
   %
   % Inputs
   %   T          - Column temperature state [K].
   %   f_ice      - Column ice fraction [-].
   %   f_liq      - Column liquid-water fraction [-].
   %   xf_liq     - Prior liquid-water fraction carried into this substep [-].
   %   d_pevp     - Potential surface vapor-driven liquid-fraction change for
   %                the top cell [-], an energy demand the SEB fixed.
   %   d_liq      - Accumulated liquid-water fraction change over the step [-].
   %   d_evp      - Accumulated top-node liquid-storage change from surface
   %                vapor exchange over the forcing step [-].
   %   d_rof      - Accumulated condensation overflow over the step [-].
   %   f_res_por  - Residual liquid-water fraction per pore volume [-].
   %   f_ice_min  - Minimum retained surface ice fraction before remeshing [-].
   %   budget     - Forcing-step budget with the phase channels and post-phase
   %                baselines already accumulated (see
   %                icemodel.column.accumulate_phase_budget).
   %   dz         - Control-volume thickness [m].
   %
   % Outputs
   %   T          - Updated column temperature state [K].
   %   f_ice      - Updated column ice fraction [-].
   %   f_liq      - Updated column liquid-water fraction [-].
   %   d_liq      - Updated cumulative liquid-water fraction change [-].
   %   d_evp      - Updated top-node liquid-storage change from surface vapor
   %                exchange [-].
   %   d_rof      - Updated condensation overflow in liquid-water fraction [-].
   %   d_applied  - The exchange the top cell realized [-], signed like
   %                d_pevp, as a liquid-water volume fraction, after every
   %                limit. Surface-budget diagnostics consume realized
   %                exchange because rejected demand never crossed the
   %                surface.
   %   budget     - Budget with the surface vapor increments accumulated and
   %                the post-exchange baselines in budget.substep, the
   %                reference state for icemodel.column.couple_vapor_step.
   %
   % Notes
   %   This routine does not merge thin layers. Call
   %   `icemodel.column.merge_thin_layers` after this bookkeeping step to keep
   %   remeshing explicit in the main model flow.
   %
   % See also: icemodel.surface.apply_surface_vapor_exchange,
   %           icemodel.column.accumulate_vapor_budget,
   %           icemodel.column.merge_thin_layers
   %
   %#codegen

   % Compute surface mass-balance term updates.
   %
   % Positive values of d_liq indicate increasing water fraction:
   % d_liq > 0 = melt.
   % d_evp > 0 = condensation.
   % d_rof > 0 = condensation which exceeds top layer porosity.
   % d_liq < 0 = freeze.
   % d_evp < 0 = evaporation.

   % Compute delta f_liq.
   %
   % The only process that affects f_liq between d_liq assignments is phase
   % change (icemodel.column.solve_column_enthalpy).
   % icemodel.column.couple_vapor_step runs after this function, so interior
   % transport never lands in d_liq and is never read as melt or refreezing.
   d_liq = d_liq + f_liq - xf_liq;

   % Reset past values for budgeting evap/condensation.
   xf_liq = f_liq;

   % Apply the surface vapor exchange to the top cell. d_sbl_err returns one
   % entry per cell; only the first entry can be nonzero here.
   %
   % Interior transport does not enter this function.
   % icemodel.column.couple_vapor_step converts the accepted solve flux and
   % calls the same icemodel.column.apply_vapor_transfer state mutator. Its
   % mass-conserving redistribution budget remains separate from this
   % energy-conserving surface budget.
   [f_ice, f_liq, d_rof, d_sbl_err, d_applied] = ...
      icemodel.surface.apply_surface_vapor_exchange( ...
      f_ice, f_liq, d_rof, d_pevp, f_ice_min, f_res_por);

   % Budget evap / cond. This differences f_liq, so it only captures liquid
   % vapor exchange. Sublimation and deposition change f_ice, so they are not
   % included in d_evp (use the budget's vapor_solid channel for those).
   d_evp = d_evp + f_liq - xf_liq;

   % Close this substep's surface vapor budget while the exchange terms are
   % in scope, so the driver never carries the shortfall record.
   budget = icemodel.column.accumulate_vapor_budget( ...
      budget, T, f_ice, f_liq, dz, d_pevp, d_rof, d_sbl_err);
end
