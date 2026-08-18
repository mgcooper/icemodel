function [T, f_ice, f_liq, d_liq, d_evp, d_rof, d_sbl_err, d_applied] = ...
      budget_surface_mass_balance(T, f_ice, f_liq, xf_liq, d_pevp, d_liq, ...
      d_evp, d_rof, f_res_por, f_ice_min, d_sbl_err)
   %BUDGET_SURFACE_MASS_BALANCE Budget surface mass-balance increments.
   %
   % budget_surface_mass_balance updates the cumulative liquid-water and
   % vapor-driven mass-change increments over the current full step. It uses
   % the phase state that the latest substep solve already updated.
   %
   % Inputs
   %   T          - Column temperature state [K].
   %   f_ice      - Column ice fraction [-].
   %   f_liq      - Column liquid-water fraction [-].
   %   xf_liq     - Prior liquid-water fraction carried into this substep [-].
   %   d_pevp     - Potential surface vapor-driven liquid-fraction change for
   %                the top cell [-], an energy demand the SEB fixed.
   %   d_liq      - Accumulated liquid-water fraction change over the step [-].
   %   d_evp      - Accumulated vapor-driven fraction change over the step [-].
   %   d_rof      - Accumulated condensation overflow over the step [-].
   %   f_res_por  - Residual liquid-water fraction per pore volume [-].
   %   f_ice_min  - Minimum retained surface ice fraction before remeshing [-].
   %   d_sbl_err  - (optional) Running per-cell unapplied vapor record [-].
   %                The surface term is added to whatever arrives rather
   %                than replacing it. Omit it and the record starts at
   %                zero, which is the production convention: the interior
   %                transport runs after this function and keeps its
   %                shortfall in its own redistribution accounting.
   %
   % Outputs
   %   T          - Updated column temperature state [K].
   %   f_ice      - Updated column ice fraction [-].
   %   f_liq      - Updated column liquid-water fraction [-].
   %   d_liq      - Updated cumulative liquid-water fraction change [-].
   %   d_evp      - Updated cumulative vapor-driven fraction change [-].
   %   d_rof      - Updated condensation overflow in liquid-water fraction [-].
   %   d_sbl_err  - Signed unapplied vapor-driven ice-fraction change per
   %                cell [-], the incoming record plus the surface term.
   %   d_applied  - The exchange the top cell realized [-], signed like
   %                d_pevp, as a liquid-water volume fraction, after every
   %                limit. Grain growth consumes this rather than the
   %                demand, because rejected demand never crossed the
   %                surface.
   %
   % Notes
   %   This routine does not merge thin layers. Call
   %   `icemodel.column.merge_thin_layers` after this bookkeeping step to keep
   %   remeshing explicit in the main model flow.
   %
   % See also: icemodel.surface.apply_surface_vapor_exchange,
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

   % A caller from before the coupled path threads no record in. Start it at
   % zero so the ten-argument convention still works.
   if nargin < 11
      d_sbl_err = 0.0;
   end

   % Compute delta f_liq.
   %
   % The only process that affects f_liq between d_liq assignments is phase
   % change (icemodel.column.solve_column_enthalpy). That holds in coupled
   % mode too: icemodel.column.couple_vapor_step runs after this function
   % and after the vapor budget close, so interior transport never lands in
   % d_liq and is never read as melt or refreezing.
   d_liq = d_liq + f_liq - xf_liq;

   % Reset past values for budgeting evap/condensation.
   xf_liq = f_liq;

   % Apply the surface vapor exchange. It reaches the top cell, and
   % d_sbl_err returns one entry per cell so the ledger integrates one shape
   % whichever path filled it. Coupled interior transport does not come
   % through here: it has its own pair,
   % icemodel.column.couple_vapor_transport and
   % icemodel.column.apply_vapor_transport, because it conserves
   % mass where this conserves energy.
   [f_ice, f_liq, d_rof, d_sbl_err_sfc, d_applied] = ...
      icemodel.surface.apply_surface_vapor_exchange( ...
      f_ice, f_liq, d_rof, d_pevp, f_ice_min, f_res_por);

   % Add rather than assign, so a caller that threads a running record in
   % keeps its accumulation.
   %
   % Keep the caller's shape. The surface applier reaches the top cell alone,
   % so its record is the first entry and every other entry is zero. A caller
   % that threaded a scalar in, which is the ten-argument convention and the
   % default path, gets a scalar back holding that top-cell term. Adding the
   % vector to a scalar would expand it and hand that caller a shape its
   % established form never returned.
   if isscalar(d_sbl_err)
      d_sbl_err = d_sbl_err + d_sbl_err_sfc(1);
   else
      d_sbl_err = d_sbl_err + d_sbl_err_sfc;
   end

   % Budget evap / cond. This differences f_liq, so it only captures liquid
   % vapor exchange. Sublimation and deposition change f_ice, so they are not
   % included in d_evp (use the ledger's vapor_solid channel for those).
   d_evp = d_evp + f_liq - xf_liq;

   % Below here:
   % f_liq = ... some modification to f_liq
   % d_XXX = f_liq - xf_liq - d_evp

   % Equivalently, re-assign xf_liq first:
   % xf_liq = f_liq;
   % f_liq = ... some modification to f_liq
   % d_XXX = f_liq - xf_liq
end
