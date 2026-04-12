function [T, f_ice, f_liq, d_liq, d_evp, x_err] = ...
      budget_surface_mass_balance(T, f_ice, f_liq, xf_liq, d_pevp, d_liq, ...
      d_evp, f_liq_res, f_ice_min)
   %BUDGET_SURFACE_MASS_BALANCE Budget surface mass-balance increments.
   %
   % budget_surface_mass_balance updates the cumulative liquid-water and
   % vapor-driven mass-change increments over the current full step using the
   % already-updated phase state from the latest substep solve.
   %
   % Inputs
   %   T          - Column temperature state [K].
   %   f_ice      - Column ice fraction [-].
   %   f_liq      - Column liquid-water fraction [-].
   %   xf_liq     - Prior liquid-water fraction carried into this substep [-].
   %   d_pevp     - Potential surface vapor-driven liquid-fraction change [-].
   %   d_liq      - Accumulated liquid-water fraction change over the step [-].
   %   d_evp      - Accumulated vapor-driven fraction change over the step [-].
   %   f_liq_res  - Residual liquid-water fraction per pore volume [-].
   %   f_ice_min  - Minimum retained surface ice fraction before remeshing [-].
   %
   % Outputs
   %   T          - Updated column temperature state [K].
   %   f_ice      - Updated column ice fraction [-].
   %   f_liq      - Updated column liquid-water fraction [-].
   %   d_liq      - Updated cumulative liquid-water fraction change [-].
   %   d_evp      - Updated cumulative vapor-driven fraction change [-].
   %   x_err      - Unsatisfied vapor-driven ice change for debugging only [-].
   %
   % Notes
   %   This routine does not merge thin layers. Call
   %   `icemodel.column.merge_thin_layers` after this bookkeeping step to keep
   %   remeshing explicit in the main model flow.
   %
   % See also: icemodel.surface.apply_surface_vapor_mass_change,
   %  icemodel.column.merge_thin_layers
   %
   %#codegen

   % Compute the latest surface mass-balance contributions from the updated
   % phase state and the prescribed surface vapor tendency.
   [T, d_liq, d_evp, ~, f_ice, f_liq, x_err] = surface_mass_balance_terms( ...
      T, d_liq, d_evp, 0, f_ice, f_liq, xf_liq, d_pevp, f_liq_res, ...
      f_ice_min);
end

function [T, d_liq, d_evp, d_con, f_ice, f_liq, x_err] = ...
      surface_mass_balance_terms(T, d_liq, d_evp, d_con, f_ice, f_liq, ...
      xf_liq, d_pevp, f_liq_res, f_ice_min)
   %SURFACE_MASS_BALANCE_TERMS Compute surface mass-balance term updates.
   %
   % Positive values of d_liq indicate increasing water fraction:
   % d_liq > 0 = melt.
   % d_evp > 0 = condensation.
   % d_con > 0 = condensation which exceeds top layer porosity.
   % d_liq < 0 = freeze.
   % d_evp < 0 = evaporation.

   persistent TL f_ell_min
   if isempty(TL)
      [TL, ~, f_ell_min] = icemodel.column.meltzone_bounds();
   end

   % Compute delta f_liq. Note, the only process which affects f_liq between
   % d_liq assignments here is phase change
   % (`icemodel.column.solve_column_enthalpy`). d_liq < 0 means
   % refreezing, d_liq > 0 melt.
   d_liq = d_liq + f_liq - xf_liq;

   % Reset past values for budgeting evap/condensation.
   xf_liq = f_liq;

   % Update the residual unfrozen water fraction defined by the phase fraction
   % characteristic function at the lower temperature boundary. This ensures
   % f_liq is not reduced below residual water for melting nodes.
   if T(1) > TL
      f_liq_min = icemodel.water_fraction(f_ice(1), f_liq(1)) * f_ell_min;
      f_liq_res = max(f_liq_res, f_liq_min / (1 - f_ice(1)));
   end

   % Wetflag controls liquid-film mass partitioning only. Treat only liquid
   % above the residual-water floor as a mobile surface film. Keep this
   % separate from the SEB surface saturation flag used for vapor pressure.
   f_res_top = f_liq_res * (1.0 - f_ice(1));
   wetflag = f_liq(1) > f_res_top;

   % Compute evaporation/condensation/sublimation.
   [f_ice, f_liq, d_con, x_err] = ...
      icemodel.surface.apply_surface_vapor_mass_change( ...
      f_ice, f_liq, d_con, d_pevp, wetflag, f_ice_min, f_liq_res);

   % Budget evap / subl.
   d_evp = d_evp + f_liq - xf_liq;

   % Below here:
   % f_liq = ... some modification to f_liq
   % d_XXX = f_liq - xf_liq - d_evp

   % Equivalently, re-assign xf_liq first:
   % xf_liq = f_liq;
   % f_liq = ... some modification to f_liq
   % d_XXX = f_liq - xf_liq
end
