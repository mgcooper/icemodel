function [Ts, T, f_ice, f_liq, dt_sum, dt_new, ...
      liqflag, ro_air_Lv, ro_sfc, f_res_por] = updatesubstep(Ts, T, f_ice, ...
      f_liq, dt_FULL_STEP, dt_sum, dt_new, TINY, snow_depth, opts)
   % Checkpoint the accepted state and allocate time within the full step.
   %
   %  [Ts, T, f_ice, f_liq, dt_sum, dt_new] = updatesubstep(...)
   %  [..., liqflag, ro_air_Lv] = updatesubstep(...)
   %  [..., liqflag, ro_air_Lv, ro_sfc] = updatesubstep(...)
   %  [..., liqflag, ro_air_Lv, ro_sfc, f_res_por] = ...
   %     updatesubstep(..., snow_depth, opts)
   %
   % State outputs (7-10) are derived from the accepted substep state and are
   % used as the initial values for the next substep:
   %   liqflag   - surface phase flag (true when top layer has enough liquid)
   %   ro_air_Lv - ro_air * latent heat for surface turbulent flux [J m-3]
   %   ro_sfc    - surface bulk density [kg m-3]
   %   f_res_por - residual pore water fraction; depends on snow_depth and opts
   %
   %#codegen

   persistent f_liq_phase_switch_threshold roLv roLs
   if isempty(f_liq_phase_switch_threshold)
      f_liq_phase_switch_threshold = icemodel.parameterLookup( ...
         'f_liq_phase_switch_threshold');
      [roLv, roLs] = icemodel.physicalConstant('roLv', 'roLs');
   end

   % Allocate this substep to the timestep
   dt_sum = dt_sum + dt_new;

   % Adjust dt to exactly complete the full step without going over.
   % The first condition is true if the full step is incomplete, the second
   % is true if the next substep will exceed the full step.
   if (dt_FULL_STEP - dt_sum) > TINY && (dt_sum + dt_new - dt_FULL_STEP) > TINY
      dt_new = dt_FULL_STEP - dt_sum; % dt_new = max(dt - dt_sum, dt_min);
   end

   if nargout > 6

      % Top node contains enough liquid water to use the liquid-phase curve.
      liqflag = f_liq(1) > f_liq_phase_switch_threshold;

      % ro_air_Lv is for surface equations - subsurface node "wetness" varies
      if liqflag
         ro_air_Lv = roLv;  % ro_air * Lv
      else
         ro_air_Lv = roLs;  % ro_air * Ls
      end
   end

   if nargout > 8
      ro_sfc = icemodel.surface.surface_bulk_density(f_ice(1), f_liq(1));
   end

   if nargout > 9
      if snow_depth > 0
         f_res_por = opts.f_res_pore_snow;
      else
         f_res_por = opts.f_res_pore_ice;
      end
   end
end
