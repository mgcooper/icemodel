function [liqflag, ro_sfc, hv_atm, H_e, f_res_por] = ...
      update_surface_state(f_ice_1, f_liq_1, ro_atm, De_e, ...
      snow_depth, opts)
   %UPDATE_SURFACE_STATE Update the surface state at substep entry.
   %
   %  [liqflag, ro_sfc, hv_atm, H_e, f_res_por] = ...
   %     icemodel.surface.update_surface_state( ...
   %     f_ice_1, f_liq_1, ro_atm, De_e, snow_depth, opts)
   %
   %  This function computes surface state variables that depend on both
   %  the current column state and the current forcing step. It's called at the
   %  start of each substep so the coupler sees state consistent with the
   %  current metstep and updated column state.
   %
   %  This function is the substep-entry complement to initialize_surface_state,
   %  which precomputes surface state values that depend only on forcing data
   %  (ea_atm, ro_atm, cv_atm, H_h, De_e, etc.) and thus can be precomputed once
   %  at model startup.
   %
   %  Inputs:
   %    f_ice_1    - top-node ice fraction (scalar) [-]
   %    f_liq_1    - top-node liquid fraction (scalar) [-]
   %    ro_atm     - moist-air density at current metstep [kg m-3]
   %    De_e       - latent exchange coefficient at current metstep
   %                 [m s-1 Pa-1]
   %    snow_depth - surface snow depth [m] (placeholder until the
   %                 snow model is added, at which point snow_depth
   %                 will evolve within substeps)
   %    opts       - model options struct (uses f_res_pore_snow,
   %                 f_res_pore_ice)
   %
   %  Outputs:
   %    liqflag   - surface phase flag; true when top node contains
   %                enough liquid to use the liquid-phase saturation
   %                vapor-pressure curve
   %    ro_sfc    - surface bulk density [kg m-3]
   %    hv_atm    - volumetric latent enthalpy [J m-3]
   %                = ro_atm * Lv (liqflag) or ro_atm * Ls
   %    H_e       - latent heat transport prefactor [W m-2 Pa-1]
   %                = hv_atm * De_e
   %    f_res_por - residual pore-water fraction [-]; selects between
   %                snow and ice values based on snow_depth
   %
   % See also:
   %   icemodel.surface.initialize_surface_state,
   %   icemodel.timestepping.acceptsubstep
   %
   %#codegen

   persistent Lv Ls f_liq_phase_switch_threshold
   if isempty(Lv)
      [Lv, Ls] = icemodel.physicalConstant('Lv', 'Ls');
      f_liq_phase_switch_threshold = icemodel.parameterLookup( ...
         'f_liq_phase_switch_threshold');
   end

   % Surface phase flag from current top-node liquid fraction.
   liqflag = f_liq_1 > f_liq_phase_switch_threshold;

   % Surface bulk density from current top-node column state.
   ro_sfc = icemodel.surface.surface_bulk_density(f_ice_1, f_liq_1);

   % Volumetric latent enthalpy from current metstep density.
   if liqflag
      hv_atm = ro_atm * Lv;
   else
      hv_atm = ro_atm * Ls;
   end

   % Latent heat transport coefficient [W m-2 Pa-1]
   H_e = hv_atm * De_e;

   % Residual pore-water fraction [-]
   if snow_depth > 0
      f_res_por = opts.f_res_pore_snow;
   else
      f_res_por = opts.f_res_pore_ice;
   end
end
