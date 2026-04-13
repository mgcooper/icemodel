function [Ts, T, f_ice, f_liq, dt_sum, dt_new, liqflag, ro_air_Lv, varargout] ...
      = updatesubstep(Ts, T, f_ice, f_liq, dt_FULL_STEP, dt_sum, dt_new, TINY)
   % Checkpoint the accepted state and allocate time within the full step.
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

   if nargout > 5

      % Top node contains enough liquid water to use the liquid-phase curve.
      liqflag = f_liq(1) > f_liq_phase_switch_threshold;

      % ro_air_Lv is for surface equations - subsurface node "wetness" varies
      if liqflag
         ro_air_Lv = roLv;  % ro_air * Lv
      else
         ro_air_Lv = roLs;  % ro_air * Ls
      end
   end

   % Legacy option to return updated state.
   if nargout > 7
      % Bulk density and specific heat capacity.
      ro_sno = icemodel.column.bulk_density(f_ice, f_liq);
      cp_sno = icemodel.column.bulk_specific_heat_capacity( ...
         f_ice, f_liq, ro_sno);
   end

   % To activate this, need k_eff and dt. But
   if nargout > 9
      % Diffusion length scale
      % zD = sqrt(k_eff(1) * dt / (ro_sno(1) * cp_sno(1)));
      % if zD > dz(1)
      %  % placeholder
      % end
   end

   switch nargout
      case 8
         varargout{1} = ro_sno;
      case 9
         varargout{1} = ro_sno;
         varargout{2} = cp_sno;
   end
end
