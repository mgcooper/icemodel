function err = subsurface_linearization_error(T, T_old, f_ice, f_liq, ...
      f_liq_old, drovdT, dHdT, Sc, dz, dt, a2, Fc, Fp, a1)
   %SUBSURFACE_LINEARIZATION_ERROR Diagnose the top-node enthalpy linearization error.
   %
   %  err = icemodel.column.subsurface_linearization_error(...)
   %
   % Returns the top-node mismatch, expressed as an equivalent temperature
   % error [K], between the accepted enthalpy update and the linearized
   % top-node energy balance used by the subsurface solve.
   %
   %#codegen

   persistent Lf Ls Lv ro_liq
   if isempty(Lf)
      [Lf, Ls, Lv, ro_liq] = icemodel.physicalConstant( ...
         'Lf', 'Ls', 'Lv', 'ro_liq');
   end

   persistent f_liq_phase_switch_threshold
   if isempty(f_liq_phase_switch_threshold)
      f_liq_phase_switch_threshold = icemodel.parameterLookup( ...
         'f_liq_phase_switch_threshold');
   end

   if f_liq(1) > f_liq_phase_switch_threshold
      Lv_top = Lv;
   else
      Lv_top = Ls;
   end

   err = (dt / dz(1) ...
      * (a2 * (T(2) - T(1)) ...
      + Fc + Fp * (Fc + a1 * T(1)) / (a1 - Fp) ...
      + Sc(1) * dz(1)) ...
      - ro_liq * Lf * (f_liq(1) - f_liq_old(1))) ...
      / (dHdT(1) + Lv_top * drovdT(1) * (1.0 - f_ice(1) - f_liq(1))) ...
      - (T(1) - T_old(1));
end
