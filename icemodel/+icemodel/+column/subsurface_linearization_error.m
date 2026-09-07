function err = subsurface_linearization_error(T, T_old, f_ice, f_liq, ...
      f_liq_old, drovdT, dHdT, Sc, q_deferred_faces, dz, dt, a2, ...
      Fc, Fp, a1)
   %SUBSURFACE_LINEARIZATION_ERROR Diagnose the top-node enthalpy
   % linearization error.
   %
   %  err = icemodel.column.subsurface_linearization_error(...)
   %
   % Returns the top-node mismatch as an equivalent temperature error [K]
   % between the accepted enthalpy update and the linearized top-node energy
   % balance used by the subsurface solve. The deferred face contribution uses
   % the downward-positive in-minus-out convergence from
   % assemble_enthalpy_system.
   %
   %#codegen

   persistent Lf ro_liq
   if isempty(Lf)
      [Lf, ro_liq] = icemodel.physicalConstant('Lf', 'ro_liq');
   end

   % Phase-aware latent heat for the top node (scalar).
   Lv_top = icemodel.vapor.latent_enthalpy_switch(f_liq(1));

   % A one-node column has no interior south face; its bottom face is closed.
   q_south = 0.0;
   if numel(T) > 1
      q_south = a2 * (T(2) - T(1));
   end

   err = (dt / dz(1) ...
      * (q_south ...
      + Fc + Fp * (Fc + a1 * T(1)) / (a1 - Fp) ...
      + Sc(1) * dz(1) ...
      + q_deferred_faces(1) - q_deferred_faces(2)) ...
      - ro_liq * Lf * (f_liq(1) - f_liq_old(1))) ...
      / (dHdT(1) + Lv_top * drovdT(1) * (1.0 - f_ice(1) - f_liq(1))) ...
      - (T(1) - T_old(1));
end
