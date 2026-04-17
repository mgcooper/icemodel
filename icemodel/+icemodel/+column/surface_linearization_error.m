function err = surface_linearization_error(T_sfc, T_top, Fc, Fp, a1)
   %SURFACE_LINEARIZATION_ERROR Diagnose the Robin surface linearization error.
   %
   %  err = icemodel.column.surface_linearization_error(T_sfc, T_top, Fc, Fp, a1)
   %
   % Returns the mismatch, expressed as an equivalent temperature error [K],
   % between the linearized surface-flux boundary condition and the conductive
   % temperature jump at the top control-volume interface.
   %
   %#codegen

   err = -(Fc + Fp * T_sfc) / a1 - (T_top - T_sfc);
end
