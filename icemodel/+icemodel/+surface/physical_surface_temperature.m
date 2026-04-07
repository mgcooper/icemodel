function T_sfc = physical_surface_temperature(T_sfc)
   %PHYSICAL_SURFACE_TEMPERATURE Cap surface temperature at the melting point.
   %
   %  T_sfc = icemodel.surface.physical_surface_temperature(T_sfc)
   %
   % Returns min(T_sfc, Tf). The solver may produce temperatures above Tf
   % as an internal iterate; this function enforces the physical constraint
   % before any surface flux evaluation.
   %
   % See also: MFENERGY, icemodel.surface.solve_surface_energy_balance
   %
   %#codegen

   persistent Tf
   if isempty(Tf)
      Tf = icemodel.physicalConstant('Tf');
   end

   if T_sfc > Tf
      T_sfc = Tf;
   end
end
