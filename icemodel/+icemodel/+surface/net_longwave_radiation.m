function Qln = net_longwave_radiation(T_sfc, Qli)
   %NET_LONGWAVE_RADIATION Compute net longwave radiation at the surface.
   %
   %  Qln = icemodel.surface.net_longwave_radiation(T_sfc, Qli)
   %
   % The sign convention is positive toward the surface:
   %
   %   Qln = emiss * Qli - emiss * SB * T_sfc^4
   %
   %#codegen

   persistent emiss SB
   if isempty(emiss)
      [SB] = icemodel.physicalConstant('SB');
      emiss = icemodel.parameterLookup('emiss');
   end

   Qln = emiss * (Qli - SB * T_sfc ^ 4);
end
