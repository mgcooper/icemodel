function [Qle, dQle_dTsfc] = outgoing_longwave_radiation(T_sfc, Qli)
   %OUTGOING_LONGWAVE_RADIATION Outgoing longwave radiation and T_sfc derivative.
   %
   %  Qle = icemodel.surface.outgoing_longwave_radiation(T_sfc)
   %  [Qle, dQle_dTsfc] = icemodel.surface.outgoing_longwave_radiation(T_sfc)
   %
   % In accordance with Kirchoff's Law, the upwelling longwave flux contains an
   % emitted component, following the Stefan-Boltzmann equation, and a reflected
   % component, which is proportional to the downwelling flux:
   %
   %   Qle = emiss * SB * T_sfc^4 + (1-emiss) * Qli
   %
   % The sign convention is positive toward the surface
   %
   %#codegen

   persistent emiss SB
   if isempty(emiss)
      SB = icemodel.physicalConstant('SB');
      emiss = icemodel.parameterLookup('emiss');
   end

   Qle = - SB * emiss * T_sfc^4 - (1 - emiss) * Qli;

   if nargout > 1
      dQle_dTsfc = -4.0 * emiss * SB * T_sfc ^ 3;
   end
end
