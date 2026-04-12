function [Qln, dQln_dTsfc] = net_longwave_radiation(T_sfc, Qli)
   %NET_LONGWAVE_RADIATION Net surface longwave radiation and T_sfc derivative.
   %
   %  Qln = icemodel.surface.net_longwave_radiation(T_sfc, Qli)
   %  [Qln, dQln_dTsfc] = icemodel.surface.net_longwave_radiation(T_sfc, Qli)
   %
   % In accordance with Kirchoff's Law, the upwelling longwave flux contains an
   % emitted component following the Stefan-Boltzmann equation, and a reflected
   % component, which is proportional to the downwelling flux:
   %
   %   Qle = - (emiss * SB * T_sfc^4 + (1-emiss) * Qli)
   %   Qln = Qli + Qle
   %       = Qli - (emiss * SB * T_sfc^4 + (1-emiss) * Qli)
   %       = -emiss * SB * T_sfc^4 + emiss * Qli
   %       = emiss * (Qli - SB * T_sfc^4)
   %
   % The sign convention is positive toward the surface.
   %
   %#codegen

   persistent emiss SB
   if isempty(emiss)
      SB = icemodel.physicalConstant('SB');
      emiss = icemodel.parameterLookup('emiss');
   end

   Qln = emiss * (Qli - SB * T_sfc ^ 4);

   if nargout > 1
      dQln_dTsfc = -4.0 * emiss * SB * T_sfc ^ 3;
   end
end
