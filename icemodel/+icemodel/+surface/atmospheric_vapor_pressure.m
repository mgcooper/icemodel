function ea = atmospheric_vapor_pressure(Tair, rh, liqflag)
   %ATMOSPHERIC_VAPOR_PRESSURE Relative humidity to atmospheric vapor pressure.
   %
   %  ea = icemodel.surface.atmospheric_vapor_pressure(Tair, rh, liqflag)
   %
   %  Computes:
   %     ea = es(Tair, liqflag) * rh / 100   [Pa]
   %
   %  where RH is supplied in percent and ea is returned in Pa.
   %
   %  Phase convention: standard synoptic RH observations follow the WMO
   %  convention of referencing saturation to liquid water for T > 0 °C
   %  and to ice for T < 0 °C. The correct phase flag for atmospheric use
   %  is therefore `Tair > Tf`, NOT the surface liqflag. Callers should pass
   %  `Tair > Tf` (or a pre-computed logical of the same shape) rather than
   %  propagating a surface-state flag to this function.
   %
   %  Supports scalar and vector inputs; liqflag may be a scalar or a
   %  logical array matching the shape of Tair.
   %
   % See also: icemodel.vapor.saturation_vapor_pressure,
   %           icemodel.surface.initialize_surface_state
   %
   %#codegen

   persistent Tf
   if isempty(Tf)
      Tf = icemodel.physicalConstant('Tf');
   end

   % Use tair > Tf as the phase selector: WMO RH data are referenced to liquid
   % water for T>0°C and to ice for T<0°C, independent of surface phase state.
   if nargin < 3
      liqflag = Tair > Tf;
   end

   ea = icemodel.vapor.saturation_vapor_pressure(Tair, liqflag) .* rh ./ 100;
end
