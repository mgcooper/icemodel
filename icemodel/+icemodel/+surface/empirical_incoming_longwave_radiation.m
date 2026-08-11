function Qli = empirical_incoming_longwave_radiation(Tair, ea_atm, emiss)
   %EMPIRICAL_INCOMING_LONGWAVE_RADIATION Estimate downwelling longwave.
   %
   %  Qli = icemodel.surface.empirical_incoming_longwave_radiation(Tair, ea_atm)
   %
   % Inputs:
   %  Tair   - air temperature [K]
   %  ea_atm - atmospheric vapor pressure [Pa]
   %
   % Output:
   %  Qli - incoming longwave radiation [W m^-2]
   %
   % This is a legacy fallback for forcing data that has no downwelling
   % longwave radiation. initialize_surface_forcings uses it to derive `lwd`
   % from raw station data. The canonical SEB longwave accounting stays in
   % `icemodel.surface.net_longwave_radiation` and
   % `icemodel.surface.outgoing_longwave_radiation`.
   %
   %#codegen

   persistent SB
   if isempty(SB)
      SB = icemodel.physicalConstant('SB');
   end

   % Estimate effective sky emissivity from vapor pressure and air
   % temperature using the legacy forcing fallback relationship.
   if nargin < 3
      emiss = 1.08 .* (1.0 - exp(-(0.01 .* ea_atm) .^ (Tair ./ 2016.0)));
   end

   % Convert emissivity to incoming longwave radiation.
   Qli = emiss .* SB .* Tair .^ 4;
end
