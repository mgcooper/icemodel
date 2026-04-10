function Qli = incoming_longwave_radiation(Tair, ea_atm)
   %incoming_longwave_radiation Estimate downwelling longwave radiation.
   %
   %  Qli = icemodel.surface.incoming_longwave_radiation(Tair, ea_atm)
   %
   % Inputs:
   %  Tair   - air temperature [K]
   %  ea_atm - atmospheric vapor pressure [Pa]
   %
   % Output:
   %  Qli - incoming longwave radiation [W m^-2]
   %
   % This is a legacy fallback used when forcing data do not provide
   % downwelling longwave radiation directly. The formula is kept explicit
   % here so initialize_surface_forcings can derive `lwd` from raw station
   % data without reintroducing root-level all-caps helpers.
   %
   %#codegen

   persistent SB
   if isempty(SB)
      % Cache the Stefan-Boltzmann constant across repeated forcing loads.
      SB = icemodel.physicalConstant('SB');
   end

   % Estimate effective sky emissivity from vapor pressure and air
   % temperature using the legacy forcing fallback relationship.
   emiss_cloud = 1.08 .* ...
      (1.0 - exp(-(0.01 .* ea_atm) .^ (Tair ./ 2016.0)));

   % Convert emissivity to incoming longwave radiation.
   Qli = emiss_cloud .* SB .* Tair .^ 4;
end
