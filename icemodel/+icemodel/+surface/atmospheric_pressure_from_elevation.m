function Pa = atmospheric_pressure_from_elevation(topo)
   %atmospheric_pressure_from_elevation Estimate pressure from elevation.
   %
   %  Pa = icemodel.surface.atmospheric_pressure_from_elevation(topo)
   %
   % Input:
   %  topo - station elevation above sea level [m]
   %
   % Output:
   %  Pa - atmospheric pressure [Pa]
   %
   % This barometric fallback serves station datasets that give elevation but
   % not surface pressure. initialize_surface_forcings uses it when the
   % forcing pressure must be derived.
   %
   %#codegen

   persistent one_atmos scale_ht
   if isempty(one_atmos)
      % These are the same constants as the legacy helper. State them here
      % so the namespaced fallback documents them.
      one_atmos = 101300.0;
      scale_ht = 8000.0;
   end

   % Apply an exponential scale-height pressure profile.
   Pa = one_atmos .* exp(-topo ./ scale_ht);
end
