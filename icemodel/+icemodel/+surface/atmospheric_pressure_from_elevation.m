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
   % This simple barometric fallback is retained for station datasets that
   % provide elevation but not surface pressure. It is intended for use in
   % initialize_surface_forcings when forcing pressure must be derived.
   %
   %#codegen

   persistent one_atmos scale_ht
   if isempty(one_atmos)
      % Use the same constants as the legacy helper, but document them
      % explicitly at the namespaced fallback boundary.
      one_atmos = 101300.0;
      scale_ht = 8000.0;
   end

   % Apply an exponential scale-height pressure profile.
   Pa = one_atmos .* exp(-topo ./ scale_ht);
end
