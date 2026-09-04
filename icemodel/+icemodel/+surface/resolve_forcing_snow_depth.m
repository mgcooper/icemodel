function snow_depth = resolve_forcing_snow_depth(forcing_snow_depth, metstep, ...
      use_forcing_snow_depth_for_thf)
   %RESOLVE_FORCING_SNOW_DEPTH Resolve scalar snow-depth for the THF scheme.
   %
   %  snow_depth = icemodel.surface.resolve_forcing_snow_depth( ...
   %     forcing_snow_depth, metstep, use_forcing_snow_depth_for_thf)
   %
   % Returns the forcing snow depth for one timestep. The current model has no
   % explicit snow model, so this hook is optional and currently only used to
   % select roughness values for the THF scheme. When the hook is disabled or
   % the requested forcing step has no finite value, this function returns the
   % bare-ice fallback 0.0.
   %
   % See also: icemodel, skinmodel,
   %  icemodel.surface.initialize_surface_forcings
   %
   %#codegen

   snow_depth = 0.0;
   if ~use_forcing_snow_depth_for_thf || isempty(forcing_snow_depth)
      return
   end

   if metstep < 1 || metstep > numel(forcing_snow_depth)
      return
   end

   value = forcing_snow_depth(metstep);
   if isfinite(value)
      snow_depth = max(0.0, value);
   end
end
