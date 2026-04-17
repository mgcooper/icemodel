function snow_depth = resolve_forcing_snow_depth(forcing_snow_depth, metstep, ...
      use_forcing_snow_depth_for_thf)
   %RESOLVE_FORCING_SNOW_DEPTH Resolve one scalar snow-depth value for THF.
   %
   %  snow_depth = icemodel.surface.resolve_forcing_snow_depth( ...
   %     forcing_snow_depth, metstep, use_forcing_snow_depth_for_thf)
   %
   % The current models have no explicit snow thermodynamics yet, so the
   % forcing-based snow-depth hook is optional and THF-only. When the hook is
   % disabled or the requested met step has no finite forcing value, the
   % returned scalar stays on the current bare-ice fallback (0.0).
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
