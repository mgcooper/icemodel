function step_opts = step_observation_heights(opts, metstep)
   %STEP_OBSERVATION_HEIGHTS Select scalar observation heights for one step.
   %
   %  step_opts = icemodel.surface.step_observation_heights(opts, metstep)
   %
   % Select element METSTEP from vector-valued PROMICE observation heights.
   % Keep scalar observation heights unchanged for forcing data that has one
   % height for the full run.
   %
   % See also: icemodel.loadmet,
   %  icemodel.surface.diagnose_turbulent_heat_fluxes
   %
   %#codegen

   step_opts = opts;
   names = {'z_tair', 'z_wind', 'z_relh'};
   for k = 1:numel(names)
      if ~isfield(opts, names{k})
         continue
      end
      values = opts.(names{k});
      if ~isscalar(values)
         step_opts.(names{k}) = values(metstep);
      end
   end
end
