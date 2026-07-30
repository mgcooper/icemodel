function step_opts = step_observation_heights(opts, metstep)
   %STEP_OBSERVATION_HEIGHTS Select scalar observation heights for one step.
   %
   % PROMICE loadmet calls carry time-varying boom heights in the three
   % observation-height options. Other forcing families retain scalar heights.
   % This boundary keeps iterative flux kernels scalar for either contract.
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
