function [rainf, snowf] = resolvePrecipPhase(ppt, tair, rainf_source, ...
      snowf_source, phase_source)
   %RESOLVEPRECIPPHASE Select the runtime rainf/snowf split.
   %
   %  [rainf, snowf] = icemodel.resolvePrecipPhase(ppt, tair, ...
   %     rainf_source, snowf_source, phase_source)
   %
   % Inputs
   %  ppt          - total precipitation rate [m s-1]
   %  tair         - air temperature [K]
   %  rainf_source - the met product's own rainfall [m s-1]
   %  snowf_source - the met product's own snowfall [m s-1]
   %  phase_source - runtime phase-source option (opts.precip_phase_source
   %     from icemodel.setopts):
   %     'source'    the met product's split as provided.
   %     'threshold' repartition PPT by air temperature with
   %                 icemodel.forcing.reconstruct.partitionPrecipitation. The
   %                 transition temperature defaults from
   %                 icemodel.forcing.reconstruct.setopts.
   %
   % Both modes enforce every finite value is nonnegative, a finite phase cannot
   % exceed a finite total, and every complete split sums to the total.
   %
   % See also: icemodel.surface.initialize_surface_forcings
   %           icemodel.forcing.reconstruct.partitionPrecipitation

   arguments
      ppt (:, 1) double
      tair (:, 1) double
      rainf_source (:, 1) double
      snowf_source (:, 1) double
      phase_source {mustBeTextScalar}
   end

   % Every input must use the same sample axis.
   if numel(tair) ~= numel(ppt) || numel(rainf_source) ~= numel(ppt) ...
         || numel(snowf_source) ~= numel(ppt)
      error('icemodel:resolvePrecipPhase:sizeMismatch', ...
         ['ppt, tair, rainf_source, and snowf_source must share one ' ...
         'sample axis']);
   end

   switch lower(char(phase_source))
      case 'source'
         % Use the met product's values.
         rainf = rainf_source;
         snowf = snowf_source;
      case 'threshold'
         % Repartition total ppt.
         [rainf, snowf] = ...
            icemodel.forcing.reconstruct.partitionPrecipitation(ppt, tair);
      otherwise
         error('icemodel:resolvePrecipPhase:invalidPhaseSource', ...
            'precip_phase_source must be ''source'' or ''threshold'': %s', ...
            char(phase_source));
   end

   % Check for bad data. The check leaves missing values in place.
   violates = ~icemodel.forcing.helpers.precipitationValidity( ...
      ppt, rainf, snowf);
   if any(violates)
      error('icemodel:resolvePrecipPhase:inconsistentSplit', ...
         ['finite precipitation values must be nonnegative, each phase ' ...
         'must not exceed ppt, and complete splits must sum to ppt: ' ...
         '%d violating samples'], nnz(violates));
   end
end
