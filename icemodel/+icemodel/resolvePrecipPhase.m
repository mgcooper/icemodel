function [rainf, snowf] = resolvePrecipPhase(ppt, tair, rainf_source, ...
      snowf_source, phase_source)
   %RESOLVEPRECIPPHASE Select the runtime rain/snow split (POLICY A10/D-18).
   %
   %  [rainf, snowf] = icemodel.resolvePrecipPhase(ppt, tair, ...
   %     rainf_source, snowf_source, phase_source)
   %
   % Inputs
   %  ppt          - canonical total precipitation rate [m s-1]
   %  tair         - air temperature [K]
   %  rainf_source - the met product's own liquid component [m s-1]
   %  snowf_source - the met product's own solid component [m s-1]
   %  phase_source - runtime phase-source option (opts.precip_phase_source
   %     from icemodel.setopts):
   %     'source'    the product's split exactly as shipped (e.g. MAR's
   %                 energy-balance split); absent or missing components
   %                 stay missing rather than being fabricated.
   %     'threshold' repartition PPT by air temperature via
   %                 icemodel.forcing.reconstruct.partitionPrecipitation;
   %                 the transition temperature defaults from
   %                 icemodel.forcing.reconstruct.setopts (single source).
   %
   % Both modes enforce the POLICY A10 validity contract: every finite value
   % is nonnegative, a finite phase cannot exceed a finite total, and every
   % complete split sums to the total.
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

   % One shared sample axis keeps the selection honest across inputs.
   if numel(tair) ~= numel(ppt) || numel(rainf_source) ~= numel(ppt) ...
         || numel(snowf_source) ~= numel(ppt)
      error('icemodel:resolvePrecipPhase:sizeMismatch', ...
         ['ppt, tair, rainf_source, and snowf_source must share one ' ...
         'sample axis']);
   end

   switch lower(char(phase_source))
      case 'source'
         % A10: the product's split is data; expose it bit-identically.
         rainf = rainf_source;
         snowf = snowf_source;
      case 'threshold'
         % A10/D-18: runtime threshold partition of the canonical total;
         % the kwargs default inside partitionPrecipitation supplies the
         % single-source transition temperature.
         [rainf, snowf] = ...
            icemodel.forcing.reconstruct.partitionPrecipitation(ppt, tair);
      otherwise
         error('icemodel:resolvePrecipPhase:invalidPhaseSource', ...
            'precip_phase_source must be ''source'' or ''threshold'': %s', ...
            char(phase_source));
   end

   % POLICY A10 validity uses the same shared helper as reconstruction and
   % artifact verification. Honest missingness remains untouched.
   violates = ~icemodel.forcing.helpers.precipitationValidity( ...
      ppt, rainf, snowf);
   if any(violates)
      error('icemodel:resolvePrecipPhase:inconsistentSplit', ...
         ['finite precipitation values must be nonnegative, each phase ' ...
         'must not exceed ppt, and complete splits must sum to ppt: ' ...
         '%d violating samples'], nnz(violates));
   end
end
