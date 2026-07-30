function opts = setopts(kwargs)
   %SETOPTS Central options for the reconstruction pipeline.
   %
   %  opts = icemodel.forcing.reconstruct.setopts()
   %  opts = icemodel.forcing.reconstruct.setopts(blend_hours=3)
   %
   % Role
   %  Single source of every scalar knob and channel list the
   %  reconstruction pipeline consumes. Engine functions default their
   %  own name-value arguments from this function and the production
   %  driver passes one opts struct down, so each parameter is defined
   %  exactly once. Per-channel physical bounds, admission bias caps,
   %  and gap-duration bucket edges stay with their dedicated
   %  single-source functions (physicalBounds, admissionGate,
   %  bucketEdges); this function owns the scalar knobs, the channel
   %  namelists, and the proxy-source catalog derived from the repo's
   %  canonical alias map.
   %
   % Name-value (all optional; defaults are the approved policy values)
   %  required_channels : channels the ready_icemodel verdict must
   %     complete (POLICY A5): the seven-channel icemodel set. Snowfall
   %     input (ppt OR snowf) is graded separately by the ready_snowmodel
   %     verdict, and swu is derived (B10), never required.
   %  core_channels : native state channels the low-coverage confidence
   %     advisory grades; they never replace the readiness verdict.
   %  plan_channels : channels the selection experiment plans methods for
   %     (precipitation is adopted per POLICY A10, never planned; swu
   %     follows albedo*swd downstream). Every required non-precipitation
   %     channel must appear here so the product carries its provenance.
   %  interp_channels : channels eligible for tier-1 bounded
   %     interpolation; unsupported channels are rejected at construction.
   %  proxy_sources : short RCM family labels, in adoption-preference
   %     order. Only MAR and MERRA are complete-meteorology fallbacks;
   %     storage directories and provenance identities are derived, never
   %     restated (see proxyCatalog).
   %  cap_hours : tier-1 interior interpolation cap (wall-clock hours),
   %     at most the adopted six-hour policy ceiling.
   %  cap_hours_by_channel : per-channel tier-1 cap overrides (POLICY
   %     B3/D-21); fields name interp channels and replace cap_hours for
   %     that channel only, within the channel ceiling. Default SWD=9 h
   %     from the observed-only D-39 holdout; every other channel uses
   %     cap_hours.
   %  jump_factor : seam-limit multiplier on the seasonal step scale.
   %  blend_hours : seam-blend taper window (hours).
   %  seam_qa_percentile : observed-step percentile used to identify
   %     anomalous post-final SWD method boundaries (POLICY D-32).
   %  seam_qa_min_reference_steps : minimum observed steps required in a
   %     season x solar-elevation band before the screen falls back to the
   %     season-wide reference.
   %  seam_qa_max_passes : maximum one-posting synthetic-side repair
   %     iterations. Two passes resolve cascaded KANL boundaries without
   %     the extra distortion caused by a two-posting window.
   %  toa_dark_wm2 : top-of-atmosphere irradiance below which a sample
   %     counts as dark for the swd CSI mask and the darkness zero-fill.
   %  max_donor_distance_km, max_donor_elev_diff_m : donor geometry gate.
   %  min_overlap_hours : concurrent finite overlap a donor/proxy fit
   %     requires.
   %  max_lag_hours, min_lag_gain : donor lag search half-width and the
   %     correlation improvement required to adopt a nonzero lag.
   %  max_extrapolation_fraction : how far beyond the fitted donor range
   %     a transfer may be applied, as a fraction of that range.
   %  rmse_improvement : fractional held-out RMSE improvement over the
   %     policy baseline (persistence for the <=6 h bucket, station
   %     climatology otherwise, POLICY B8) that admission requires.
   %  min_variability_ratio, max_variability_ratio : admitted prediction
   %     standard deviation divided by held-out truth standard deviation
   %     must stay inside this interval when truth has nonzero spread.
   %  min_coverage : admission usefulness floor on reconstructed
   %     fraction of drawn validation samples.
   %  max_validation_duration_factor : a real gap may not exceed this
   %     multiple of the longest held-out gap that admitted its method.
   %  z0_bulk : aerodynamic roughness threshold below which boom-height
   %     geometry is not runtime-valid.
   %  min_native_core_coverage : low-confidence advisory floor on native
   %     core-channel coverage per year.
   %  lapse_rate : fallback temperature lapse rate (K per m) for donor
   %     elevation adjustment when an overlap fit is underdetermined.
   %  elevation_threshold_m : elevation difference above which donor
   %     elevation adjustment applies.
   %  tair_for_pressure : representative air temperature (K) for the
   %     barometric pressure elevation adjustment.
   %  min_season_samples : overlap samples a season needs for its own
   %     proxy calibration before falling back to the annual fit.
    %  rain_snow_transition_temperature_k : air-temperature threshold the
    %     RUNTIME 'threshold' phase option uses to partition the product's
    %     total precipitation (POLICY A10/D-18; reconstruction itself
    %     never partitions).
    %  native_winter_albedo, native_winter_months : the native PROMICE
   %     builder's winter-albedo stamp (from icemodel.parameterLookup);
   %     samples carrying exactly this constant in these months are
   %     legacy fills, not observations, and re-enter the engine as
   %     missing so methods fill them with honest provenance.
   %  plan_n_gaps : synthetic validation gaps the planner draws per
   %     stratum (duration bucket x season) and split.
   %  knot_candidates : monotone piecewise knot counts the donor-transfer
   %     hyperparameter search tries (0 = pure linear).
   %  max_donors : donor stations the planner fits per channel.
   %  climatology_window_days : day-of-year half-window for the
   %     climatology lookup pool.
   %  climatology_min_support : samples the climatology pool needs before
   %     it produces an estimate.
   %  synthetic_n_gaps : synthetic gaps the missingness sampler draws.
   %  synthetic_context_hours : observed context required around each
   %     synthetic gap.
   %  selection_fraction : fraction of observed years assigned to the
   %     selection split (the rest are evaluation years).
   %  last_resort_proxies : adopt aligned proxy values for residual
   %     missing required-channel samples after composition (bounds
   %     checked, stamped with the proxy code, audited as last resort).
   %  use_modis_albedo : adopt staged GEUS MODIS daily albedo ahead of
   %     the RCM proxies in albedo's last-resort order (POLICY A11/B12).
   %  seed : deterministic seed for the selection experiment.
   %
   % Returns
   %  opts : struct of the fields above plus proxy_catalog, a struct
   %     array with one row per proxy source — label (short family
   %     label), storage (versioned artifact directory token), and
   %     code_name (provenance-registry field).
   %
   % See also: icemodel.forcing.reconstruct.fillPromiceStation,
   %  icemodel.forcing.reconstruct.reconstructSeries,
   %  icemodel.forcing.reconstruct.admissionGate,
   %  icemodel.verification.namelists.rcmProductIds

   arguments
      kwargs.required_channels (1, :) string = ...
         icemodel.forcing.reconstruct.icemodelRequiredChannels()
      kwargs.core_channels (1, :) string = ["tair", "rh", "psfc"]
      kwargs.plan_channels (1, :) string = ...
         ["tair", "rh", "wspd", "psfc", "swd", "swu", "lwd", "albedo"]
      kwargs.interp_channels (1, :) string = interpolableChannels()
      kwargs.proxy_sources (1, :) string = supportedProxySources()
      kwargs.cap_hours (1, 1) double {mustBePositive} = ...
         defaultInterpolationCap()
      % Per-channel tier-1 cap overrides (POLICY B3/D-21): a struct whose
      % fields name interp channels and whose values replace cap_hours for
      % that channel only. Defaults empty = every channel uses cap_hours.
      % The default SWD override is supported by the D-39 held-out
      % evidence. Any future extension still needs its own evidence.
      kwargs.cap_hours_by_channel (1, 1) struct = struct()
      kwargs.jump_factor (1, 1) double {mustBePositive} = 3
      kwargs.blend_hours (1, 1) double {mustBeNonnegative} = 6
      kwargs.seam_qa_percentile (1, 1) double ...
         {mustBeInRange(kwargs.seam_qa_percentile, 50, 100)} = 99.9
      kwargs.seam_qa_min_reference_steps (1, 1) double ...
         {mustBeInteger, mustBePositive} = 100
      kwargs.seam_qa_max_passes (1, 1) double ...
         {mustBeInteger, mustBeNonnegative} = 2
      kwargs.toa_dark_wm2 (1, 1) double {mustBePositive} = 10
      kwargs.max_donor_distance_km (1, 1) double {mustBePositive} = 60
      kwargs.max_donor_elev_diff_m (1, 1) double {mustBePositive} = 600
      kwargs.min_overlap_hours (1, 1) double {mustBePositive} = 8760
      kwargs.max_lag_hours (1, 1) double {mustBeNonnegative} = 18
      kwargs.min_lag_gain (1, 1) double {mustBeNonnegative} = 0.02
      kwargs.max_extrapolation_fraction (1, 1) double ...
         {mustBeNonnegative} = 0.10
      kwargs.rmse_improvement (1, 1) double {mustBeNonnegative} = 0.10
      kwargs.min_variability_ratio (1, 1) double {mustBePositive} = 0.67
      kwargs.max_variability_ratio (1, 1) double {mustBePositive} = 1.50
      kwargs.min_coverage (1, 1) double = 0.10
      kwargs.max_validation_duration_factor (1, 1) double ...
         {mustBePositive} = 1
      kwargs.z0_bulk (1, 1) double {mustBePositive} = ...
         icemodel.parameterLookup('thf_z0_bulk')
      kwargs.min_native_core_coverage (1, 1) double = 0.30
      kwargs.lapse_rate (1, 1) double = -0.0060
      kwargs.elevation_threshold_m (1, 1) double {mustBePositive} = 100
       kwargs.tair_for_pressure (1, 1) double {mustBePositive} = 255
       kwargs.min_season_samples (1, 1) double {mustBePositive} = 300
       kwargs.rain_snow_transition_temperature_k (1, 1) double ...
          {mustBePositive} = icemodel.physicalConstant('Tf')
       kwargs.native_winter_albedo (1, 1) double = ...
         icemodel.parameterLookup('promice_winter_albedo')
      kwargs.native_winter_months (1, :) double = ...
         icemodel.parameterLookup('promice_winter_albedo_months')
      kwargs.plan_n_gaps (1, 1) double {mustBeInteger, mustBePositive} = 12
      kwargs.knot_candidates (1, :) double = [0 6]
      kwargs.max_donors (1, 1) double {mustBeInteger, mustBePositive} = 3
      kwargs.climatology_window_days (1, 1) double {mustBePositive} = 7
      kwargs.climatology_min_support (1, 1) double {mustBePositive} = 5
      kwargs.synthetic_n_gaps (1, 1) double ...
         {mustBeInteger, mustBePositive} = 25
      kwargs.synthetic_context_hours (1, 1) double {mustBePositive} = 24
      kwargs.selection_fraction (1, 1) double = 0.7
      kwargs.last_resort_proxies (1, 1) logical = true
      kwargs.use_modis_albedo (1, 1) logical = true
      kwargs.seed (1, 1) double {mustBeInteger, mustBeNonnegative} = ...
         20260724
   end

   if kwargs.min_variability_ratio > kwargs.max_variability_ratio
      error('icemodel:reconstruct:setopts:variabilityRange', ...
         'min_variability_ratio must not exceed max_variability_ratio');
   end
   coverage = [kwargs.min_coverage, kwargs.min_native_core_coverage];
   if any(coverage < 0 | coverage > 1)
      error('icemodel:reconstruct:setopts:coverageRange', ...
         'coverage thresholds must lie in the inclusive interval [0, 1]');
   end
   caps = icemodel.forcing.reconstruct.interpolationCapHours();
   if kwargs.cap_hours > caps.default
      error('icemodel:reconstruct:setopts:capHours', ...
         'cap_hours must not exceed the adopted six-hour policy ceiling');
   end
   % Per-channel overrides obey the same ceiling and must name channels
   % that tier 1 can actually interpolate.
   for f = string(fieldnames(kwargs.cap_hours_by_channel)).'
      if ~ismember(f, kwargs.interp_channels)
         error('icemodel:reconstruct:setopts:capChannel', ...
            'cap override names a non-interpolable channel: %s', f);
      end
      v = kwargs.cap_hours_by_channel.(f);
      ceiling = caps.default;
      if isfield(caps, f)
         ceiling = caps.(f);
      end
      if ~(isscalar(v) && isfinite(v) && v > 0 && v <= ceiling)
         error('icemodel:reconstruct:setopts:capHours', ...
            'cap override for %s must lie in (0, %g]', f, ceiling);
      end
   end
   precipitation = icemodel.forcing.helpers.precipitationVariables();
   planned_precipitation = intersect(kwargs.plan_channels, precipitation, ...
      'stable');
   if ~isempty(planned_precipitation)
      error('icemodel:reconstruct:setopts:plannedPrecipitation', ...
         ['precipitation is adopted and partitioned by its dedicated ' ...
         'policy, never planned statistically: %s'], ...
         strjoin(planned_precipitation, ', '));
   end
   unsupported_interpolation = setdiff(kwargs.interp_channels, ...
      interpolableChannels(), 'stable');
   if ~isempty(unsupported_interpolation)
      error('icemodel:reconstruct:setopts:unsupportedInterpolationChannel', ...
         'channels have no tier-1 interpolation policy: %s', ...
         strjoin(unsupported_interpolation, ', '));
   end
   % Materialize evidenced defaults only after the caller's interpolation
   % schema is valid. A deliberately narrower schema must not inherit an
   % irrelevant override, and an explicit caller value wins.
   defaults = defaultInterpolationOverrides();
   for f = string(fieldnames(defaults)).'
      if ismember(f, kwargs.interp_channels) ...
            && ~isfield(kwargs.cap_hours_by_channel, f)
         kwargs.cap_hours_by_channel.(f) = defaults.(f);
      end
   end
   proxy_sources = lower(kwargs.proxy_sources);
   unsupported_proxies = setdiff(proxy_sources, ...
      supportedProxySources(), 'stable');
   if ~isempty(unsupported_proxies)
      error('icemodel:reconstruct:setopts:unsupportedProxySource', ...
         'sources have no complete-meteorology fallback policy: %s', ...
         strjoin(unsupported_proxies, ', '));
   end
   unplanned_required = setdiff(kwargs.required_channels, ...
      [kwargs.plan_channels, precipitation], 'stable');
   if ~isempty(unplanned_required)
      error('icemodel:reconstruct:setopts:requiredChannelNotPlanned', ...
         ['required reconstruction channels must also be planned; ' ...
         'only adopted total precipitation is exempt: %s'], ...
         strjoin(unplanned_required, ', '));
   end
   opts = kwargs;
   opts.proxy_sources = proxy_sources;
   % The catalog is derived, never restated: storage tokens come from the
   % repo's canonical alias map and provenance identities from the
   % published code registry.
   opts.proxy_catalog = proxyCatalog(opts.proxy_sources);
end

function cap = defaultInterpolationCap()
   %DEFAULTINTERPOLATIONCAP Return the ordinary tier-1 ceiling.
   caps = icemodel.forcing.reconstruct.interpolationCapHours();
   cap = caps.default;
end

function overrides = defaultInterpolationOverrides()
   %DEFAULTINTERPOLATIONOVERRIDES Return evidence-backed channel overrides.
   overrides = rmfield( ...
      icemodel.forcing.reconstruct.interpolationCapHours(), 'default');
end

function channels = interpolableChannels()
   %INTERPOLABLECHANNELS Channels with an approved tier-1 implementation.
   channels = ["tair", "rh", "wspd", "psfc", "lwd", "swd", "albedo"];
end

function labels = supportedProxySources()
   %SUPPORTEDPROXYSOURCES Complete-meteorology fallback source labels.
   labels = ["mar", "merra"];
end

function catalog = proxyCatalog(labels)
   %PROXYCATALOG Resolve proxy labels to storage and provenance identity.
   storage = icemodel.verification.namelists.rcmProductIds(labels);
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   slots = cell(numel(labels), 1);
   for k = 1:numel(labels)
      % Provenance identity: the family label itself when the registry
      % carries it, else the storage token (merra resolves to merra2).
      if isfield(codes, labels(k))
         code_name = labels(k);
      elseif isfield(codes, storage(k))
         code_name = storage(k);
      else
         error('icemodel:reconstruct:setopts:unknownProxySource', ...
            'no provenance code registered for proxy source %s', labels(k));
      end
      slots{k} = struct('label', labels(k), 'storage', storage(k), ...
         'code_name', code_name);
   end
   catalog = vertcat(slots{:});
end
