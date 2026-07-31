function plan = stationMethodPlan(target, donors, proxies, kwargs)
   %STATIONMETHODPLAN Select and fit admitted fill methods for one station.
   %
   %  plan = icemodel.forcing.reconstruct.stationMethodPlan( ...
   %     target, donors, proxies, channels=["tair", "rh"], ...
   %     split_manifest="...", seed=20260724)
   %
   % Role
   %  The per-station selection experiment (gap-fill policy validation
   %  design): for every requested channel it fits the candidate methods on
   %  SELECTION years, grades each on held-out blocked synthetic draws
   %  through the harness, admits per stratum via admissionGate against
   %  persistence for <=6 h or climatology otherwise, ranks admitted methods
   %  by common-support fractional improvement, and
   %  re-grades the ranked winners on EVALUATION-year draws for the
   %  reported final numbers. Everything is role-typed: the same function
   %  plans a PROMICE target from K-transect donors or a K-transect target
   %  from PROMICE donors (the role-reversal acceptance test).
   %
   % Inputs
   %  target : struct — series (timetable), station (string), location
   %     (struct lat_wgs84, lon_wgs84, elev_m).
   %  donors : struct array (possibly empty) — series, station, family,
   %     location, and optional observed_mask (logical vector or timetable
   %     by channel; filled samples are never donors). Source sensor-height
   %     and AWS-generation metadata in series.Properties.UserData is carried
   %     into fitted transfer provenance.
   %  proxies : struct array (possibly empty) — series (model met on the
   %     target axis), name ("mar"|"merra2"), code_name (provenanceCodes
   %     field name).
   %
   % Name-value
   %  channels : channels to plan (default the policy's required set).
   %  seed : deterministic seed for split and draws (required).
   %  split_manifest : persisted split path ("" disables persistence).
   %  n_gaps : synthetic gaps per exact stratum (duration bucket x season).
   %  knot_candidates : spline hyperparameter set.
   %  max_donors : distinct donors retained per channel after held-out
   %     validation ranking.
    %  max_donor_distance_km, max_donor_elev_diff_m : geometry gate.
    %  cap_hours, cap_hours_by_channel : wall-clock tier-1 interpolation
    %     caps used by the census.
    %  Remaining name-values forward the identically named policy controls
    %     from setopts to split, fitting, calibration, validation, and gates.
   %  Defaults come from the central
   %  icemodel.forcing.reconstruct.setopts contract.
   %
   % Returns
   %  plan : struct — station, split, per-channel struct array (channel,
   %     methods: ordered admitted methods with fitted parameters, codes,
   %     per-stratum admissions, selection and evaluation metrics), and
   %     the census the readiness ledger reuses.
   %
   % See also: icemodel.forcing.reconstruct.reconstructSeries,
   %  icemodel.forcing.reconstruct.admissionGate

   arguments
      target (1, 1) struct
      donors (1, :) struct
      proxies (1, :) struct
       kwargs.channels (1, :) string = ...
          icemodel.forcing.reconstruct.setopts().plan_channels
       kwargs.core_channels (1, :) string = ...
          icemodel.forcing.reconstruct.setopts().core_channels
      kwargs.seed (1, 1) double {mustBeInteger, mustBeNonnegative}
      kwargs.split_manifest (1, 1) string = ""
      kwargs.n_gaps (1, 1) double {mustBeInteger, mustBePositive} = ...
         icemodel.forcing.reconstruct.setopts().plan_n_gaps
      kwargs.knot_candidates (1, :) double = ...
         icemodel.forcing.reconstruct.setopts().knot_candidates
      kwargs.max_donors (1, 1) double {mustBeInteger, mustBePositive} = ...
         icemodel.forcing.reconstruct.setopts().max_donors
      kwargs.max_donor_distance_km (1, 1) double {mustBePositive} = ...
         icemodel.forcing.reconstruct.setopts().max_donor_distance_km
       kwargs.max_donor_elev_diff_m (1, 1) double {mustBePositive} = ...
          icemodel.forcing.reconstruct.setopts().max_donor_elev_diff_m
       kwargs.cap_hours (1, 1) double {mustBePositive, ...
          icemodel.forcing.reconstruct.mustBeCapHours(kwargs.cap_hours)} = ...
          icemodel.forcing.reconstruct.setopts().cap_hours
       kwargs.cap_hours_by_channel (1, 1) struct = ...
          icemodel.forcing.reconstruct.setopts().cap_hours_by_channel
       kwargs.selection_fraction (1, 1) double = ...
          icemodel.forcing.reconstruct.setopts().selection_fraction
       kwargs.min_overlap_hours (1, 1) double {mustBePositive} = ...
          icemodel.forcing.reconstruct.setopts().min_overlap_hours
       kwargs.max_lag_hours (1, 1) double {mustBeNonnegative} = ...
          icemodel.forcing.reconstruct.setopts().max_lag_hours
       kwargs.min_lag_gain (1, 1) double {mustBeNonnegative} = ...
          icemodel.forcing.reconstruct.setopts().min_lag_gain
       kwargs.max_extrapolation_fraction (1, 1) double ...
          {mustBeNonnegative} = ...
          icemodel.forcing.reconstruct.setopts().max_extrapolation_fraction
       kwargs.rmse_improvement (1, 1) double {mustBeNonnegative} = ...
          icemodel.forcing.reconstruct.setopts().rmse_improvement
       kwargs.min_variability_ratio (1, 1) double {mustBePositive} = ...
          icemodel.forcing.reconstruct.setopts().min_variability_ratio
       kwargs.max_variability_ratio (1, 1) double {mustBePositive} = ...
          icemodel.forcing.reconstruct.setopts().max_variability_ratio
       kwargs.min_coverage (1, 1) double = ...
          icemodel.forcing.reconstruct.setopts().min_coverage
       kwargs.lapse_rate (1, 1) double = ...
          icemodel.forcing.reconstruct.setopts().lapse_rate
       kwargs.elevation_threshold_m (1, 1) double {mustBePositive} = ...
          icemodel.forcing.reconstruct.setopts().elevation_threshold_m
       kwargs.tair_for_pressure (1, 1) double {mustBePositive} = ...
          icemodel.forcing.reconstruct.setopts().tair_for_pressure
       kwargs.min_season_samples (1, 1) double {mustBePositive} = ...
          icemodel.forcing.reconstruct.setopts().min_season_samples
       kwargs.climatology_window_days (1, 1) double {mustBePositive} = ...
          icemodel.forcing.reconstruct.setopts().climatology_window_days
       kwargs.climatology_min_support (1, 1) double {mustBePositive} = ...
          icemodel.forcing.reconstruct.setopts().climatology_min_support
       kwargs.synthetic_context_hours (1, 1) double {mustBeNonnegative} = ...
          icemodel.forcing.reconstruct.setopts().synthetic_context_hours
       kwargs.jump_factor (1, 1) double {mustBePositive} = ...
          icemodel.forcing.reconstruct.setopts().jump_factor
       kwargs.toa_dark_wm2 (1, 1) double {mustBePositive} = ...
          icemodel.forcing.reconstruct.setopts().toa_dark_wm2
    end

   % Validate the split control before the no-core path can return a blocked plan.
   mustBeInRange(kwargs.selection_fraction, 0, 1, 'exclusive')

   series = target.series;
   times = series.Properties.RowTimes;

   % One split and one census govern every channel of this station. Only
   % jointly finite core-channel years may enter the split; a peripheral
   % channel outside the station's core record cannot fabricate evaluation
   % support for an otherwise disjoint record.
   precipitation = intersect(kwargs.channels, ...
      icemodel.forcing.helpers.precipitationVariables(), 'stable');
   if ~isempty(precipitation)
      error('icemodel:reconstruct:stationMethodPlan:plannedPrecipitation', ...
         ['precipitation is adopted as one total source and partitioned ' ...
         'afterward; statistical planning is forbidden for: %s'], ...
         strjoin(precipitation, ", "));
   end
   channels = intersect(kwargs.channels, ...
      string(series.Properties.VariableNames), 'stable');
   try
       census = icemodel.forcing.reconstruct.gapCensus(series, ...
          channels=channels, core_channels=kwargs.core_channels, ...
          cap_hours=kwargs.cap_hours, ...
          cap_hours_by_channel=kwargs.cap_hours_by_channel, ...
          latitude=target.location.lat_wgs84, ...
          longitude=target.location.lon_wgs84, ...
          toa_dark_wm2=kwargs.toa_dark_wm2);
   catch err
      if strcmp(err.identifier, ...
            'icemodel:reconstruct:gapCensus:noFiniteCore')
         plan = blockedPlan(target.station, channels, kwargs.seed, ...
            kwargs.selection_fraction, ...
            "no jointly finite core-channel sample");
         return
      end
      rethrow(err)
   end
   in_record = times >= census.record_start & times <= census.record_end;
   core = unique(kwargs.core_channels, 'stable');
   if isempty(core)
      split_support = in_record;
   else
      split_support = in_record;
      for channel = reshape(core, 1, [])
         split_support = split_support & isfinite(series.(channel));
      end
   end
   observed_years = unique(year(times(split_support)));
   split = icemodel.forcing.reconstruct.validationSplit( ...
      observed_years.', station=target.station, seed=kwargs.seed, ...
      selection_fraction=kwargs.selection_fraction, ...
      manifest_file=kwargs.split_manifest);

   % Geometry gate once per donor (POLICY B4), nearest first; GC-Net
   % donors are pre-masked to origin-observed samples. The planner enforces
    % target exclusion itself so direct callers cannot leak held-out truth.
    if ~isempty(donors)
       donor_station = icemodel.forcing.helpers.normalizedFileToken( ...
          icemodel.forcing.helpers.gcnetVandecruxStation( ...
          string({donors.station})));
       target_station = icemodel.forcing.helpers.normalizedFileToken( ...
          icemodel.forcing.helpers.gcnetVandecruxStation(target.station));
       donors = donors(donor_station ~= target_station);
    end
   eligible = donorGeometry(target.location, donors, ...
      kwargs.max_donor_distance_km, kwargs.max_donor_elev_diff_m);

   channel_plans = cell(numel(channels), 1);
   for c = 1:numel(channels)
      channel = channels(c);
      channel_plans{c} = planChannel(channel, series, times, target, ...
         eligible, proxies, split, census, kwargs);
   end

   plan = struct('station', target.station, 'split', split, ...
      'census_summary', census.summary, ...
      'fixed_methods', fixedMethodRecords(), ...
      'channels', vertcat(channel_plans{:}));
end

function plan = blockedPlan(station, channels, seed, selection_fraction, reason)
   %BLOCKEDPLAN Preserve a triageable plan when no core record exists.
   denials = denialTable({{'all_methods', NaN, 'all', char(reason)}}, 1);
   channel_plans = cell(numel(channels), 1);
   for c = 1:numel(channels)
       channel_plans{c} = struct('channel', channels(c), ...
          'methods', emptyMethods(), 'n_candidates', 0, ...
          'denials', denials, ...
          'proxy_calibrations', emptyProxyCalibrations());
   end
   split = struct('station', station, 'seed', seed, ...
      'selection_fraction', selection_fraction, ...
      'years_selection', zeros(1, 0), ...
      'years_evaluation', zeros(1, 0));
   plan = struct('station', station, 'split', split, ...
      'census_summary', table(), 'fixed_methods', fixedMethodRecords(), ...
      'channels', vertcat(channel_plans{:}));
end

function methods = fixedMethodRecords()
   %FIXEDMETHODRECORDS Persist contexts for non-fitted policy tiers.
   methods = repmat(struct('name', "", 'audit_context_id', "", ...
      'policy_basis', ""), 2, 1);
   methods(1) = struct('name', "darkness_zero", ...
      'audit_context_id', "darkness_zero", 'policy_basis', "POLICY B2");
   methods(2) = struct('name', "bounded_interp", ...
      'audit_context_id', "bounded_interp", 'policy_basis', "POLICY B3");
end

function eligible = donorGeometry(location, donors, max_km, max_delev)
   %DONORGEOMETRY Apply the policy geometry gate and sort by distance.
   n = numel(donors);
   km = nan(n, 1);
   delev = nan(n, 1);
   for k = 1:n
      loc = donors(k).location;
      km(k) = haversineKm(location.lat_wgs84, location.lon_wgs84, ...
         loc.lat_wgs84, loc.lon_wgs84);
      delev(k) = location.elev_m - loc.elev_m;
   end
   pass = km <= max_km & abs(delev) <= max_delev;
   [~, order] = sort(km);
   order = order(pass(order));
   eligible = donors(order);
   for k = 1:numel(eligible)
      eligible(k).distance_km = km(order(k));
      eligible(k).delev_m = delev(order(k));
   end
end

function d = haversineKm(lat1, lon1, lat2, lon2)
   %HAVERSINEKM Great-circle distance in kilometers.
   d = 2 * 6371 * asin(sqrt(sind((lat2 - lat1) / 2).^2 + ...
      cosd(lat1) .* cosd(lat2) .* sind((lon2 - lon1) / 2).^2));
end

function channel_plan = planChannel(channel, series, times, target, ...
      eligible, proxies, split, census, kwargs)
   %PLANCHANNEL Fit, grade, admit, and rank candidate methods for a channel.
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   x = series.(channel);

   % Upward shortwave is an algebraic dependency, not a statistically
   % admitted donor/proxy candidate. The production driver derives its
   % missing samples only after albedo and swd have completed every tier.
   if channel == "swu"
      channel_plan = struct('channel', channel, ...
         'methods', emptyMethods(), 'n_candidates', 0, ...
         'denials', denialTable(cell(0, 1), 0), ...
         'proxy_calibrations', emptyProxyCalibrations());
      return
   end

    % Draw each selection stratum before fitting, then withhold the union of
    % those samples from every fitted candidate. This makes the synthetic
    % observations genuinely held out instead of leaking into calibration.
    buckets = 1:(numel(icemodel.forcing.reconstruct.bucketEdges()) - 1);
    strata = validationStrata(buckets);
    selection_draws = drawSynthetic(series, channel, census, ...
       split.years_selection, strata, kwargs, target.location);
    withheld = false(height(series), 1);
    for q = 1:numel(selection_draws)
       if ~isempty(selection_draws{q})
          withheld = withheld | selection_draws{q}.mask;
       end
    end
    x_fit = x;
    x_fit(withheld) = NaN;

    % Candidate estimates span the full axis, but every fitted parameter and
    % the climatology pool see selection years minus the withheld draws.
    [candidates, candidate_denials, proxy_calibrations] = ...
       buildCandidates(channel, x_fit, times, target, eligible, proxies, ...
       split, kwargs, codes);

    % Grade every candidate on the same precomputed selection draws; the
    % climatology candidate also supplies the longer-gap policy baseline.
    % The calendar year of every axis sample is grading-pass-invariant, so
    % one full-axis year() evaluation here replaces the per-candidate
    % recomputation gradeCandidate used to perform (an O(axis) datetime
    % pass per candidate that showed up in the production driver profile).
    axis_years = year(times);
    graded = cell(numel(candidates), 1);
    for m = 1:numel(candidates)
       graded{m} = gradeCandidate(candidates(m), channel, series, ...
          strata, selection_draws, kwargs, "selection", target.location, ...
          split.years_selection, axis_years);
    end

    % Grade the policy's short-gap persistence baseline on the same blocked
    % draws. Longer strata use the climatology candidate as their baseline.
    persistence = emptyCandidate();
    persistence.name = "persistence";
    persistence.code = codes.bounded_interp;
    persistence_graded = gradeCandidate(persistence, channel, series, ...
       strata, selection_draws, kwargs, "selection", target.location, ...
       split.years_selection, axis_years);

    % Admit per stratum against its policy baseline, then rank admitted
    % methods by held-out RMSE within each bucket.
    baseline_idx = find(string({candidates.name}) == "climatology", 1);
    [methods, gate_denials] = admitAndRank(candidates, graded, ...
       baseline_idx, persistence_graded, strata, channel, kwargs);
    [methods, donor_limit_denials] = limitDonorMethods( ...
       methods, kwargs.max_donors);
    denials = [candidate_denials; gate_denials; donor_limit_denials];

    % Selection evidence cannot authorize production use when the record
    % has no disjoint whole evaluation year.
    if isempty(split.years_evaluation)
       rows = cell(numel(candidates), 1);
       for m = 1:numel(candidates)
          rows{m} = {char(candidates(m).name), NaN, 'all', ...
             'no disjoint evaluation year'};
       end
       denials = [denials; denialTable(rows, numel(candidates))];
        methods = methods([]);
        channel_plan = struct('channel', channel, 'methods', methods, ...
           'n_candidates', numel(candidates), 'denials', denials, ...
           'proxy_calibrations', proxy_calibrations);
        return
    end

    % Final numbers for the report come from evaluation-year draws of the
    % ranked winners only (never the selection draws).
    evaluation_draws = drawSynthetic(series, channel, census, ...
       split.years_evaluation, strata, kwargs, target.location);
    for m = 1:numel(methods)
       use_stratum = strata.bucket == methods(m).buckets ...
          & strata.season == methods(m).seasons;
        methods(m).evaluation = gradeCandidate( ...
           candidates(methods(m).candidate_index), channel, series, ...
           strata(use_stratum, :), evaluation_draws(use_stratum), kwargs, ...
           "evaluation", target.location, split.years_evaluation, ...
           axis_years);
   end

    channel_plan = struct('channel', channel, 'methods', methods, ...
       'n_candidates', numel(candidates), 'denials', denials, ...
       'proxy_calibrations', proxy_calibrations);
end

function [candidates, denials, proxy_calibrations] = buildCandidates( ...
      channel, x_fit, times, target, eligible, proxies, split, kwargs, codes)
   %BUILDCANDIDATES Assemble fitted candidate estimates for one channel.
    % Every geometry-eligible channel-bearing donor must reach held-out
    % validation; max_donors applies only after that scientific ranking.
    has_channel = arrayfun(@(d) ismember(channel, ...
       string(d.series.Properties.VariableNames)), eligible);
    channel_donors = eligible(has_channel);
    n_pairs = numel(channel_donors) * (numel(channel_donors) - 1) / 2;

   % Worst case: donor/knot fits, every donor pair, proxies, estimator,
   % and climatology.
   max_candidates = numel(channel_donors) * ...
      numel(kwargs.knot_candidates) + n_pairs + numel(proxies) + 2;
   candidates = repmat(emptyCandidate(), max_candidates, 1);
   n = 0;
   denial_rows = cell(max_candidates, 1);
   n_denials = 0;
   proxy_calibrations = repmat( ...
      struct('source', "", 'parameters', struct()), numel(proxies), 1);
   n_proxy_calibrations = 0;

    % Donor transfers: geometry order is deterministic, but validation skill
    % decides which donors survive, with one candidate per knot setting.
   for k = 1:numel(channel_donors)
      donor = channel_donors(k);
      d = donorChannel(donor, channel, times);
      [d, elevation_adjustment] = adjustDonorElevation(channel, x_fit, ...
         d, times, split.years_selection, donor.delev_m, kwargs);
      for knots = kwargs.knot_candidates
         try
             transfer = icemodel.forcing.reconstruct.fitDonorTransfer( ...
                times, x_fit, d, channel, fit_years=split.years_selection, ...
                 knots=knots, max_lag_hours=kwargs.max_lag_hours, ...
                 min_lag_gain=kwargs.min_lag_gain, ...
                 min_overlap_hours=kwargs.min_overlap_hours, ...
                 target_location=target.location, ...
                 donor_location=donor.location, ...
                 toa_dark_wm2=kwargs.toa_dark_wm2);
         catch err
            if strcmp(err.identifier, ...
                  'icemodel:reconstruct:fitDonorTransfer:insufficientOverlap')
               n_denials = n_denials + 1;
               denial_rows{n_denials} = { ...
                  char("donor:" + donor.station), NaN, 'all', ...
                  char(string(err.message))};
               break
            end
            rethrow(err)
         end
         transfer.donor_station = donor.station;
         transfer.elevation_adjustment = elevation_adjustment;
         % K-transect visit heights and AWS generations do not define an
         % approved transfer correction, but they must survive in the fitted
         % candidate so mixed sensor geometries remain reviewable.
         metadata = donor.series.Properties.UserData;
         if isstruct(metadata)
            transfer.height_provenance = struct( ...
               'sensor_heights', metadataField(metadata, 'sensor_heights'), ...
               'aws_types', metadataField(metadata, 'aws_types'));
         end
         n = n + 1;
         code_name = "station_transfer";
         if knots > 0
            code_name = "spline_adjustment";
         end
         candidates(n) = struct( ...
            'name', sprintf("donor:%s:k%d", donor.station, knots), ...
            'code', codes.(char(code_name)), ...
            'estimate', icemodel.forcing.reconstruct.applyDonorTransfer( ...
             transfer, times, d, max_extrapolation_fraction= ...
             kwargs.max_extrapolation_fraction), ...
            'parameters', transfer);
      end
   end

   % A target bracketed in elevation by two donors gets an explicit
   % linearly interpolated candidate before overlap calibration. The
   % target TOA reconversion factor is pair-invariant, so it is evaluated
   % once here instead of inside the O(donors^2) pair loop below.
   bracket_target_toa = zeros(0, 1);
   if channel == "swd" && numel(channel_donors) > 1
      bracket_target_toa = icemodel.forcing.reconstruct.toaIrradiance( ...
         times, target.location.lat_wgs84, target.location.lon_wgs84);
   end
   for a = 1:(numel(channel_donors) - 1)
      for b = (a + 1):numel(channel_donors)
         lower = channel_donors(a);
         upper = channel_donors(b);
         if lower.delev_m * upper.delev_m >= 0
            continue
         end
         z1 = target.location.elev_m - lower.delev_m;
         z2 = target.location.elev_m - upper.delev_m;
         weight = (target.location.elev_m - z1) / (z2 - z1);
         lower_values = donorChannel(lower, channel, times);
         upper_values = donorChannel(upper, channel, times);
         donor_location = target.location;
         if channel == "swd"
            lower_values = icemodel.forcing.reconstruct.clearSkyIndex( ...
               times, lower_values, lower.location, ...
               toa_dark_wm2=kwargs.toa_dark_wm2);
            upper_values = icemodel.forcing.reconstruct.clearSkyIndex( ...
               times, upper_values, upper.location, ...
               toa_dark_wm2=kwargs.toa_dark_wm2);
            bracket = ((1 - weight) * lower_values ...
               + weight * upper_values) .* bracket_target_toa;
         else
            bracket = (1 - weight) * lower_values + weight * upper_values;
         end
         pair_name = lower.station + "+" + upper.station;
         try
            transfer = icemodel.forcing.reconstruct.fitDonorTransfer( ...
               times, x_fit, bracket, channel, ...
               fit_years=split.years_selection, knots=0, ...
                max_lag_hours=kwargs.max_lag_hours, ...
                min_lag_gain=kwargs.min_lag_gain, ...
                min_overlap_hours=kwargs.min_overlap_hours, ...
                target_location=target.location, ...
                donor_location=donor_location, ...
                toa_dark_wm2=kwargs.toa_dark_wm2);
         catch err
            if strcmp(err.identifier, ...
                  'icemodel:reconstruct:fitDonorTransfer:insufficientOverlap')
               n_denials = n_denials + 1;
               denial_rows{n_denials} = { ...
                  char("bracket:" + pair_name), NaN, 'all', ...
                  char(string(err.message))};
               continue
            end
            rethrow(err)
         end
         transfer.donor_station = pair_name;
         transfer.donor_stations = [lower.station, upper.station];
         transfer.elevation_adjustment = struct( ...
            'mode', "elevation_bracket", 'station_elevations_m', ...
            [z1, z2], 'target_elevation_m', target.location.elev_m, ...
            'upper_weight', weight);
         n = n + 1;
         candidates(n) = struct('name', "bracket:" + pair_name, ...
            'code', codes.station_transfer, ...
            'estimate', ...
            icemodel.forcing.reconstruct.applyDonorTransfer( ...
            transfer, times, bracket, max_extrapolation_fraction= ...
            kwargs.max_extrapolation_fraction), ...
            'parameters', transfer);
      end
   end

   % Calibrated proxies (MAR default, MERRA-2 alternative) use target solar
   % geometry for shortwave overlap; the lwd empirical estimator follows the
   % same calibration path with no solar screen.
   proxy_target_toa = zeros(0, 1);
   proxy_target_elevation = zeros(0, 1);
   if ismember(channel, ["swd", "swu"])
      proxy_target_toa = icemodel.forcing.reconstruct.toaIrradiance( ...
         times, target.location.lat_wgs84, target.location.lon_wgs84);
   end
   if channel == "swd"
      % D-28: the elevation-binned swd calibration needs the SIGNED solar
      % elevation (TOA clamps twilight to zero, hiding the very band the
      % bins correct), from the same solarElevation SSOT that
      % toaIrradiance scales.
      proxy_target_elevation = icemodel.forcing.helpers.solarElevation( ...
         times, target.location.lat_wgs84, target.location.lon_wgs84);
   end
   for p = 1:numel(proxies)
      proxy = proxies(p);
      if ~ismember(channel, string(proxy.series.Properties.VariableNames))
         continue
      end
      model = alignToAxis(proxy.series, channel, times);
      calibration = icemodel.forcing.reconstruct.fitProxyCalibration( ...
         times, x_fit, model, channel, fit_years=split.years_selection, ...
         min_season_samples=kwargs.min_season_samples, ...
         min_light_wm2=kwargs.toa_dark_wm2, ...
         target_toa=proxy_target_toa, ...
         target_elevation=proxy_target_elevation);
      calibration.source = proxy.name;
      n_proxy_calibrations = n_proxy_calibrations + 1;
      proxy_calibrations(n_proxy_calibrations) = struct( ...
         'source', proxy.name, 'parameters', calibration);
      if ~isscalar(calibration.n_overlap) ...
            || ~isfinite(calibration.n_overlap) || calibration.n_overlap <= 0
         n_denials = n_denials + 1;
         denial_rows{n_denials} = {char("proxy:" + proxy.name), ...
            NaN, 'all', 'no positive finite target overlap'};
         continue
      end
      n = n + 1;
      candidates(n) = struct( ...
         'name', "proxy:" + proxy.name, ...
         'code', codes.(char(proxy.code_name)), ...
         'estimate', icemodel.forcing.reconstruct.applyProxyCalibration( ...
         calibration, times, model, ...
         target_elevation=proxy_target_elevation), ...
         'parameters', calibration);
   end
   if channel == "lwd" && all(ismember(["tair", "rh"], ...
         string(target.series.Properties.VariableNames)))
      raw = icemodel.forcing.reconstruct.lwdEstimator( ...
         target.series.tair, target.series.rh);
      calibration = icemodel.forcing.reconstruct.fitProxyCalibration( ...
         times, x_fit, raw, channel, fit_years=split.years_selection, ...
         min_season_samples=kwargs.min_season_samples, ...
         min_light_wm2=kwargs.toa_dark_wm2, ...
         target_toa=proxy_target_toa);
      calibration.source = "brutsaert";
      if ~isscalar(calibration.n_overlap) ...
            || ~isfinite(calibration.n_overlap) || calibration.n_overlap <= 0
         n_denials = n_denials + 1;
         denial_rows{n_denials} = {'estimator:brutsaert', ...
            NaN, 'all', 'no positive finite target overlap'};
      else
         n = n + 1;
         candidates(n) = struct('name', "estimator:brutsaert", ...
            'code', codes.empirical_estimator, ...
            'estimate', icemodel.forcing.reconstruct.applyProxyCalibration( ...
            calibration, times, raw), ...
            'parameters', calibration);
      end
   end

    % Climatology is always a candidate AND the admission baseline. Its pool
    % is restricted to non-withheld selection-year observations.
    climatology_pool = x_fit;
    climatology_pool(~ismember(year(times), split.years_selection)) = NaN;
    n = n + 1;
    candidates(n) = struct('name', "climatology", ...
       'code', codes.climatology, ...
      'estimate', icemodel.forcing.reconstruct.climatologyFill( ...
       times, climatology_pool, times, ...
       window_days=kwargs.climatology_window_days, ...
       min_support=kwargs.climatology_min_support, ...
      diurnal=(channel ~= "albedo")), ...
      'parameters', struct('kind', "doy_median", ...
         'fit_years', split.years_selection));
   candidates = candidates(1:n);
   proxy_calibrations = proxy_calibrations(1:n_proxy_calibrations);
   denials = denialTable(denial_rows, n_denials);
end

function value = metadataField(metadata, name)
   %METADATAFIELD Return optional donor metadata without inventing a value.
   value = [];
   if isfield(metadata, name)
      value = metadata.(name);
   end
end

function candidate = emptyCandidate()
   %EMPTYCANDIDATE Prototype candidate record.
   candidate = struct('name', "", 'code', uint8(255), ...
      'estimate', [], 'parameters', struct());
end

function [adjusted, record] = adjustDonorElevation(channel, target, donor, ...
      times, fit_years, dz, kwargs)
   %ADJUSTDONORELEVATION Fit the approved correction or use its fallback.
   adjusted = donor;
   record = struct('mode', "none", 'dz_m', dz, ...
      'fit_years', fit_years, 'n_overlap', 0, 'n_overlap_hours', 0);
   if abs(dz) <= kwargs.elevation_threshold_m ...
         || ~ismember(channel, ["tair", "psfc"])
      return
   end

   overlap = ismember(year(times), fit_years) ...
      & isfinite(target) & isfinite(donor);
   % Support-held donor repetitions integrate to elapsed support only when
   % multiplied by the target interval represented by each aligned row.
   dt_target_hours = hours(median(diff(times)));
   record.n_overlap = nnz(overlap);
   record.n_overlap_hours = record.n_overlap * dt_target_hours;
   enough = record.n_overlap_hours >= kwargs.min_overlap_hours;
   if channel == "tair"
      lapse_rate = kwargs.lapse_rate;
      record.mode = "fallback";
      if enough
         lapse_rate = median((target(overlap) - donor(overlap)) / dz, ...
            'omitnan');
         record.mode = "overlap_fit";
      end
      adjusted = icemodel.forcing.reconstruct.elevationAdjust( ...
         channel, donor, dz, lapse_rate=lapse_rate, ...
         threshold_m=kwargs.elevation_threshold_m);
      record.lapse_rate = lapse_rate;
      return
   end

   % Pressure overlap directly identifies the barometric multiplicative
   % factor (equivalently a fitted scale height) without an unstable
   % logarithmic temperature inversion.
   record.mode = "fallback";
   record.tair_for_pressure = kwargs.tair_for_pressure;
   if enough
      factor = median(target(overlap) ./ donor(overlap), 'omitnan');
      adjusted = donor * factor;
      record.mode = "overlap_fit";
      record.factor = factor;
   else
      adjusted = icemodel.forcing.reconstruct.elevationAdjust( ...
         channel, donor, dz, ...
         threshold_m=kwargs.elevation_threshold_m, ...
         tair_for_pressure=kwargs.tair_for_pressure);
   end
end

function strata = validationStrata(buckets)
   %VALIDATIONSTRATA Enumerate every required season x duration stratum.
   seasons = ["DJF"; "MAM"; "JJA"; "SON"];
   bucket_grid = repelem(buckets(:), numel(seasons));
   season_grid = repmat(seasons, numel(buckets), 1);
   strata = table(bucket_grid, season_grid, ...
      'VariableNames', {'bucket', 'season'});
end

function draws = drawSynthetic(series, channel, census, years, strata, ...
      kwargs, location)
    %DRAWSYNTHETIC Create one deterministic held-out draw per exact stratum.
    draws = cell(height(strata), 1);
    if isempty(years)
       % A one-year station cannot furnish disjoint selection/evaluation
       % years, so the caller will refuse every production method.
       return
    end
    seasons = ["DJF", "MAM", "JJA", "SON"];
    for q = 1:height(strata)
       try
          season_index = find(seasons == strata.season(q), 1);
          draws{q} = icemodel.forcing.reconstruct.syntheticMissingness( ...
             series, channel, census.runs, years=years, ...
             seed=kwargs.seed + 100 * strata.bucket(q) + ...
             1000 * season_index, n_gaps=kwargs.n_gaps, ...
             bucket=strata.bucket(q), season=strata.season(q), ...
             context_hours=kwargs.synthetic_context_hours, ...
             latitude=location.lat_wgs84, ...
             longitude=location.lon_wgs84, ...
             toa_dark_wm2=kwargs.toa_dark_wm2);
       catch err
          if strcmp(err.identifier, ...
                'icemodel:reconstruct:syntheticMissingness:emptyDurationPool')
             % This station has no empirical duration pool for the stratum.
             continue
          end
          rethrow(err)
       end
    end
end

function grade = gradeCandidate(candidate, channel, series, ...
      strata, draws, kwargs, label, location, years, axis_years)
   %GRADECANDIDATE Held-out metrics of one candidate per exact stratum.
   % axis_years carries year(series.Properties.RowTimes) precomputed once
   % by the caller: the axis never changes between grading calls, so the
   % hoist is exact while removing a per-candidate O(axis) datetime pass.
   grade = repmat(struct('bucket', NaN, 'metrics', [], 'inserted', 0, ...
       'max_duration_hours', NaN, 'by_stratum', table(), ...
       'truth', [], 'filled', [], 'season', "", 'label', label), ...
       height(strata), 1);
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   validation_series = series;
   % A validation pass may derive seam scales only from its own split years.
   outside_split = ~ismember(axis_years, years);
   validation_series.(channel)(outside_split) = NaN;
   for q = 1:numel(draws)
      % No held-out truth may set the native-observation seam scale used to
      % decide whether any candidate in this validation pass is admitted.
      if ~isempty(draws{q})
         validation_series.(channel)(draws{q}.mask) = NaN;
      end
   end
    for q = 1:height(strata)
       grade(q).bucket = strata.bucket(q);
       grade(q).season = strata.season(q);
       draw = draws{q};
       if isempty(draw) || draw.inserted == 0
          continue
       end
       truth = series.(channel)(draw.mask);
       if candidate.name == "persistence"
          % The shared helper masks every draw before finding the active gap's
          % anchor, so no caller can expose held-out truth to this baseline.
          filled = icemodel.forcing.reconstruct.persistenceEstimate( ...
             series.(channel), series.Properties.RowTimes, draw, draws);
       else
          filled = candidate.estimate(draw.mask);
       end
       provenance = repmat(double(candidate.code), numel(filled), 1);
       provenance(~isfinite(filled)) = double(codes.missing);
       metrics = icemodel.forcing.reconstruct.validationMetrics( ...
          truth, filled, draw.gaps, validation_series, channel, ...
          provenance=provenance, ...
          jump_factor=kwargs.jump_factor, ...
          latitude=location.lat_wgs84, longitude=location.lon_wgs84);
       grade(q).metrics = metrics.overall;
       grade(q).inserted = draw.inserted;
       grade(q).max_duration_hours = metrics.max_complete_gap_hours;
       grade(q).truth = truth;
       grade(q).filled = filled;
       row = [strata(q, :), metrics.overall];
       row.max_duration_hours = grade(q).max_duration_hours;
       grade(q).by_stratum = row;
   end
end

function [methods, denials] = admitAndRank(candidates, graded, ...
      baseline_idx, persistence_graded, strata, channel, kwargs)
   %ADMITANDRANK Admit and rank candidates within each exact stratum.
   % One record per season x duration bucket prevents performance from
   % another regime from changing the method order for the active gap.
   method_slots = cell(numel(candidates) * height(strata), 1);
   denial_slots = cell(numel(candidates) * height(strata), 1);
    n_denials = 0;
    n_methods = 0;
    mar_idx = find(lower(string({candidates.name})) == "proxy:mar", 1);
    for m = 1:numel(candidates)
      for q = 1:height(strata)
         metrics = graded{m}(q).metrics;
         if isempty(metrics)
            n_denials = n_denials + 1;
            denial_slots{n_denials} = {char(candidates(m).name), ...
               strata.bucket(q), char(strata.season(q)), ...
               'no held-out draw'};
            continue
         end
         if graded{m}(q).max_duration_hours <= 0
            n_denials = n_denials + 1;
            denial_slots{n_denials} = {char(candidates(m).name), ...
               strata.bucket(q), char(strata.season(q)), ...
               'no completely reconstructed held-out gap'};
            continue
         end
         baseline_name = "climatology";
         baseline_grade = graded{baseline_idx}(q);
         if strata.bucket(q) == 1
            baseline_name = "persistence";
            baseline_grade = persistence_graded(q);
         end
         comparison = icemodel.forcing.reconstruct.commonSupportSkill( ...
            graded{m}(q).truth, graded{m}(q).filled, ...
            baseline_grade.filled);
         min_common = max(2, ceil(kwargs.min_coverage * metrics.n));
         if comparison.n < min_common
            n_denials = n_denials + 1;
            denial_slots{n_denials} = {char(candidates(m).name), ...
               strata.bucket(q), char(strata.season(q)), ...
               sprintf('common baseline support %d below required %d', ...
               comparison.n, min_common)};
            continue
         end
         candidate_rmse = comparison.candidate_rmse;
         baseline_rmse = comparison.baseline_rmse;
          is_wind_donor = channel == "wspd" ...
             && (startsWith(candidates(m).name, "donor:") ...
             || startsWith(candidates(m).name, "bracket:"));
          if is_wind_donor
             if isempty(mar_idx) || isempty(graded{mar_idx}(q).metrics) ...
                   || graded{mar_idx}(q).max_duration_hours <= 0 ...
                   || ~isfinite(graded{mar_idx}(q).metrics.rmse)
                n_denials = n_denials + 1;
                denial_slots{n_denials} = {char(candidates(m).name), ...
                   strata.bucket(q), char(strata.season(q)), ...
                   'calibrated MAR wind benchmark unavailable'};
                continue
             end
              mar_comparison = ...
                 icemodel.forcing.reconstruct.commonSupportSkill( ...
                 graded{m}(q).truth, graded{m}(q).filled, ...
                 graded{mar_idx}(q).filled);
              if mar_comparison.n < min_common
                 n_denials = n_denials + 1;
                 denial_slots{n_denials} = {char(candidates(m).name), ...
                    strata.bucket(q), char(strata.season(q)), ...
                    sprintf(['common calibrated-MAR support %d below ' ...
                    'required %d'], mar_comparison.n, min_common)};
                 continue
              end
              mar_rmse = mar_comparison.baseline_rmse;
              if ~isfinite(mar_comparison.candidate_rmse) ...
                    || mar_comparison.candidate_rmse >= mar_rmse
                 n_denials = n_denials + 1;
                 denial_slots{n_denials} = {char(candidates(m).name), ...
                    strata.bucket(q), char(strata.season(q)), ...
                    sprintf(['held-out RMSE %.3g does not beat calibrated ' ...
                    'MAR wind %.3g'], mar_comparison.candidate_rmse, ...
                    mar_rmse)};
                 continue
              end
              % The baseline with the larger common-support error ratio is
              % the stricter fractional-improvement comparison.
              if ~isfinite(comparison.rmse_ratio) ...
                    || mar_comparison.rmse_ratio > comparison.rmse_ratio
                 baseline_name = "proxy:mar";
                 comparison = mar_comparison;
                 candidate_rmse = comparison.candidate_rmse;
                 baseline_rmse = mar_rmse;
              end
           end
           gate_baseline_rmse = baseline_rmse;
         if m == baseline_idx && baseline_name == "climatology"
            % Long-gap climatology admits against itself on
            % bias/bounds/coverage only; short-gap climatology must still
             % improve on persistence.
             gate_baseline_rmse = NaN;
          end
          gate_metrics = metrics;
          gate_metrics.rmse = candidate_rmse;
          gate = icemodel.forcing.reconstruct.admissionGate(channel, ...
             gate_metrics, gate_baseline_rmse, ...
            rmse_improvement=kwargs.rmse_improvement, ...
            min_variability_ratio=kwargs.min_variability_ratio, ...
            max_variability_ratio=kwargs.max_variability_ratio, ...
            min_coverage=kwargs.min_coverage, ...
            typical_magnitude=metrics.typical_magnitude);
         if ~gate.admit
            n_denials = n_denials + 1;
            denial_slots{n_denials} = {char(candidates(m).name), ...
               strata.bucket(q), char(strata.season(q)), ...
               char(strjoin(gate.reasons, "; "))};
            continue
         end
          selection = graded{m}(q).by_stratum;
          selection.rmse = candidate_rmse;
          selection.baseline = baseline_name;
          selection.baseline_rmse = baseline_rmse;
          selection.comparison_n = comparison.n;
          selection.fractional_improvement = ...
             comparison.fractional_improvement;
          n_methods = n_methods + 1;
         method_slots{n_methods} = struct('name', candidates(m).name, ...
            'code', candidates(m).code, 'candidate_index', m, ...
            'audit_context_id', channel + ":candidate-" + string(m) + ...
            ":" + strata.season(q) + ":b" + string(strata.bucket(q)), ...
            'parameters', candidates(m).parameters, ...
            'estimate', candidates(m).estimate, ...
            'buckets', strata.bucket(q), 'seasons', strata.season(q), ...
             'max_validated_hours', selection.max_duration_hours, ...
             'selection_rmse', candidate_rmse, ...
            'selection', selection, 'evaluation', [], ...
            'uncertainty', "not_provided");
      end
   end
   denials = denialTable(denial_slots, n_denials);
   methods = vertcat(method_slots{1:n_methods});
   if isempty(methods)
      methods = emptyMethods();
      return
   end
   skill = arrayfun(@(method) ...
      method.selection.fractional_improvement, methods);
   skill(~isfinite(skill)) = 0;
   [~, order] = sortrows([-skill(:), [methods.selection_rmse].']);
   methods = methods(order);
   % D-29: for swd, calibrated RCM proxies rank ahead of climatology
   % whenever a stratum admits both — the day-of-year median structurally
   % inserts clearer-than-context days into cloudy weeks (cen fills
   % +0.09 median CSI vs observed neighbors), a context error the skill
   % ranking above cannot see. Deterministic post-ranking swap only;
   % climatology's POLICY B8 admission-baseline role is untouched.
   if channel == "swd"
      methods = preferProxyOverClimatology(methods);
   end
end

function methods = preferProxyOverClimatology(methods)
   %PREFERPROXYOVERCLIMATOLOGY Apply the D-29 swd candidate preference.
   % Within each validation stratum, the rank slots occupied by the
   % climatology/calibrated-proxy family are re-filled with the proxies
   % first (in their existing skill order) and climatology after them.
   % Every other method keeps its exact slot, so donor ordering and all
   % cross-stratum ranks are untouched; strata where only one of the two
   % is admitted are left alone.
   names = string({methods.name}).';
   is_proxy = startsWith(names, "proxy:");
   is_climatology = names == "climatology";
   family = is_proxy | is_climatology;
   if ~any(is_proxy) || ~any(is_climatology)
      return
   end
   buckets = [methods.buckets].';
   seasons = string({methods.seasons}).';
   strata = unique(table(buckets(family), seasons(family), ...
      'VariableNames', {'bucket', 'season'}), 'rows', 'stable');
   for q = 1:height(strata)
      slots = find(family & buckets == strata.bucket(q) ...
         & seasons == strata.season(q));
      if any(is_proxy(slots)) && any(is_climatology(slots))
         % Same slots, reordered occupants: proxies keep their relative
         % skill order; climatology entries follow them.
         methods(slots) = methods( ...
            [slots(is_proxy(slots)); slots(is_climatology(slots))]);
      end
   end
end

function [methods, denials] = limitDonorMethods(methods, max_donors)
   %LIMITDONORMETHODS Apply the donor cap within each validation stratum.
   rows = cell(numel(methods), 1);
   n_rows = 0;
   keep = true(numel(methods), 1);
   strata = unique(table([methods.buckets].', ...
      string({methods.seasons}).', ...
      'VariableNames', {'bucket', 'season'}), 'rows', 'stable');
   for q = 1:height(strata)
      % Globally RMSE-sorted methods retain that order inside this stratum,
      % while unrelated seasons and duration buckets get independent budgets.
      selected = strings(1, 0);
      members = find([methods.buckets].' == strata.bucket(q) ...
         & string({methods.seasons}).' == strata.season(q));
      for m = reshape(members, 1, [])
         donors = methodDonors(methods(m));
         if isempty(donors)
            continue
         end
         new_donors = setdiff(donors, selected, 'stable');
         if numel(selected) + numel(new_donors) > max_donors
            keep(m) = false;
            n_rows = n_rows + 1;
            rows{n_rows} = {char(methods(m).name), methods(m).buckets, ...
               char(methods(m).seasons), ...
               'excluded by max_donors after held-out ranking'};
            continue
         end
         selected = unique([selected, donors], 'stable');
      end
   end
   methods = methods(keep);
   denials = denialTable(rows, n_rows);
end

function donors = methodDonors(method)
   %METHODDONORS Return the station identities consumed by one method.
   donors = strings(1, 0);
   parameters = method.parameters;
   if ~isstruct(parameters)
      return
   end
   if isfield(parameters, 'donor_stations') ...
         && ~isempty(parameters.donor_stations)
      donors = string(parameters.donor_stations);
   elseif isfield(parameters, 'donor_station') ...
         && ~isempty(parameters.donor_station)
      donors = string(parameters.donor_station);
   end
   donors = unique(donors(:).', 'stable');
end

function methods = emptyMethods()
   %EMPTYMETHODS Return the schema-stable empty admitted-method array.
   methods = repmat(struct('name', "", 'code', uint8(255), ...
      'candidate_index', NaN, 'audit_context_id', "", ...
      'parameters', struct(), 'estimate', [], ...
      'buckets', [], 'seasons', "all", 'selection_rmse', NaN, ...
      'max_validated_hours', [], 'selection', table(), ...
      'evaluation', [], 'uncertainty', "not_provided"), 0, 1);
end

function calibrations = emptyProxyCalibrations()
   %EMPTYPROXYCALIBRATIONS Return the schema-stable calibration registry.
   calibrations = repmat(struct('source', "", 'parameters', struct()), 0, 1);
end

function denials = denialTable(rows, n_rows)
   %DENIALTABLE Build the schema-stable candidate-decision table.
   if n_rows == 0
      denials = table('Size', [0 4], 'VariableTypes', ...
         {'cellstr', 'double', 'cellstr', 'cellstr'}, ...
         'VariableNames', {'candidate', 'bucket', 'season', 'reasons'});
   else
      denials = cell2table(vertcat(rows{1:n_rows}), 'VariableNames', ...
         {'candidate', 'bucket', 'season', 'reasons'});
   end
end

function values = donorChannel(donor, channel, times)
   %DONORCHANNEL Align one channel and enforce its origin-observed mask.
   values = alignToAxis(donor.series, channel, times);
   if ~isfield(donor, 'observed_mask') || isempty(donor.observed_mask)
      return
   end

   % A timetable preserves independent origin flags for each variable.
   if istimetable(donor.observed_mask)
      if ~ismember(channel, ...
            string(donor.observed_mask.Properties.VariableNames))
         values(:) = NaN;
         return
      end
      observed = alignToAxis(donor.observed_mask, channel, times);
      values(observed ~= 1) = NaN;
      return
   end

   % Legacy vector masks are accepted on either donor or target cadence.
   observed = donor.observed_mask(:);
   if numel(observed) == height(donor.series)
      mask_series = timetable(donor.series.Properties.RowTimes, ...
         logical(observed), 'VariableNames', {char(channel)});
      observed = alignToAxis(mask_series, channel, times);
   elseif numel(observed) ~= numel(times)
      error('icemodel:reconstruct:stationMethodPlan:badObservedMask', ...
         'observed_mask must match the donor or target time axis');
   end
   values(observed ~= 1) = NaN;
end

function aligned = alignToAxis(donor_series, channel, times)
   %ALIGNTOAXIS Sample one donor/proxy channel onto the target time axis.
   % Support-held alignment: every target sample takes the donor posting
   % whose support interval covers it — the same interval-support
   % convention the target's own 15-minute axis uses (resampleMetTimestep
   % holds each hourly posting across its four quarter-hour samples). An
   % exact join would leave a coarser donor covering only one target
   % sample per posting, which starves the coverage gate for reasons that
   % have nothing to do with method skill. No sample is invented: postings
   % never stretch past their own cadence.
   donor_times = donor_series.Properties.RowTimes;
   aligned = nan(numel(times), 1);
   if isscalar(donor_times)
      [tf, loc] = ismember(times, donor_times);
      aligned(tf) = donor_series.(channel)(loc(tf));
      return
   end
   dt_donor = median(diff(donor_times));
   snapped = donor_times(1) + ...
      floor((times - donor_times(1)) / dt_donor) * dt_donor;
   [tf, loc] = ismember(snapped, donor_times);
   tf = tf & times >= donor_times(1) ...
      & times < donor_times(end) + dt_donor;
   aligned(tf) = donor_series.(channel)(loc(tf));
end
