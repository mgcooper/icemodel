function calibration = fitProxyCalibration(times, x_obs, x_model, channel, kwargs)
   %FITPROXYCALIBRATION Calibrate a model proxy on its observed overlap.
   %
   %  calibration = icemodel.forcing.reconstruct.fitProxyCalibration( ...
   %     times, x_obs, x_model, "lwd", fit_years=split.years_selection)
   %
   % Role
   %  The proxy-calibration tier's overlap bias fit (POLICY B5):
   %  additive median bias per season for state channels and a
   %  multiplicative ratio for shortwave and wind speed. Multiplication
   %  preserves shortwave shape and wind's nonnegative support. Fits use
   %  only the caller's fit years. Seasons with fewer than the minimum
   %  overlap samples inherit the annual correction; with no usable
   %  overlap at all the calibration is the recorded identity — B5
   %  requires calibration and evidence, not a mandatory nonzero
   %  adjustment.
   %
   % Name-value
   %  fit_years : calendar years eligible for fitting (required).
   %  min_season_samples : per-season overlap needed for a seasonal
   %     correction (default 300, a Section-C parameter).
   %  target_toa : target-station top-of-atmosphere irradiance; required
   %     for multiplicative shortwave calibration.
   %  target_elevation : target-station solar elevation (degrees, signed;
   %     TOA clamps twilight to zero, which would hide the very band the
   %     bins correct). When supplied for swd it additionally fits the
   %     D-28 elevation-binned ratios; omitted, the record keeps the
   %     legacy single-ratio shape.
   %
   % Returns
   %  calibration : struct — channel, mode ("additive"|"multiplicative"),
   %     per-season correction, n_overlap, identity flag. Binned swd
   %     records add: version (2), bin_edges_deg, per-season
   %     binned_corrections and binned_counts (one entry per elevation
   %     bin), and min_bin_samples — all read by applyProxyCalibration
   %     through a field-presence guard so legacy single-ratio records
   %     stay applicable (D-28 backward compatibility).
   %
   % See also: icemodel.forcing.reconstruct.applyProxyCalibration,
   %  icemodel.forcing.reconstruct.solarElevationBands

   arguments
      times datetime
      x_obs (:, 1) double
      x_model (:, 1) double
      channel (1, 1) string
      kwargs.fit_years (1, :) double {mustBeInteger, mustBeNonempty}
      kwargs.min_season_samples (1, 1) double {mustBePositive} = ...
         icemodel.forcing.reconstruct.setopts().min_season_samples
      kwargs.min_light_wm2 (1, 1) double {mustBePositive} = ...
         icemodel.forcing.reconstruct.setopts().toa_dark_wm2
      kwargs.target_toa (:, 1) double = zeros(0, 1)
      kwargs.target_elevation (:, 1) double = zeros(0, 1)
   end

   % A supplied elevation vector must cover the axis exactly; a silent
   % length mismatch would disable the binned fit without any signal.
   if ~isempty(kwargs.target_elevation) ...
         && numel(kwargs.target_elevation) ~= numel(times)
      error('icemodel:reconstruct:fitProxyCalibration:targetElevationSize', ...
         ['target_elevation must contain one target-station solar ' ...
         'elevation value per time when supplied']);
   end

   % Shortwave and wind speed correct multiplicatively: scaling preserves
   % the solar cycle and prevents a negative wind estimate. Only shortwave
   % ratios need the meaningful-light denominator screen below.
   multiplicative = ismember(channel, ["swd", "swu", "wspd"]);
   in_years = ismember(year(times), kwargs.fit_years);
   overlap = in_years & isfinite(x_obs) & isfinite(x_model);
   if ismember(channel, ["swd", "swu"])
      if numel(kwargs.target_toa) ~= numel(times)
         error('icemodel:reconstruct:fitProxyCalibration:targetToaRequired', ...
            ['target_toa must contain one target-station irradiance ' ...
            'value per time for shortwave calibration']);
      end
      overlap = overlap & kwargs.target_toa >= kwargs.min_light_wm2 ...
         & x_model > 0;
   end

   season = icemodel.forcing.reconstruct.seasonOf(times);
   annual = identityCorrection(multiplicative);
   if any(overlap)
      annual = oneCorrection(x_obs(overlap), x_model(overlap), ...
         multiplicative);
   end
   corrections = struct();
   for name = ["DJF", "MAM", "JJA", "SON"]
      in_season = overlap & season == name;
      corrections.(char(name)) = annual;
      if nnz(in_season) >= kwargs.min_season_samples
         corrections.(char(name)) = oneCorrection(x_obs(in_season), ...
            x_model(in_season), multiplicative);
      end
   end

   calibration = struct('channel', channel, ...
      'mode', modeName(multiplicative), ...
      'fit_years', kwargs.fit_years, ...
      'corrections', corrections, 'n_overlap', nnz(overlap), ...
      'identity', ~any(overlap));

   % D-28 twilight-shape fix, swd only: one seasonal ratio mixes regimes
   % with opposite biases (RCM proxies run 4-20x LOW near solar midnight
   % and 1.5-2x HIGH on the morning shoulder), so swd additionally fits a
   % ratio per solar-elevation band when the caller supplies the signed
   % elevation. The record keeps the seasonal scalar alongside the bins
   % so elevation-less consumers (and legacy readers) stay correct.
   if channel == "swd" && ~isempty(kwargs.target_elevation)
      bands = icemodel.forcing.reconstruct.solarElevationBands();
      edges = bands.calibration_bin_edges_deg;
      n_bins = numel(edges) - 1;
      bin_index = discretize(kwargs.target_elevation, edges);
      % Band membership screens on a positive proxy denominator only: the
      % meaningful-sun TOA screen guarding the single seasonal ratio
      % would empty the twilight bands whose bias this fit corrects
      % (stations measure 15-38 W m^-2 there while TOA reads zero).
      band_overlap = in_years & isfinite(x_obs) & isfinite(x_model) ...
         & x_model > 0;
      binned_corrections = struct();
      binned_counts = struct();
      for name = ["DJF", "MAM", "JJA", "SON"]
         % Bands without enough overlap inherit the season's single ratio
         % (which itself already fell back to the annual fit if thin), so
         % application degrades gracefully instead of extrapolating a
         % noisy band median.
         ratios = repmat(corrections.(char(name)), 1, n_bins);
         counts = zeros(1, n_bins);
         for b = 1:n_bins
            members = band_overlap & season == name & bin_index == b;
            counts(b) = nnz(members);
            if counts(b) >= bands.min_bin_samples
               ratios(b) = oneCorrection(x_obs(members), ...
                  x_model(members), multiplicative);
            end
         end
         binned_corrections.(char(name)) = ratios;
         binned_counts.(char(name)) = counts;
      end
      % Version stamp plus explicit edges make the persisted record
      % self-describing: applyProxyCalibration keys on field presence, so
      % an old single-ratio record (no version, no bins) still applies.
      calibration.version = 2;
      calibration.bin_edges_deg = edges;
      calibration.binned_corrections = binned_corrections;
      calibration.binned_counts = binned_counts;
      calibration.min_bin_samples = bands.min_bin_samples;
   end
end

function name = modeName(multiplicative)
   %MODENAME Correction mode label.
   name = "additive";
   if multiplicative
      name = "multiplicative";
   end
end

function value = identityCorrection(multiplicative)
   %IDENTITYCORRECTION The no-op correction for each mode.
   value = 0;
   if multiplicative
      value = 1;
   end
end

function value = oneCorrection(obs, model, multiplicative)
   %ONECORRECTION Robust bias (or ratio) of one overlap subset.
   if multiplicative
      value = median(obs ./ model, 'omitnan');
   else
      value = median(obs - model, 'omitnan');
   end
end
