function [estimate, clamped] = applyProxyCalibration( ...
      calibration, times, x_model, kwargs)
   %APPLYPROXYCALIBRATION Apply a fitted proxy calibration to model samples.
   %
   %  estimate = icemodel.forcing.reconstruct.applyProxyCalibration( ...
   %     calibration, times, x_model)
   %
   % Role
   %  Application half of the proxy-calibration tier: evaluates the
   %  per-season correction fitProxyCalibration estimated. Missing model
   %  samples stay missing; the identity calibration passes the model
   %  through unchanged (recorded, per POLICY B5). Version-2 swd records
   %  carry D-28 elevation-binned ratios that apply when the caller
   %  supplies the matching solar elevation.
   %
   % Name-value
   %  target_elevation : target-station solar elevation (degrees) per
   %     sample. Required to use the bands of a binned swd record. If you
   %     omit it, the per-season scalar applies. Callers without station
   %     geometry omit it.
   %
   % See also: icemodel.forcing.reconstruct.fitProxyCalibration,
   %  icemodel.forcing.reconstruct.solarElevationBands

   arguments
      calibration (1, 1) struct
      times datetime
      x_model (:, 1) double
      kwargs.target_elevation (:, 1) double = zeros(0, 1)
   end

   % A supplied elevation vector must cover the axis exactly; a length
   % mismatch would disable the binned correction.
   if ~isempty(kwargs.target_elevation) ...
         && numel(kwargs.target_elevation) ~= numel(x_model)
      error('icemodel:reconstruct:applyProxyCalibration:targetElevationSize', ...
         ['target_elevation must contain one solar elevation value per ' ...
         'model sample when supplied']);
   end

   % Schema guard (D-28): only version-2 swd records carry elevation-binned
   % ratios, detected by field presence rather than a version comparison.
   % A single-ratio record, or a binned record applied by a caller without
   % station geometry, falls back to the per-season scalar.
   use_bins = isfield(calibration, 'binned_corrections') ...
      && ~isempty(kwargs.target_elevation);
   if use_bins
      bin_index = discretize(kwargs.target_elevation, ...
         calibration.bin_edges_deg);
   end

   season = icemodel.forcing.reconstruct.seasonOf(times);
   estimate = nan(numel(x_model), 1);
   for name = ["DJF", "MAM", "JJA", "SON"]
      in_season = season == name & isfinite(x_model);
      if ~any(in_season)
         continue
      end

      % A saved or hand-edited record can carry a nonfinite correction. That
      % correction turns finite proxy input into Inf or NaN, and the later
      % validity checks then reject the samples as if the PROXY were
      % unusable. Skip the season instead, so those samples stay missing and
      % the denial names the real cause.
      if ~use_bins && ~isfinite(calibration.corrections.(char(name)))
         continue
      end
      if use_bins
         % Each sample scales by its own elevation band's ratio, so the
         % twilight and shoulder regimes get their fitted corrections
         % instead of one season-wide multiplier (D-28).
         ratios = calibration.binned_corrections.(char(name));
         members = find(in_season);
         estimate(members) = x_model(members) ...
            .* reshape(ratios(bin_index(members)), [], 1);
      elseif calibration.mode == "multiplicative"
         estimate(in_season) = x_model(in_season) ...
            * calibration.corrections.(char(name));
      else
         estimate(in_season) = x_model(in_season) ...
            + calibration.corrections.(char(name));
      end
   end

   % D-27 (user ruling 2026-07-27): a correction that pushes rh past its
   % physical bounds is calibration arithmetic, not physics. A
   % near-saturation source and a positive ratio can exceed 100%. Without a
   % clamp the check rejects the candidate (SWC lost 4.5% of rh). The clamp
   % here covers every consumer, both the method tiers and the last resort.
   % D-51 (2026-08-05) extends the same rule to wspd: a wind ratio below one
   % can move a valid 0.1 m/s proxy posting below the runtime floor (TAS_L
   % 2010 lost four samples that way). The second output reports every
   % clamped sample.
   clamped = false(size(estimate));
   bounded_calibration = isfield(calibration, 'channel') ...
      && ismember(string(calibration.channel), ["rh", "wspd"]);
   if bounded_calibration
      bounds = icemodel.forcing.reconstruct.physicalBounds( ...
         string(calibration.channel));
      finite = isfinite(estimate);
      clamped = finite ...
         & (estimate < bounds(1) | estimate > bounds(2));
      estimate(finite) = min(max(estimate(finite), bounds(1)), bounds(2));
   end
end
