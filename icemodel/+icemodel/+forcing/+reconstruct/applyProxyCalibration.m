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
   %     sample. Required to honor a binned swd record's bands; omitted,
   %     the per-season scalar applies for callers without station geometry.
   %
   % See also: icemodel.forcing.reconstruct.fitProxyCalibration,
   %  icemodel.forcing.reconstruct.solarElevationBands

   arguments
      calibration (1, 1) struct
      times datetime
      x_model (:, 1) double
      kwargs.target_elevation (:, 1) double = zeros(0, 1)
   end

   % A supplied elevation vector must cover the axis exactly; a silent
   % length mismatch would quietly disable the binned correction.
   if ~isempty(kwargs.target_elevation) ...
         && numel(kwargs.target_elevation) ~= numel(x_model)
      error('icemodel:reconstruct:applyProxyCalibration:targetElevationSize', ...
         ['target_elevation must contain one solar elevation value per ' ...
         'model sample when supplied']);
   end

   % Backward-compatible schema guard (D-28): only version-2 swd records
   % carry elevation-binned ratios, detected by field presence rather
   % than a version comparison so any legacy single-ratio record — and
   % any binned record applied by an elevation-less caller — falls back
   % to the per-season scalar exactly as before.
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
   % physical bounds is calibration arithmetic, not physics — near-saturation
   % sources plus a positive ratio exceed 100% and previously refused
   % adoption (SWC lost 4.5% of rh). Clamping ONCE here covers every
   % consumer (method tier and last resort alike).
   clamped = false(size(estimate));
   if isfield(calibration, 'channel') && string(calibration.channel) == "rh"
      bounds = icemodel.forcing.reconstruct.physicalBounds("rh");
      finite = isfinite(estimate);
      clamped = finite ...
         & (estimate < bounds(1) | estimate > bounds(2));
      estimate(finite) = min(max(estimate(finite), bounds(1)), bounds(2));
   end
end
