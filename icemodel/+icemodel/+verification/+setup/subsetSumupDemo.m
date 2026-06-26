function [observations, provenance] = subsetSumupDemo(observations, point, kwargs)
   %SUBSETSUMUPDEMO Reduce a SUMup obs bundle to the committed demo fixture.
   %
   %  [observations, provenance] = ...
   %     icemodel.verification.setup.subsetSumupDemo(observations, [lat lon])
   %  [observations, provenance] = ...
   %     icemodel.verification.setup.subsetSumupDemo(observations, [lat lon], ...
   %        max_rows=500)
   %
   %  Derives the small committed demo fixture from the full per-point SUMup
   %  observation bundle (icemodel.verification.setup.buildSumupObservations)
   %  by an EXPLICIT, REPRODUCIBLE, two-step rule so the committed fixture is
   %  not opaque (the prior 200-row tables had no recorded derivation):
   %
   %    MINIMAL-FIXTURE RULE
   %      1. SINGLE NEAREST PROFILE - per sub-table, select the `name_key`
   %         whose mean record location is closest to the case point
   %         (EPSG:3413 distance) and drop every other profile. This keeps one
   %         self-consistent core/station series per channel.
   %      2. HEAD CAP - keep the first `max_rows` rows of that profile in
   %         native record order (time-ascending for the temperature timetable;
   %         depth-ascending for the density profile). A single station's
   %         thermistor series spans many years of hourly samples, so the cap
   %         bounds the committed fixture to a small, contiguous slice.
   %
   %    Both steps are deterministic, so the fixture is fully reproducible from
   %    the full staged bundle + the point + max_rows.
   %
   %  Inputs
   %    observations : struct  full SUMup obs bundle from buildSumupObservations
   %                   (format="subsurface_profile_bundle"; density TABLE,
   %                   subsurface_temperature TIMETABLE, accumulation TABLE).
   %    point        : [lat lon] WGS84 query point.
   %
   %  Name-value
   %    max_rows : double (default 500)  per-channel head cap (step 2).
   %
   %  Returns
   %    observations : struct  same bundle, each sub-table reduced to the
   %                   single nearest profile capped at max_rows rows.
   %    provenance   : struct  records the rule and, per channel, the selected
   %                   name_key, the resulting row count, and the datetime span,
   %                   so the manifest documents exactly how the fixture was
   %                   derived.
   %
   % See also: icemodel.verification.setup.buildSumupObservations,
   %  icemodel.verification.setup.importSumup

   arguments
      observations (1, 1) struct
      point (1, 2) double
      kwargs.max_rows (1, 1) double {mustBePositive} = 500
   end

   proj = icemodel.forcing.helpers.psnProjection();
   [px, py] = projfwd(proj, point(1), point(2));

   provenance = struct( ...
      'rule', ['(1) single nearest profile per channel (minimum EPSG:3413 ' ...
      'distance of the profile mean location to the case point); ' ...
      '(2) first max_rows rows of that profile in native record order'], ...
      'max_rows', kwargs.max_rows, ...
      'point_lat_wgs84', point(1), ...
      'point_lon_wgs84', point(2), ...
      'channels', struct([]));

   channels = ["density", "subsurface_temperature", "smb"];
   for k = 1:numel(channels)
      ch = channels(k);
      if ~isfield(observations, ch) || isempty(observations.(ch))
         continue
      end
      [observations.(ch), rec] = nearestProfile( ...
         observations.(ch), proj, px, py, ch, kwargs.max_rows);
      if isempty(provenance.channels)
         provenance.channels = rec;
      else
         provenance.channels(end + 1) = rec;
      end
   end
end

function [reduced, rec] = nearestProfile(tbl, proj, px, py, channel, max_rows)
   %NEARESTPROFILE Keep first max_rows rows of the single nearest name_key.
   istt = istimetable(tbl);
   keys = tbl.name_key;
   uniq = unique(keys);

   % Mean projected location of each candidate profile, then nearest.
   d = zeros(numel(uniq), 1);
   for i = 1:numel(uniq)
      m = keys == uniq(i);
      [rx, ry] = projfwd(proj, mean(tbl.latitude(m)), mean(tbl.longitude(m)));
      d(i) = hypot(rx - px, ry - py);
   end
   [~, best] = min(d);
   keep = keys == uniq(best);
   reduced = tbl(keep, :);

   % Step 2: head cap to the first max_rows native-order rows.
   if height(reduced) > max_rows
      reduced = reduced(1:max_rows, :);
   end

   if istt
      span = [min(reduced.Time), max(reduced.Time)];
   elseif ismember("datetime", string(reduced.Properties.VariableNames))
      span = [min(reduced.datetime), max(reduced.datetime)];
   elseif ismember("start_date", string(reduced.Properties.VariableNames))
      span = [min(reduced.start_date), max(reduced.end_date)];
   else
      span = [NaT, NaT];
   end

   rec = struct( ...
      'channel', char(channel), ...
      'name_key', uniq(best), ...
      'name', char(firstName(reduced)), ...
      'n_rows', height(reduced), ...
      'datetime_start', char(string(span(1))), ...
      'datetime_end', char(string(span(2))));
end

function nm = firstName(tbl)
   %FIRSTNAME First resolved profile name (already deblanked upstream).
   if ismember("name", string(tbl.Properties.VariableNames)) && height(tbl) > 0
      nm = string(tbl.name(1));
   else
      nm = "";
   end
end
