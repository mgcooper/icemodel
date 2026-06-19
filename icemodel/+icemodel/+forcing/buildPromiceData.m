function [Data, metadata] = buildPromiceData(site, kwargs)
   %BUILDPROMICEDATA Build a PROMICE evaluation/userdata Data timetable.
   %
   %  [Data, metadata] = icemodel.forcing.buildPromiceData(site)
   %  [Data, metadata] = ... buildPromiceData(site, source_dir=..., ...
   %     startdate=..., enddate=..., frequency="daily")
   %
   % Reads the station's PROMICE v3 hourly NetCDF and assembles the
   % observational channels used to evaluate icemodel output and to feed
   % the met-swap (userdata) mechanism: the surface energy balance
   % terms, derived net fluxes, cumulative ice ablation from the
   % pressure-transducer record, snow depth derived from the sonic boom
   % height, and the ice-temperature string. Location metadata attaches
   % as table CustomProperties so the result can be written with
   % icemodel.forcing.helpers.writeuserdata.
   %
   % Ablation: the corrected pressure-transducer depth is despiked with
   % site-specific correction recipes (service-visit windows curated for
   % KAN_M and KAN_L from the legacy processing; other stations use the
   % record as-is) and converted to cumulative ablation, in meters of
   % surface lowering, positive downward, zeroed at the first valid
   % sample of the requested window.
   %
   % Snow depth: per accumulation year (Sep 1 to Aug 31), the
   % late-summer surface standoff is the September median boom height;
   % snow depth is that reference minus the boom height, clamped at
   % zero. This is a first-order method (station maintenance and summer
   % melt mix into the boom record); treat the channel as approximate.
   %
   % Gap policy: observational channels are NOT gap-filled (missing data
   % stays missing so evaluations are honest); the physical-range clamps
   % of metchecks are still applied.
   %
   % Inputs
   %  site - station name ("KAN_M" or compact alias "kanm")
   %
   % Name-value
   %  source_dir : PROMICE NetCDF directory (see readPromiceAws)
   %  startdate, enddate : optional window; default full station record
   %  frequency : "hourly" (default) or "daily" (daily means)
   %
   % Outputs
   %  Data     - timetable with CustomProperties (X, Y, Lat, Lon, Elev,
   %             Slope, ScalarUnits). Channels: tair, tsfc [K]; swd, swu,
   %             lwd, lwu, swn, lwn, netr, shf, lhf, thf [W m-2]; albedo,
   %             cfrac [-]; rh [%]; wspd [m s-1]; wdir [deg]; psfc [Pa];
   %             ablation, snow_depth, boom_height, transducer_depth [m];
   %             tice1..tice8 [K].
   %  metadata - provenance: source file, station, correction recipe
   %             applied, QA/QC gap counts
   %
   % See also: icemodel.forcing.readPromiceAws,
   %  icemodel.forcing.buildPromiceMet,
   %  icemodel.forcing.helpers.writeuserdata

   arguments
      site (1, 1) string
      kwargs.source_dir (1, 1) string = ""
      kwargs.startdate = ""
      kwargs.enddate = ""
      kwargs.frequency (1, 1) string ...
         {mustBeMember(kwargs.frequency, ["hourly", "daily"])} = "hourly"
   end

   [aws, source_meta] = icemodel.forcing.readPromiceAws(site, ...
      source_dir=kwargs.source_dir, timescale="hourly", ...
      startdate=kwargs.startdate, enddate=kwargs.enddate);

   % Derived net fluxes.
   aws.swn = aws.swd - aws.swu;
   aws.lwn = aws.lwd - aws.lwu;
   aws.netr = aws.swn + aws.lwn;
   aws.thf = aws.shf + aws.lhf;

   % Cumulative ice ablation from the despiked transducer record.
   [aws.ablation, recipe_applied] = cumulativeAblation( ...
      aws.transducer_depth, aws.Time, site);

   % Snow depth from the sonic boom height.
   aws.snow_depth = snowDepthFromBoom(aws.boom_height, aws.Time);

   % Order the output channels.
   channels = ["tair", "tsfc", "swd", "swu", "lwd", "lwu", "swn", ...
      "lwn", "netr", "shf", "lhf", "thf", "albedo", "cfrac", "rh", ...
      "wspd", "wdir", "psfc", "ablation", "snow_depth", "boom_height", ...
      "transducer_depth", "tice1", "tice2", "tice3", "tice4", "tice5", ...
      "tice6", "tice7", "tice8"];
   channels = channels(ismember(channels, ...
      string(aws.Properties.VariableNames)));
   Data = aws(:, cellstr(channels));

   % Physical-range clamps only; observational gaps stay missing.
   [Data, checks] = icemodel.forcing.helpers.metchecks(Data, ...
      fillgaps=false);

   if kwargs.frequency == "daily"
      Data = retime(Data, 'daily', 'mean');
   end

   % Location metadata for the Data-file (userdata) contract.
   proj = icemodel.forcing.helpers.psnProjection();
   [x, y] = projfwd(proj, source_meta.lat, source_meta.lon);
   Data = addprop(Data, ...
      {'X', 'Y', 'Lat', 'Lon', 'Elev', 'Slope', 'ScalarUnits'}, ...
      repmat({'table'}, 1, 7));
   Data.Properties.CustomProperties.X = x;
   Data.Properties.CustomProperties.Y = y;
   Data.Properties.CustomProperties.Lat = source_meta.lat;
   Data.Properties.CustomProperties.Lon = source_meta.lon;
   Data.Properties.CustomProperties.Elev = source_meta.elev;
   Data.Properties.CustomProperties.Slope = NaN;
   Data.Properties.CustomProperties.ScalarUnits = ...
      ["m", "m", "degrees N", "degrees W", "m asl", "m/m"];

   metadata = source_meta;
   metadata.checks = checks;
   metadata.frequency = kwargs.frequency;
   metadata.ablation_recipe = recipe_applied;
   metadata.snow_depth_method = ...
      "September-median boom-height reference minus boom height, clamped >= 0";
   metadata.gap_policy = "no gap fill (observational); clamps applied";
end

%% Local functions
function [ablation, recipe_name] = cumulativeAblation(depth, Time, site)
   %CUMULATIVEABLATION Despiked cumulative surface ablation [m].
   %
   % The corrected transducer depth DECREASES as the surface ablates
   % (the sensor hangs in a drilled hole; surface lowering shortens the
   % column above it - verified on the v3 record: KAN_L JJA 2010 trend
   % -4.3 m). Station service re-drills the assembly deeper, which
   % appears as a large POSITIVE depth jump. Real ablation can only
   % decrease the depth, so the pipeline is:
   %
   %  1. Apply the site recipe: trim the unusable record start and
   %     blank the curated noisy service windows (KAN_M / KAN_L windows
   %     from the legacy GREENLAND get_PROMICE_ablation processing).
   %  2. Remove every positive jump larger than the reset threshold
   %     between consecutive finite samples (service resets at any
   %     date, on any station - this generalizes the legacy per-event
   %     step-shift bookkeeping).
   %  3. ablation = -(depth - first finite depth), positive = surface
   %     lowering (the legacy sign convention).
   [recipe, recipe_name] = ablationRecipe(site);
   for n = 1:numel(recipe)
      step = recipe(n);
      switch step.op
         case "trim_before"
            depth(Time < step.t1) = NaN;
         case "blank"
            depth(Time > step.t1 & Time < step.t2) = NaN;
      end
   end

   depth = removeResets(depth, 0.5);   % resets are meters; melt is mm/h

   first = find(isfinite(depth), 1);
   if isempty(first)
      ablation = nan(size(depth));
   else
      ablation = -(depth - depth(first));
   end
end

function depth = removeResets(depth, threshold)
   %REMOVERESETS Remove positive level jumps (service re-drills).
   finite = find(isfinite(depth));
   jumps = diff(depth(finite));
   resets = jumps > threshold;
   offsets = cumsum(jumps .* resets);
   depth(finite(2:end)) = depth(finite(2:end)) - offsets;
end

function [recipe, recipe_name] = ablationRecipe(site)
   %ABLATIONRECIPE Site-keyed curated correction windows.
   %
   % The KAN_M and KAN_L windows are the curated service-visit/spike
   % windows from the legacy GREENLAND get_PROMICE_ablation processing
   % (timestamps preserved verbatim); resets themselves are handled
   % generically by removeResets. Stations without a recipe use the
   % corrected transducer record as-is.
   step = @(op, t1, t2) struct('op', string(op), 't1', utc(t1), 't2', utc(t2));
   switch lower(erase(site, "_"))
      case 'kanm'
         recipe_name = "kanm legacy service windows + generic reset removal";
         recipe = [ ...
            step("trim_before", datetime(2012, 5, 4), NaT)
            step("blank", datetime(2014, 3, 13, 14, 0, 0), datetime(2014, 3, 14, 16, 0, 0))
            step("blank", datetime(2014, 8, 3, 12, 0, 0), datetime(2014, 8, 3, 21, 0, 0))
            step("blank", datetime(2014, 9, 16, 5, 0, 0), datetime(2014, 10, 6, 11, 0, 0))
            step("blank", datetime(2015, 2, 25, 7, 0, 0), datetime(2015, 2, 26, 8, 0, 0))
            step("blank", datetime(2015, 6, 22, 22, 0, 0), datetime(2015, 7, 11, 19, 0, 0))
            ];
      case 'kanl'
         recipe_name = "kanl legacy service windows + generic reset removal";
         recipe = [ ...
            step("blank", datetime(2009, 9, 22), datetime(2009, 10, 31))
            step("blank", datetime(2010, 10, 2), datetime(2010, 10, 13))
            step("blank", datetime(2014, 10, 14), datetime(2014, 11, 29))
            step("blank", datetime(2015, 10, 4), datetime(2015, 10, 7))
            step("blank", datetime(2016, 9, 27), datetime(2016, 10, 13))
            step("blank", datetime(2011, 6, 4, 15, 0, 0), datetime(2011, 6, 4, 19, 0, 0))
            step("blank", datetime(2012, 8, 22, 12, 0, 0), datetime(2012, 8, 22, 22, 0, 0))
            step("blank", datetime(2012, 10, 7, 22, 0, 0), datetime(2012, 10, 25, 17, 0, 0))
            step("blank", datetime(2015, 4, 28, 1, 0, 0), datetime(2015, 4, 28, 15, 0, 0))
            step("blank", datetime(2015, 7, 6, 17, 0, 0), datetime(2015, 7, 6, 20, 0, 0))
            step("blank", datetime(2016, 7, 16, 15, 0, 0), datetime(2016, 7, 16, 18, 0, 0))
            step("blank", datetime(2017, 9, 1, 17, 0, 0), datetime(2017, 9, 2, 9, 0, 0))
            step("blank", datetime(2018, 8, 28, 13, 0, 0), datetime(2018, 8, 29, 5, 0, 0))
            ];
      otherwise
         recipe_name = "generic reset removal only";
         recipe = struct('op', {}, 't1', {}, 't2', {});
   end
end

function t = utc(t)
   %UTC Tag a timezone-naive datetime as UTC (PROMICE time base).
   if ~isnat(t)
      t.TimeZone = 'UTC';
   end
end

function snow_depth = snowDepthFromBoom(boom, Time)
   %SNOWDEPTHFROMBOOM First-order snow depth from sonic boom height [m].
   snow_depth = nan(size(boom));
   for yyyy = unique(year(Time))'
      ref_window = Time >= datetime(yyyy, 9, 1, 'TimeZone', 'UTC') & ...
         Time < datetime(yyyy, 10, 1, 'TimeZone', 'UTC');
      ref = median(boom(ref_window), 'omitnan');
      if ~isfinite(ref)
         continue
      end
      season = Time >= datetime(yyyy, 9, 1, 'TimeZone', 'UTC') & ...
         Time < datetime(yyyy + 1, 9, 1, 'TimeZone', 'UTC');
      snow_depth(season) = max(0, ref - boom(season));
   end
end
