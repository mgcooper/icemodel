function [state, alive, skipped] = reuseDatasetFamilyCases( ...
      manifest_file, requested_ids, prototype, kwargs)
   %REUSEDATASETFAMILYCASES Load staged cases for forcing-only attachment.
   %
   %  [state, alive, skipped] = ...
   %     icemodel.verification.setup.reuseDatasetFamilyCases( ...
   %     manifest_file, requested_ids, prototype, ...
   %     forcing_sources=["mar","merra"], coverage=coverage)
   %
   % This is the guarded build_observations=false path shared by dataset-family
   % importers. It requires every requested case to exist, preserves the decoded
   % case entry and colocation graph, and creates only the state needed to attach
   % explicitly requested RCM forcing. An explicit observation window must fit
   % inside the staged case period. forcing_startdate/forcing_enddate may select
   % an independent forcing window (SUMup years); otherwise the observation
   % request or staged case period supplies the forcing probe. A forcing probe
   % never relabels the staged observation period held in state.

   arguments
      manifest_file (1, 1) string
      requested_ids (1, :) string
      prototype (1, 1) struct
      kwargs.forcing_sources (1, :) string ...
         {icemodel.verification.validators.mustBeRcmSourceSelection( ...
         kwargs.forcing_sources)} = strings(1, 0)
      kwargs.coverage (1, 1) struct = struct()
      kwargs.startdate = ""
      kwargs.enddate = ""
      kwargs.forcing_startdate = ""
      kwargs.forcing_enddate = ""
   end

   % Normalize both optional windows before any manifest existence/read check.
   % Malformed public input therefore wins without touching staged state.
   [request_start, request_end, has_request] = ...
      icemodel.internal.pairedWindow( ...
      kwargs.startdate, kwargs.enddate);
   [forcing_start, forcing_end, has_forcing] = ...
      icemodel.internal.pairedWindow( ...
      kwargs.forcing_startdate, kwargs.forcing_enddate);

   % The fast path cannot infer cases from raw-source catalogs: its durable
   % contract is the already-staged family manifest.
   if ~isfile(manifest_file)
      error( ...
         'icemodel:verification:reuseDatasetFamilyCases:missingManifest', ...
         'build_observations=false requires an existing manifest: %s', ...
         manifest_file)
   end
   manifest = jsondecode(fileread(manifest_file));
   if ~isfield(manifest, 'cases') || isempty(manifest.cases)
      cases = struct([]);
      case_ids = strings(1, 0);
   else
      cases = manifest.cases;
      case_ids = string({cases.case_id});
   end

   requested_ids = reshape(requested_ids, 1, []);
   state = repmat(prototype, 1, numel(requested_ids));

   for k = 1:numel(requested_ids)
      hit = find(case_ids == requested_ids(k), 1);
      if isempty(hit)
         error('icemodel:verification:reuseDatasetFamilyCases:missingCase', ...
            ['build_observations=false requires existing case "%s" in ' ...
            '%s'], requested_ids(k), manifest_file)
      end
      entry = cases(hit);
      [period_start, period_end] = ...
         icemodel.verification.setup.periodBounds(entry.period);

      % A narrower forcing-only request may reuse a broader observation bundle,
      % but it may not advertise observations outside the staged case period.
      if has_request
         if isnat(period_start) || isnat(period_end)
            error( ...
               'icemodel:verification:reuseDatasetFamilyCases:unboundedPeriod', ...
               ['build_observations=false requires a bounded period for case ' ...
               '"%s" when a window is requested'], requested_ids(k))
         end
         if request_start < period_start || request_end > period_end
            error( ...
               'icemodel:verification:reuseDatasetFamilyCases:windowConflict', ...
               ['requested window %s to %s is outside existing case "%s" ' ...
               'period %s to %s'], string(request_start), string(request_end), ...
               requested_ids(k), string(period_start), string(period_end))
         end
      end

      % SUMup can request an RCM year span independent of an unbounded
      % observation period. Other importers use the explicit observation window
      % when present, then fall back to the staged case period.
      if has_forcing
         leg_start = forcing_start;
         leg_end = forcing_end;
      elseif has_request
         leg_start = request_start;
         leg_end = request_end;
      else
         leg_start = period_start;
         leg_end = period_end;
      end
      leg = struct();
      if ~isempty(kwargs.forcing_sources)
         leg = icemodel.verification.setup.resolveLegWindows( ...
            kwargs.forcing_sources, kwargs.coverage, leg_start, leg_end);
      end

      % Preserve the prior entry verbatim. Family-specific callbacks update only
      % colocation and its derived source lists after requested forcing stages.
      state(k).case_id = string(entry.case_id);
      state(k).alias = string(entry.case_id);
      state(k).point = [entry.site_location.lat_wgs84, ...
         entry.site_location.lon_wgs84];
      state(k).entry = entry;
      state(k).colocation = entry.colocation;
      state(k).leg = leg;
      % The case period belongs to the staged observations. Keep it stable while
      % the independently resolved leg carries any narrower or wider forcing
      % request used by the RCM builder.
      state(k).period = entry.period;
      state(k).reuse_entry = true;
      state(k).dry_run = false;
      if isfield(state(k), 'site_location')
         state(k).site_location = entry.site_location;
      end
   end

   alive = true(1, numel(state));
   skipped = struct('site', {}, 'reason', {});
end
