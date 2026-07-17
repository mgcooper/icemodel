function state = stageDatasetRcmForcing(state, alive, kwargs)
   %STAGEDATASETRCMFORCING Delegate one RCM source at a time to stageRcmForcing.
   %
   %  state = icemodel.verification.setup.stageDatasetRcmForcing(state, alive, ...)
   %
   % Shared orchestration for family importers that already staged their native
   % observations. It builds a per-source legspec through a family-supplied
   % callback, calls stageRcmForcing for exactly one source at a time, merges the
   % returned colocation legs into each state record, and invokes an optional
   % persist callback after every source for kill-safe manifests.
   % forcing_sources selects only sources requested by this call; duplicate
   % tokens collapse in stable order and omitted fields remain unchanged.
   % method is restricted to the builders' supported nearest/natural choices.
   % Model met defaults to dt_out="15m"; pass dt_out="" for native cadence.
   % Data/userdata defaults to hourly at the shared writer boundary.
   %
   % See also: icemodel.verification.setup.stageRcmForcing,
   %  icemodel.verification.setup.mergeColocation

   arguments
      state (1, :) struct
      alive (1, :) logical
      kwargs.dataset_family (1, 1) string = ""
      kwargs.forcing_sources (1, :) string ...
         {icemodel.verification.validators.mustBeRcmSourceSelection( ...
         kwargs.forcing_sources)} = ...
         icemodel.verification.namelists.rcmsources()
      kwargs.leg_callback = []
      kwargs.after_source_callback = []
      kwargs.persist_callback = []
      kwargs.met_outdir (1, 1) string = ""
      kwargs.userdata_outdir (1, 1) string = ""
      kwargs.mar_dir (1, 1) string = ""
      kwargs.merra_dir (1, 1) string = ""
      kwargs.racmo_dir (1, 1) string = ""
      kwargs.modis_dir (1, 1) string = ""
      kwargs.method (1, 1) string ...
         {mustBeMember(kwargs.method, ["nearest", "natural"])} = "nearest"
      kwargs.dt_out (1, 1) string ...
         {mustBeMember(kwargs.dt_out, ["", "15m"])} = "15m"
      kwargs.overwrite (1, 1) logical = false
   end

   % Treat a caller's selection as an ordered set so duplicate tokens cannot
   % repeat expensive builders, warnings, or persistence callbacks.
   forcing_sources = unique(kwargs.forcing_sources, 'stable');
   alive_idx = find(alive);
   if isempty(alive_idx) || isempty(forcing_sources)
      return
   end

   % A callback is required only when at least one live source will be staged.
   if isempty(kwargs.leg_callback)
      error('icemodel:verification:stageDatasetRcmForcing:missingCallback', ...
         'leg_callback is required')
   end

   % Resolve writer defaults once so discovery and prior fallback inspect the
   % same roots even when the caller omits explicit output directories.
   [met_outdir, userdata_outdir] = ...
      icemodel.verification.setup.rcmArtifactOutputDirs( ...
      kwargs.met_outdir, kwargs.userdata_outdir);

   points = vertcat(state(alive_idx).point);
   for src = reshape(forcing_sources, 1, [])
      fprintf('[staging] %s forcing for %d %s case(s)...\n', ...
         upper(char(src)), numel(alive_idx), kwargs.dataset_family);

      legspec = repmat(legProto(src), 1, numel(alive_idx));
      for j = 1:numel(alive_idx)
         s = state(alive_idx(j));
         % A manifest id can name a different physical cache point in another
         % family. Keep that scientific id in the manifest while letting the
         % importer provide a collision-safe RCM-only storage identity.
         legspec(j).alias = string(s.alias);
         if isfield(s, 'storage_alias') ...
               && strlength(string(s.storage_alias)) > 0
            legspec(j).alias = string(s.storage_alias);
         end
         legspec(j).(char(src)) = legWithStateDiscoveryWindow( ...
            kwargs.leg_callback(s, src), s);
      end

      colocation = icemodel.verification.setup.stageRcmForcing(points, ...
         legspec=legspec, forcing_sources=src, ...
         met_outdir=met_outdir, userdata_outdir=userdata_outdir, ...
         mar_dir=kwargs.mar_dir, merra_dir=kwargs.merra_dir, ...
         racmo_dir=kwargs.racmo_dir, modis_dir=kwargs.modis_dir, ...
         method=kwargs.method, dt_out=kwargs.dt_out, ...
         overwrite=kwargs.overwrite);

      n_staged = 0;
      for j = 1:numel(alive_idx)
         idx = alive_idx(j);
         period = stateDiscoveryPeriod(state(idx));
         if isempty(period)
            period = struct('start', '', 'end', '');
         end
         colocation{j} = icemodel.verification.setup.preserveRcmLegs( ...
            state(idx).colocation, colocation{j}, src, period, ...
            met_outdir=met_outdir, userdata_outdir=userdata_outdir, ...
            method=kwargs.method, point=state(idx).point, ...
            overwrite=kwargs.overwrite);
         state(idx).colocation = ...
            icemodel.verification.setup.mergeColocation( ...
            state(idx).colocation, colocation{j});
         if isfield(state(idx).colocation, char(src)) ...
               && state(idx).colocation.(char(src)).staged
            n_staged = n_staged + 1;
         end
      end
      fprintf('[staging] %s: %d staged, %d skipped\n', upper(char(src)), ...
         n_staged, numel(alive_idx) - n_staged);

      % Family-specific optional products may attach to the completed source
      % leg, but cache discovery and primary RCM writes stay owned above.
      if ~isempty(kwargs.after_source_callback)
         state = kwargs.after_source_callback(state, alive_idx, src);
      end

      if ~isempty(kwargs.persist_callback)
         kwargs.persist_callback(state);
      end
   end
end

function L = legWithStateDiscoveryWindow(L, state)
   %LEGWITHSTATEDISCOVERYWINDOW Search the full case beyond a clipped raw leg.
   period = stateDiscoveryPeriod(state);
   if isempty(period) || strlength(string(period.start)) == 0 ...
         || strlength(string(period.end)) == 0
      return
   end
   L.discovery_start = icemodel.verification.setup.ensureUtc(period.start);
   L.discovery_end = icemodel.verification.setup.ensureUtc(period.end);
end

function period = stateDiscoveryPeriod(state)
   %STATEDISCOVERYPERIOD Return case period for cached-artifact discovery.
   period = [];
   if isfield(state, 'period') && isstruct(state.period) ...
         && isfield(state.period, 'start') && isfield(state.period, 'end')
      period = state.period;
   elseif isfield(state, 'entry') && isstruct(state.entry) ...
         && isfield(state.entry, 'period') && isstruct(state.entry.period) ...
         && isfield(state.entry.period, 'start') ...
         && isfield(state.entry.period, 'end')
      period = state.entry.period;
   end
end

function L = legProto(sources)
   %LEGPROTO Prototype legspec element for one or more RCM source fields.
   L = struct('alias', "");
   proto = struct('staged', false, 'years', [], 'start', NaT, 'end', NaT, ...
      'reason', "");
   for src = sources
      L.(char(src)) = proto;
   end
end
