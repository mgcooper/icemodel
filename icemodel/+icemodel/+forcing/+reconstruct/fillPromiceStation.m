function result = fillPromiceStation(site, kwargs)
   %FILLPROMICESTATION Produce the gap-filled met product for one station.
   %
   %  result = icemodel.forcing.reconstruct.fillPromiceStation("kanm")
   %
   % Role
   %  The per-station production driver for the <family>_filled product —
   %  by default promice_filled, the canonical runnable PROMICE forcing
   %  (native promice met is
   %  incomplete for most station-years and cannot force the model; it
   %  is retained unmodified as the provenance source). The driver loads
   %  the staged native met, assembles the donor pool (nearby PROMICE
   %  stations, staged K-transect cases, GC-Net origin-observed
   %  samples) and the calibrated RCM proxies, runs the
   %  stationMethodPlan selection experiment, composes the admitted
   %  methods with reconstructSeries, adopts proxy precipitation per the
   %  approved policy (PROMICE supplies no total or solid precipitation;
   %  corrected liquid observations remain native constraints — proxy total
   %  is a recorded policy adoption, not an admission), adopts aligned proxy
   %  values for any residual missing
   %  required-channel samples as the audited last resort, and writes
   %  the met_<site>_<family>_filled artifact with inline per-sample
   %  provenance channels plus the plan, audit, and readiness records.
   %
   % Name-value
   %  opts : reconstruction options struct from
   %     icemodel.forcing.reconstruct.setopts — the single source of the
   %     channel lists, thresholds, proxy-source order, and seed this
   %     driver consumes.
   %  family : staged native product family token (default "promice").
   %     Derives the met_dir default data/input/met/<family>, the native
   %     target grammar met_<site>_<family>_*_15m.mat, the producer
   %     manifest data/eval/<family>/manifest.json (colocation.<family>
   %     leg), and the shipped product token <family>_filled. The
   %     promice-only builder machinery — raw-source replay, raw-fallback
   %     shortwave QC, and the winter-albedo stamp mask — runs only for
   %     family "promice" (bead icemodel-g1n.49); station donors always
   %     draw from the promice donor family.
   %  met_dir : staged native met directory (default
   %     data/input/met/<family> under the repo). Proxies are read from
   %     the sibling canonical staged directories named by the proxy
   %     catalog (widest 15-minute window file per site) — no side cache
   %     exists. K-transect, GC-Net, and the donor inventory resolve from
   %     that same selected data root; no repository-root fallback occurs.
   %  modis_dir : staged MODIS daily-albedo userdata directory used by
   %     the A11/B12 albedo tier (default data/input/userdata/modis).
   %  out_dir : filled-product met directory (default
   %     data/input/met/<family>_filled).
   %  qa_dir : plan/audit/ledger directory (default
   %     data/preview/qa/gapfill).
   %  donor_sites : PROMICE donor site tokens to consider (default "auto":
   %     every other staged station, geometry-gated inside the plan).
   %  use_ktransect, use_gcnet : donor-family toggles (default true).
   %  write : write artifacts (default true); false returns results only.
   %
   % Returns
   %  result : struct — site, plan, filled (timetable), provenance,
   %     audit, readiness (per-year table), seam_quality (SWD boundary-QA
   %     table, empty when SWD is absent), flat_run_findings (native
   %     buried/rimed-sensor exclusions), and met_file (written path or
   %     "").
   %
   % See also: icemodel.forcing.reconstruct.setopts,
   %  icemodel.forcing.reconstruct.stationMethodPlan,
   %  icemodel.forcing.reconstruct.reconstructSeries

   arguments
      site (1, 1) string ...
         {icemodel.forcing.reconstruct.mustBeStationToken}
      kwargs.opts (1, 1) struct = icemodel.forcing.reconstruct.setopts()
      kwargs.family (1, 1) string = "promice"
      kwargs.met_dir (1, 1) string = ""
      kwargs.modis_dir (1, 1) string = ""
      kwargs.out_dir (1, 1) string = ""
      kwargs.qa_dir (1, 1) string = ""
      kwargs.donor_sites (1, :) string = "auto"
      kwargs.use_ktransect (1, 1) logical = true
      kwargs.use_gcnet (1, 1) logical = true
      kwargs.write (1, 1) logical = true
   end
   opts = kwargs.opts;
   % The family token names the staged native product this driver fills;
   % every family-derived default below (met_dir, out_dir, manifest leg,
   % product token) flows from this single value (bead icemodel-g1n.49).
   family = kwargs.family;
   if ~(isscalar(kwargs.donor_sites) && kwargs.donor_sites == "auto")
      if any(kwargs.donor_sites == "auto")
         error('icemodel:reconstruct:mustBeStationToken:invalidToken', ...
            'auto must be the sole donor_sites sentinel');
      end
      icemodel.forcing.reconstruct.mustBeStationToken(kwargs.donor_sites);
   end

   repo = icemodel.internal.fullpath;
   met_dir = defaultPath(kwargs.met_dir, repo, ...
      fullfile('data', 'input', 'met', char(family)));
   out_dir = defaultPath(kwargs.out_dir, repo, ...
      fullfile('data', 'input', 'met', char(family + "_filled")));
   qa_dir = defaultPath(kwargs.qa_dir, repo, ...
      fullfile('data', 'preview', 'qa', 'gapfill'));
   modis_dir = defaultPath(kwargs.modis_dir, repo, ...
      fullfile('data', 'input', 'userdata', 'modis'));

   % Target: staged native met plus its point from the artifact metadata.
   % The winter-albedo mask marks the native builder's constant stamp as
   % missing so methods fill those samples with honest provenance.
   [staged_series, location, winter_albedo_mask, native_file, ...
      staged_provenance, flat_run_findings] = ...
      loadStationMet(met_dir, site, opts, family);
   [series, winter_albedo_mask, restore_15m, native_provenance] = ...
      reconstructionAxis(staged_series, winter_albedo_mask, ...
      staged_provenance);

    % Validate and select the acceptance-window proxy files first, then load
    % exactly that pinned set in catalog adoption-preference order. A missing
    % staged proxy is recorded, not fatal — other tiers still plan.
    try
       [acceptance_window, proxy_window_files] = ...
          icemodel.forcing.reconstruct.acceptanceWindow( ...
          site, met_dir=met_dir, location=location, opts=opts);
    catch exception
       % A malformed or internally disjoint proxy inventory invalidates the
       % old A6 product just as surely as an empty or record-disjoint window.
       if kwargs.write && startsWith(string(exception.identifier), ...
             "icemodel:reconstruct:acceptanceWindow:")
          retirePublishedArtifacts(site, out_dir, qa_dir, family);
       end
       rethrow(exception)
    end

   % POLICY A6/D-17: validated staged proxy met defines the complete product
   % span. An empty proxy inventory defines no product span; it is not
   % permission to reconstruct the full native record from donors or
   % climatology. Production refusals retire any stale prior publication.
   if any(isnat(acceptance_window))
      if kwargs.write
         retirePublishedArtifacts(site, out_dir, qa_dir, family);
      end
      error('icemodel:reconstruct:fillPromiceStation:noProxyWindow', ...
         ['no validated staged MAR/MERRA proxy window exists for %s; ' ...
         'no product span exists (POLICY A6/D-17)'], site);
   end
   staged_in = staged_series.Properties.RowTimes >= acceptance_window(1) ...
      & staged_series.Properties.RowTimes <= acceptance_window(2);
   if ~any(staged_in)
      if kwargs.write
         retirePublishedArtifacts(site, out_dir, qa_dir, family);
      end
      error( ...
         'icemodel:reconstruct:fillPromiceStation:windowRecordDisjoint', ...
         ['the %s record lies wholly outside its staged proxy ' ...
         'window; no product span exists (POLICY A6/D-17)'], site);
   end

   % A guarded quarter-hour product can begin inside its owning hourly or
   % half-hourly posting. Retain every source posting that owns at least one
   % delivered row, then clip the delivered axis to the exact proxy window.
   if restore_15m
      owner = postingLocations(series.Properties.RowTimes, ...
         staged_series.Properties.RowTimes(staged_in));
      in_window = false(height(series), 1);
      in_window(unique(owner)) = true;
   else
      in_window = series.Properties.RowTimes >= acceptance_window(1) ...
         & series.Properties.RowTimes <= acceptance_window(2);
   end
   series = series(in_window, :);
   winter_albedo_mask = winter_albedo_mask(in_window);
   for pv = string(fieldnames(native_provenance)).'
      native_provenance.(pv) = native_provenance.(pv)(in_window);
   end
   staged_series = staged_series(staged_in, :);

   target = struct('series', series, 'station', site, ...
      'location', location);

   % Donor pool under the role contract; the plan applies the geometry gate.
   donors = assembleDonors(site, met_dir, kwargs);

    proxies = loadStagedProxies(site, location, opts.proxy_catalog, ...
       proxy_window_files);

    % Keep split creation private until the complete artifact transaction
    % publishes. Seed a temporary replay copy from the prior committed split.
    split_manifest = "";
    split_file = "";
    if kwargs.write
       split_file = fullfile(qa_dir, 'splits', site + "-split.json");
       split_manifest = string(tempname) + ".json";
       split_cleanup = onCleanup(@() removeArtifactPath(split_manifest));
       if isfile(split_file)
          copyfile(split_file, split_manifest);
       end
    end

    % Selection experiment and admitted-method plan (albedo plans with the
   % non-interpolating method set; swu follows albedo*swd downstream).
   plan = icemodel.forcing.reconstruct.stationMethodPlan(target, donors, ...
       proxies, seed=opts.seed, channels=opts.plan_channels, ...
       core_channels=opts.core_channels, ...
      n_gaps=opts.plan_n_gaps, knot_candidates=opts.knot_candidates, ...
      max_donors=opts.max_donors, ...
       max_donor_distance_km=opts.max_donor_distance_km, ...
       max_donor_elev_diff_m=opts.max_donor_elev_diff_m, ...
       selection_fraction=opts.selection_fraction, ...
       min_overlap_hours=opts.min_overlap_hours, ...
       max_lag_hours=opts.max_lag_hours, ...
       min_lag_gain=opts.min_lag_gain, ...
       max_extrapolation_fraction=opts.max_extrapolation_fraction, ...
       rmse_improvement=opts.rmse_improvement, ...
       min_variability_ratio=opts.min_variability_ratio, ...
       max_variability_ratio=opts.max_variability_ratio, ...
       min_coverage=opts.min_coverage, lapse_rate=opts.lapse_rate, ...
       elevation_threshold_m=opts.elevation_threshold_m, ...
       tair_for_pressure=opts.tair_for_pressure, ...
       min_season_samples=opts.min_season_samples, ...
       climatology_window_days=opts.climatology_window_days, ...
       climatology_min_support=opts.climatology_min_support, ...
       synthetic_context_hours=opts.synthetic_context_hours, ...
       cap_hours=opts.cap_hours, ...
       cap_hours_by_channel=opts.cap_hours_by_channel, ...
       jump_factor=opts.jump_factor, ...
       toa_dark_wm2=opts.toa_dark_wm2, ...
       split_manifest=split_manifest);
    % Persist central options plus resolved runtime donor controls so later
    % reports describe this run, not live defaults or a different donor pool.
    producer_options = opts;
    producer_options.donor_sites = kwargs.donor_sites;
    producer_options.use_ktransect = kwargs.use_ktransect;
    producer_options.use_gcnet = kwargs.use_gcnet;
    plan.reconstruction_options = producer_options;
   if kwargs.write
      % Blocked triage plans bypass validationSplit, while ordinary plans may
      % have loaded an existing replay. In both cases stage the exact split
      % carried by this plan for the later artifact transaction.
      fid = fopen(split_manifest, 'w');
      cleaner = onCleanup(@() fclose(fid));
      fprintf(fid, '%s', jsonencode(plan.split));
      clear cleaner
   end

   % Compose admitted methods; reconstructSeries stamps provenance and
   % leaves unresolvable gaps missing.
   channel_methods = plannedMethods(plan);
   composed = icemodel.forcing.reconstruct.reconstructSeries(series, ...
      channel_methods, latitude=location.lat_wgs84, ...
      longitude=location.lon_wgs84, ...
      interp_channels=opts.interp_channels, cap_hours=opts.cap_hours, ...
      cap_hours_by_channel=opts.cap_hours_by_channel, ...
      jump_factor=opts.jump_factor, blend_hours=opts.blend_hours, ...
       toa_dark_wm2=opts.toa_dark_wm2, ...
       max_validation_duration_factor= ...
       opts.max_validation_duration_factor, ...
       native_provenance=native_provenance);
   filled = composed.series;
   provenance = composed.provenance;
   audit = composed.audit;

    codes = icemodel.forcing.reconstruct.provenanceCodes();
    % Instrument geometry has no donor/proxy tier. Staged station metadata
    % does not contain a complete maintenance-visit registry, so every gap
    % remains missing in the product; the runtime fallback chain (measured
    % -> interpolated -> nominal, POLICY A3) owns geometry, and gaps never
    % grade readiness.
    [filled, provenance, audit] = fillBoomHeight( ...
       filled, provenance, audit, codes);

    % Staged MODIS daily albedo holds first position in albedo's
    % last-resort order (POLICY A11/B12): residual albedo gaps adopt the
    % satellite observation before any RCM value, attached through the
    % single daily->met-cadence helper, bounds-checked, seam-blended, and
    % stamped with the modis provenance code.
    if opts.use_modis_albedo
       [filled, provenance, audit] = adoptModisAlbedo( ...
          filled, provenance, audit, site, modis_dir, series, codes, opts);
    end

    % Last-resort proxy adoption: residual missing required-channel
   % samples take aligned proxy values in catalog order, bounds-checked
   % and audited. Adopting whole outage spans from one proxy source
   % keeps thermodynamically coupled channels (tair/rh/lwd and their
   % kin) physically consistent with each other, which per-channel
   % method composition cannot guarantee.
   % Final-tier denial notes stay empty when the tier is disabled; the
   % audit reconciliation then keeps the provisional reasons untouched.
   last_resort_denials = struct();
    if opts.last_resort_proxies
        [filled, provenance, audit, last_resort_denials] = ...
             icemodel.forcing.reconstruct.lastResortProxies(filled, ...
            provenance, audit, proxies, codes, opts, ...
            latitude=location.lat_wgs84, longitude=location.lon_wgs84, ...
             native=series, plan=plan);
     end

   % Re-run the same physically aware bounded interpolation after the
   % source ladder. Long native outages can leave only a few refused seam
   % samples after donor/proxy adoption; those samples had no finite
   % anchors during tier 1 but are now ordinary short interior gaps. The
   % second pass keeps the original per-channel caps, native-only seam
   % scale, SWD clear-sky geometry, provenance code, and exact audit
   % (POLICY B1/B3/D-32).
   [filled, provenance, audit] = closeResidualShortGaps( ...
      filled, provenance, audit, series, codes, opts, location);

   % Publication-quality SWD seam pass (POLICY D-32): compare method
   % boundaries with this station's observed steps in the same season and
   % solar-elevation band, then replace only the reconstructed posting
   % beside an empirical outlier. Native samples are immutable and the
   % resulting diagnostics ship with the station plan.
   seam_quality = table();
   if ismember("swd", string(filled.Properties.VariableNames)) ...
         && ismember("swd", string(provenance.Properties.VariableNames))
      seam_cap = opts.cap_hours;
      if isfield(opts.cap_hours_by_channel, 'swd')
         seam_cap = opts.cap_hours_by_channel.swd;
      end
      [filled.swd, provenance.swd, seam_rows, seam_quality] = ...
         icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
         filled.Properties.RowTimes, filled.swd, provenance.swd, ...
         location.lat_wgs84, location.lon_wgs84, ...
         percentile=opts.seam_qa_percentile, ...
         min_reference_steps=opts.seam_qa_min_reference_steps, ...
         max_passes=opts.seam_qa_max_passes, cap_hours=seam_cap);
      if ~isempty(seam_rows)
         audit = [audit; cell2table(vertcat(seam_rows{:}), ...
            'VariableNames', audit.Properties.VariableNames)];
      end
   end

    % Precipitation policy adoption follows the state-channel last resort;
    % phases come from native components, exact complements, or the proxy's
    % own split — never a reconstruction-time partition (POLICY A10/D-18).
     [filled, provenance, audit] = adoptPrecip( ...
        filled, provenance, audit, proxies, codes, opts, series);

   % Winter-albedo fallback (POLICY B13/D-15a): masked winter samples no
   % method, MODIS, or proxy filled take a bounded seasonal BRIDGE — a
   % straight line between the gap-edge observed albedo, floored at the
   % dry-snow value so the dark months never dip below fresh snow. Where
   % no two anchors exist the floor itself applies. Winter sunlight is
   % near zero, so the forcing impact is small; the bridge removes the
   % unphysical step edges the bare constant produced.
   if any(winter_albedo_mask) && ismember("albedo", ...
         string(provenance.Properties.VariableNames))
       x = filled.albedo;
       restore = winter_albedo_mask & ~isfinite(x);
       if any(restore)
          candidate = x;
          % Bridge across every remaining missing sample using the finite
          % neighbors as anchors; applying only at restore keeps other
          % residual gaps untouched by this fallback.
          bridged = fillmissing(x, 'linear', 'EndValues', 'none');
          candidate(restore) = max(bridged(restore), ...
             opts.native_winter_albedo);
          % Single-sided winters with no usable bridge anchor fall back to
          % the floor alone.
          candidate(restore & ~isfinite(candidate)) = ...
             opts.native_winter_albedo;
          [candidate, seam_note] = ...
             icemodel.forcing.reconstruct.blendFallbackSeams( ...
             filled.Properties.RowTimes, series.albedo, x, candidate, ...
             restore, jump_factor=opts.jump_factor, ...
             blend_hours=opts.blend_hours);
          x(restore) = candidate(restore);
          filled.albedo = x;
         code = provenance.albedo;
         code(restore) = codes.constant;
         provenance.albedo = code;
          t = filled.Properties.RowTimes;
           rows = icemodel.forcing.reconstruct.auditSegments(t, restore, ...
              "albedo", "winter_albedo_bridge", ...
              "seasonal bridge floored at the dry-snow value" + ...
              string(seam_note));
          audit = [audit; cell2table(vertcat(rows{:}), ...
             'VariableNames', audit.Properties.VariableNames)];
      end
   end

   % Upward shortwave is a dependency-derived channel. Preserve native swu
   % and fill only its missing samples after swd and albedo finish every
   % admitted, proxy, and constant tier. A staged series with no swu
   % column at all still ships the derived product (POLICY B10): the
   % channel is created empty so the derivation owns every sample instead
   % of being silently skipped and branding swu wholly missing. SWU is not
   % a forcing-plan channel, so an existing native column must explicitly
   % carry its native provenance into this dependent-channel pass.
   names_now = string(filled.Properties.VariableNames);
   provenance_names = string(provenance.Properties.VariableNames);
   if ismember("swu", names_now) && ~ismember("swu", provenance_names)
      if isfield(native_provenance, 'swu')
         if numel(native_provenance.swu) ~= height(provenance)
            error('icemodel:reconstruct:fillPromiceStation:swuProvenance', ...
               'native swu provenance is not aligned to its source axis');
         end
         provenance.swu = native_provenance.swu;
      else
         % Families or fixtures without PROMICE raw-source selection
         % metadata still own their finite staged SWU values. Stamp those
         % exact copies as observed and leave only missing samples for B10.
         swu_code = repmat(codes.missing, height(provenance), 1);
         swu_code(isfinite(filled.swu)) = codes.observed;
         provenance.swu = swu_code;
      end
   elseif ~ismember("swu", names_now) ...
         && all(ismember(["swd", "albedo"], names_now))
      filled.swu = nan(height(filled), 1);
      provenance.swu = repmat(codes.missing, height(provenance), 1);
   end
   if ismember("swu", string(provenance.Properties.VariableNames))
      [filled, provenance, audit] = ...
         icemodel.forcing.reconstruct.deriveUpwardShortwave( ...
         filled, provenance, audit, codes);
   end

   % Later tiers can fill samples that reconstructSeries had already
   % recorded as unfilled. Trim those provisional rows to the final
   % residual mask so the report never overstates missing intervals,
   % and append the last-resort tier's actual denial cause so the
   % shipped reason is never a stale provisional one.
   audit = reconcileUnfilledAudit(audit, filled, last_resort_denials);
   % Reconstruct on the source postings, then restore the original
   % 15-minute support only after every scientific decision is complete.
   % FILLED hourly or half-hourly postings disaggregate mean-preservingly
   % over their quarter-hour slots while OBSERVED postings stay exact held
   % copies (POLICY D-30/A7); provenance codes always repeat.
   if restore_15m
      filled = expandPostingSupport(filled, ...
         staged_series.Properties.RowTimes, provenance, codes);
      provenance = expandPostingSupport(provenance, ...
         staged_series.Properties.RowTimes);
      series = staged_series;
      % B10 on the shipped axis: the derived QC channel must satisfy
      % swu = albedo * swd per DELIVERED sample, so code-12 samples
      % recompute from the expanded operands instead of keeping an
      % independent disaggregation that would contradict the channel's
      % own definition. Where the operands are held copies this is
      % itself an exact mean-preserving expansion of the source-posting
      % derivation (a constant albedo factors out of the posting mean).
      names_out = string(filled.Properties.VariableNames);
      if ismember("swu", string(provenance.Properties.VariableNames)) ...
            && all(ismember(["swu", "swd", "albedo"], names_out))
         rederive = provenance.swu == codes.derived_shortwave;
         filled.swu(rederive) = filled.albedo(rederive) ...
            .* filled.swd(rederive);
      end

      % One coarse source posting can remain unfilled when its two hourly
      % anchors cannot satisfy the low-sun seam envelope. Once restored to
      % the delivered 15-minute axis, the same existing D-32 interpolation
      % physics can prove a smooth scalar-valid bridge over the four (or
      % eight) genuinely missing samples. Re-run only the shared residual
      % closer; observed held copies remain immutable and every adoption
      % retains bounded-interpolation provenance and audit.
      n_audit_before_delivery = height(audit);
      [filled, provenance, audit] = closeResidualShortGaps( ...
         filled, provenance, audit, staged_series, codes, opts, location);
      % The closer can supply a previously missing final SWD or albedo
      % operand. Re-run the existing B10 derivation before audit
      % reconciliation so swu cannot remain missing or stale on the
      % delivered axis.
      if ismember("swu", string(provenance.Properties.VariableNames))
         [filled, provenance, audit] = ...
            icemodel.forcing.reconstruct.deriveUpwardShortwave( ...
            filled, provenance, audit, codes);
      end
      % Distinguish this final delivered-axis proof in the durable audit;
      % it uses the same method code but closes a gap the source-axis pass
      % could not certify. This includes any dependent SWU samples the
      % newly finite operand made derivable.
      for r = n_audit_before_delivery + 1:height(audit)
         audit.detail{r} = ['delivered-axis; ' audit.detail{r}];
      end
      % The first reconciliation already attached the last-resort denial
      % evidence on the source axis. Trim those rows once more against the
      % delivered result so a newly closed posting cannot remain mislabeled
      % as unfilled; no second denial vector exists on this expanded axis.
      audit = reconcileUnfilledAudit(audit, filled, struct());
   end

   % A7: register every emitted context after cadence restoration and the
   % delivered-axis residual pass. This gives every final audit row one
   % persisted plan record without maintaining a second method-name list.
   plan.audit_contexts = auditContextRecords(plan, audit);

   % Publication boundary: incomplete phase data may remain missing, but
   % every finite component must be nonnegative and every complete split
   % must satisfy the exact A10 mass balance on the delivered axis.
   assertPrecipitationProduct(filled);

   % Readiness ledger rows: a forcing verdict and separate native-support
   % confidence advisory for every calendar year of the station record.
   readiness = readinessLedger(site, series, filled, plan, opts, location);

   met_file = "";
   if kwargs.write
      % A6 clips the artifact to the fillable proxy window, so every row in
      % a published product must satisfy both forcing consumers. Exploratory
      % write=false calls still return the diagnostic ledger for repair.
      ready = string(readiness.verdict_icemodel) == "ready" ...
         & string(readiness.verdict_snowmodel) == "ready";
      if ~all(ready)
         retirePublishedArtifacts(site, out_dir, qa_dir, family);
         blocked = readiness(~ready, :);
         details = strings(height(blocked), 1);
         for k = 1:height(blocked)
            details(k) = sprintf( ...
               '%d [IceModel: %s; snow model: %s]', blocked.year(k), ...
               string(blocked.reason_icemodel(k)), ...
               string(blocked.reason_snowmodel(k)));
         end
         error( ...
            'icemodel:reconstruct:fillPromiceStation:notForcingReady', ...
            ['refusing to publish %s because station-year readiness ' ...
            'failed: %s'], site, strjoin(details, "; "));
      end
      met_file = writeArtifacts(site, filled, provenance, audit, ...
          plan, readiness, seam_quality, flat_run_findings, out_dir, qa_dir, ...
          codes, native_file, ...
          acceptance_window, proxy_window_files, split_manifest, split_file, ...
          unique(string({donors.station}), 'stable'), family);
      clear split_cleanup
   end
   result = struct('site', site, 'plan', plan, 'filled', filled, ...
      'provenance', provenance, 'audit', audit, ...
      'readiness', readiness, 'seam_quality', seam_quality, ...
      'flat_run_findings', flat_run_findings, ...
      'met_file', met_file);
end

function assertPrecipitationProduct(filled)
   %ASSERTPRECIPITATIONPRODUCT Refuse an invalid delivered phase split.
   names = string(filled.Properties.VariableNames);
   precip = icemodel.forcing.helpers.precipitationVariables();
   if ~all(ismember(precip, names))
      return
   end
   invalid = ~icemodel.forcing.helpers.precipitationValidity( ...
      filled.ppt, filled.rainf, filled.snowf);
   if any(invalid)
      error('icemodel:reconstruct:fillPromiceStation:invalidPrecipitation', ...
         ['finite delivered precipitation must be nonnegative, each phase ' ...
         'must not exceed ppt, and complete splits must sum to ppt: ' ...
         '%d samples'], nnz(invalid));
   end
end

function pathname = defaultPath(value, repo, relative)
   %DEFAULTPATH Resolve one directory default under the repo root.
   pathname = value;
   if pathname == ""
      pathname = string(fullfile(repo, relative));
   end
end

function [series, winter_mask, restore_15m, native_provenance] = ...
      reconstructionAxis(staged, staged_winter_mask, staged_provenance)
   %RECONSTRUCTIONAXIS Recover source postings from guarded 15-minute support.
   series = staged;
   winter_mask = staged_winter_mask;
   restore_15m = false;
   native_provenance = staged_provenance;
   metadata = staged.Properties.UserData;
   is_guarded_source = false;
   cadence_seconds = NaN;
   if isstruct(metadata) && isscalar(metadata) ...
         && isfield(metadata, 'met_resample_policy') ...
         && isfield(metadata, 'met_resample_source_cadence_seconds')
      policy = string(metadata.met_resample_policy);
      cadence_seconds = metadata.met_resample_source_cadence_seconds;
      % Guarded artifacts stamp their source cadence. Reconstruction
      % accepts hourly sources (PROMICE, IMAU) and half-hourly sources
      % (K-transect) with cadence_seconds / 900 rows per posting (bead
      % icemodel-g1n.49 family generalization).
      is_guarded_source = isscalar(policy) ...
         && policy == "interval_start_zero_order_hold" ...
         && isnumeric(cadence_seconds) && isscalar(cadence_seconds) ...
         && isfinite(cadence_seconds) ...
         && ismember(cadence_seconds, [3600, 1800]);
   end
   if ~is_guarded_source
      times = staged.Properties.RowTimes;
      is_native_hourly = height(staged) >= 2 ...
         && all(diff(times) == hours(1));
      if is_native_hourly
         return
      end
      error('icemodel:reconstruct:fillPromiceStation:unverifiedNativeCadence', ...
         ['native reconstruction requires an hourly axis or a guarded ' ...
         '15-minute artifact stamped from hourly or half-hourly support']);
   end

   times = staged.Properties.RowTimes;
   if height(staged) == 0 || any(diff(times) ~= minutes(15))
      error('icemodel:reconstruct:fillPromiceStation:invalidGuardedAxis', ...
         'guarded source met must have a complete 15-minute axis');
   end
   source_row = mod(seconds(times - times(1)), cadence_seconds) == 0;
   if nnz(source_row) * (cadence_seconds / 900) ~= height(staged)
      error('icemodel:reconstruct:fillPromiceStation:invalidGuardedSupport', ...
         ['guarded source support must contain cadence_seconds / 900 ' ...
         'rows per posting']);
   end
   series = staged(source_row, :);
   winter_mask = staged_winter_mask(source_row);
   names = string(fieldnames(staged_provenance));
   for k = 1:numel(names)
      native_provenance.(names(k)) = staged_provenance.(names(k))(source_row);
   end
   restore_15m = true;
end

function expanded = expandPostingSupport(postings, target_times, ...
      provenance, codes)
   %EXPANDPOSTINGSUPPORT Restore source postings onto their 15-minute axis.
   % Every posting repeats over its quarter-hour slots by default — the
   % exact held-copy behavior provenance codes and observed samples keep.
   % When the posting provenance table and code registry are supplied, the
   % numeric channels' FILLED postings (any code other than observed and
   % the raw/clamped PROMICE pyranometer copies, which are honest copies
   % of native sensor data — A7) instead receive the D-30 mean-preserving
   % disaggregation so the delivered product is smooth while each posting's
   % mean still equals its source value exactly.
   [loc, posting_seconds] = postingLocations( ...
      postings.Properties.RowTimes, target_times);
   n_slots = round(posting_seconds / 900);
   if ~ismember(posting_seconds, [1800, 3600]) ...
         || n_slots * 900 ~= posting_seconds
      error('icemodel:reconstruct:fillPromiceStation:invalidExpansionCadence', ...
         'posting expansion requires 1800- or 3600-second source cadence');
   end
   expanded = postings(loc, :);
   expanded.Properties.RowTimes = target_times;
   if nargin < 4
      return
   end

   % Honest-copy codes: native observations and the raw/clamped PROMICE
   % shortwave fallbacks replay actual sensor readings, so their support
   % stays a zero-order hold (POLICY D-30/A7). Everything else is an
   % estimate the product may disaggregate.
   honest = [codes.observed, codes.raw_shortwave, ...
      codes.clamped_shortwave, codes.darkness];
   % Quarter-hour slot of each product row within its source posting.
   slot = round(seconds(target_times ...
      - postings.Properties.RowTimes(loc)) / 900) + 1;
   for name = string(postings.Properties.VariableNames)
      if ~ismember(name, string(provenance.Properties.VariableNames))
         continue
      end
      values = postings.(name);
      if ~isnumeric(values) || ~iscolumn(values)
         continue
      end
      % Disaggregate only finite filled postings; missing postings hold
      % their NaN copies and honest codes hold their sensor copies.
      posting_filled = ~ismember(provenance.(name), honest) ...
         & isfinite(values);
      if ~any(posting_filled)
         continue
      end
      bounds = disaggregationBounds(name);
      if name == "swd"
         % D-44 is an interval-twilight estimate with a hard 50 W m-2
         % delivered ceiling. Give only its dedicated provenance rows the
         % tighter vector bound; daytime SWD fills keep the scalar registry.
         twilight = provenance.swd == codes.twilight_climatology;
         if any(twilight)
            bands = ...
               icemodel.forcing.reconstruct.solarElevationBands();
            lower = bounds(1) + zeros(height(postings), 1);
            upper = bounds(2) + zeros(height(postings), 1);
            upper(twilight) = bands.twilight_ceiling_wm2;
            bounds = [lower, upper];
         end
      end
      if name == "ppt" && all(ismember(["rainf", "snowf"], ...
            string(postings.Properties.VariableNames)))
         % A native phase is an immovable lower bound on the reconstructed
         % total. This lets the total stay smooth without making its exact
         % complement negative; when both phases are native, their sum pins
         % the total to the honest held copy.
         rain_honest = ismember(provenance.rainf, honest) ...
            & isfinite(postings.rainf);
         snow_honest = ismember(provenance.snowf, honest) ...
            & isfinite(postings.snowf);
         lower = bounds(1) + zeros(height(postings), 1);
         lower(rain_honest) = max(lower(rain_honest), ...
            postings.rainf(rain_honest));
         lower(snow_honest) = max(lower(snow_honest), ...
            postings.snowf(snow_honest));
         both_honest = rain_honest & snow_honest;
         lower(both_honest) = postings.rainf(both_honest) ...
            + postings.snowf(both_honest);
         bounds = [lower, bounds(2) + zeros(height(postings), 1)];
      end
      quarters = disaggregatePostings(values, posting_filled, bounds, n_slots);
      column = expanded.(name);
      % Replace only rows whose SOURCE posting is filled: disaggregation
      % never crosses into an adjacent observed posting's support.
      replace = posting_filled(loc);
      flat = sub2ind(size(quarters), loc(replace), slot(replace));
      column(replace) = quarters(flat);
      expanded.(name) = column;
   end

   % Precipitation phases must satisfy ppt = rainf + snowf per DELIVERED
   % sample (the A10/D-18 phase identity the snowmodel verdict grades).
   % Independent per-channel disaggregation cannot guarantee that when a
   % neighbor posting is missing for one channel but finite for another
   % (the curve degenerates flat on different sides), so FILLED phase
   % samples instead carry the posting's phase FRACTION onto the expanded
   % total. The fraction is constant over the posting, so the phase mean
   % stays exact and the per-sample identity holds wherever both phases
   % are filled.
   precip = icemodel.forcing.helpers.precipitationVariables();
   if all(ismember(precip, string(postings.Properties.VariableNames))) ...
         && all(ismember(precip, ...
         string(provenance.Properties.VariableNames)))
      ppt_posting = postings.ppt;
      for name = setdiff(precip, "ppt", 'stable')
         phase_posting = postings.(name);
         rescale = ~ismember(provenance.(name), honest) ...
            & isfinite(phase_posting) & isfinite(ppt_posting);
         if ~any(rescale)
            continue
         end
         other = precip(precip ~= name & precip ~= "ppt");
         other_honest = ismember(provenance.(other), honest) ...
            & isfinite(postings.(other));
         frac = zeros(size(ppt_posting));
         positive = ppt_posting > 0;
         frac(positive) = phase_posting(positive) ...
            ./ ppt_posting(positive);
         column = expanded.(name);
         % When the opposite phase is native, preserve that held
         % observation and derive this filled phase as its exact
         % complement. Otherwise both filled phases retain the posting
         % source fractions.
         complement = rescale(loc) & other_honest(loc);
         column(complement) = expanded.ppt(complement) ...
            - expanded.(other)(complement);
         fraction = rescale(loc) & ~other_honest(loc);
         column(fraction) = frac(loc(fraction)) .* expanded.ppt(fraction);
         expanded.(name) = column;
      end
   end
end

function [loc, posting_seconds] = postingLocations(posting_times, target_times)
   %POSTINGLOCATIONS Map delivered rows to their owning source postings.
   posting_seconds = mode(seconds(diff(posting_times)));
   offset_seconds = seconds(target_times - posting_times(1));
   snapped = posting_times(1) ...
      + seconds(floor(offset_seconds / posting_seconds) * posting_seconds);
   [present, loc] = ismember(snapped, posting_times);
   if ~all(present)
      error('icemodel:reconstruct:fillPromiceStation:outputAxisMismatch', ...
         'posting reconstruction does not cover the staged output axis');
   end
end

function quarters = disaggregatePostings(values, disaggregate, bounds, n_slots)
   %DISAGGREGATEPOSTINGS Mean-preserving quarter-hour split of postings.
   % Standard conservative disaggregation (POLICY D-30): within each
   % selected posting the quarter-hour samples lie on straight lines from
   % the posting value toward its neighbors (piecewise linear between
   % posting midpoints), then receive one constant mean-restoring shift.
   %
   % Mean-preservation proof: sample centers are evenly spaced across the
   % posting; the first half approaches v from the left neighbor and the
   % second approaches the right neighbor from v. Subtracting each row's
   % base-sample mean offset makes its delivered mean exactly v for either
   % two or four slots. Bound redistribution is zero-sum; rows it cannot
   % settle fall back to held copies of v.
   %
   % Inputs: values (K x 1 postings), disaggregate (K x 1 logical),
   % bounds (1 x 2 or K x 2 inclusive clamp limits), n_slots (2 or 4).
   % Returns K x n_slots samples, held except at disaggregated postings.
   quarters = repmat(values, 1, n_slots);
   rows = find(disaggregate);
   if isempty(rows)
      return
   end
   v = values(rows);
   % Neighbor postings anchor the within-posting slope; a record edge or a
   % missing neighbor degenerates that side to a flat segment so NaN can
   % never leak into the curve.
   left = nan(numel(rows), 1);
   ok = rows > 1;
   left(ok) = values(rows(ok) - 1);
   right = nan(numel(rows), 1);
   ok = rows < numel(values);
   right(ok) = values(rows(ok) + 1);
   left(~isfinite(left)) = v(~isfinite(left));
   right(~isfinite(right)) = v(~isfinite(right));
   % Base curve plus its exact per-posting mean-restoring shift.
   centers = ((1:n_slots) - 0.5) / n_slots;
   left_weight = max(0, 0.5 - centers);
   right_weight = max(0, centers - 0.5);
   S = v + (left - v) .* left_weight + (right - v) .* right_weight;
   S = S - (mean(S, 2) - v);
   if size(bounds, 1) == 1
      lower = bounds(1);
      upper = bounds(2);
   else
      lower = bounds(rows, 1);
      upper = bounds(rows, 2);
   end
   % Clamp into the channel bounds, then restore each posting's mean by
   % spreading the residual over samples still free to move; the loop
   % converges in one pass unless both bounds bind.
   for iteration = 1:8
      S = min(max(S, lower), upper);
      err = v - mean(S, 2);
      moving = abs(err) > eps(max(1, abs(v)));
      if ~any(moving)
         break
      end
      free = (S < upper & err > 0) | (S > lower & err < 0);
      n_free = sum(free, 2);
      adjust = n_slots * err ./ max(n_free, 1);
      adjust(~moving | n_free == 0) = 0;
      S = S + free .* adjust;
   end
   % Postings the redistribution cannot settle (every sample pinned at a
   % bound, or a posting itself outside bounds) fall back to exact held
   % copies — trivially mean-preserving and never worse than the
   % pre-D-30 repetition behavior.
   bad = ~(abs(v - mean(S, 2)) <= 1e-9 * max(1, abs(v))) ...
      | any(~isfinite(S), 2);
   S(bad, :) = repmat(v(bad), 1, n_slots);
   quarters(rows, :) = S;
end

function bounds = disaggregationBounds(channel)
   %DISAGGREGATIONBOUNDS Clamp limits for one disaggregated channel.
   % The scalar registry is the single source (A15); precipitation
   % components share the total's entry because they are the same
   % nonnegative accumulation quantity and hold no registry row of
   % their own.
   if ismember(channel, icemodel.forcing.helpers.precipitationVariables())
      channel = "ppt";
   end
   bounds = icemodel.forcing.reconstruct.physicalBounds(channel);
end

function [filled, provenance, audit] = closeResidualShortGaps( ...
      filled, provenance, audit, native, codes, opts, location)
   %CLOSERESIDUALSHORTGAPS Bridge short slivers exposed by later tiers.
   % The public fillShortGaps implementation remains the single source of
   % interpolation physics. This orchestration pass only supplies the
   % completed series as candidate anchors while freezing seam scales to
   % the untouched native record.
   names = string(filled.Properties.VariableNames);
   native_names = string(native.Properties.VariableNames);
   channels = intersect(opts.interp_channels, names, 'stable');
   times = filled.Properties.RowTimes;
   for channel = channels
      x = filled.(channel);
      if all(isfinite(x))
         continue
      end

      % Synthetic fills must not inflate the seasonal seam tolerance.
      native_x = nan(size(x));
      if ismember(channel, native_names)
         native_x = native.(channel);
      end
      scale = icemodel.forcing.reconstruct.stepScale(times, native_x);
      channel_cap = opts.cap_hours;
      if isfield(opts.cap_hours_by_channel, channel)
         channel_cap = opts.cap_hours_by_channel.(channel);
      end
      [x, interp_adopted, rows] = ...
         icemodel.forcing.reconstruct.fillShortGaps(times, x, channel, ...
         cap_hours=channel_cap, latitude=location.lat_wgs84, ...
         longitude=location.lon_wgs84, jump_factor=opts.jump_factor, ...
         blend_hours=opts.blend_hours, ...
         toa_dark_wm2=opts.toa_dark_wm2, ...
         allow_swd_flux_fallback=channel == "swd", step_scale=scale);
      twilight_adopted = false(size(x));
      twilight_rows = cell(0, 1);
      if channel == "swd" && any(~isfinite(x))
         [x, twilight_adopted, twilight_rows] = ...
            icemodel.forcing.reconstruct.fillTwilightClimatology( ...
            times, x, native_x, location.lat_wgs84, ...
            location.lon_wgs84);
      end
      if ~any(interp_adopted | twilight_adopted)
         continue
      end

      filled.(channel) = x;
      code = provenance.(channel);
      code(interp_adopted) = codes.bounded_interp;
      code(twilight_adopted) = codes.twilight_climatology;
      provenance.(channel) = code;
      % Method identities stay exact; the detail distinguishes this
      % post-final application from the initial composition tiers.
      rows = [rows; twilight_rows]; %#ok<AGROW>
      for r = 1:numel(rows)
         rows{r}{6} = ['post-final residual; ' rows{r}{6}];
      end
      audit = [audit; cell2table(vertcat(rows{:}), ...
         'VariableNames', audit.Properties.VariableNames)]; %#ok<AGROW>
   end
end

function audit = reconcileUnfilledAudit(audit, filled, denials)
   %RECONCILEUNFILLEDAUDIT Restrict provisional rows to final missing samples.
   % Later tiers can both FILL provisionally-unfilled samples and DENY
   % them for a new reason; each residual row therefore appends the
   % final-tier denial notes lastResortProxies recorded so the shipped
   % detail states the actual cause instead of a stale provisional one.
   is_unfilled = string(audit.method) == "unfilled";
   if ~any(is_unfilled)
      return
   end

   provisional = audit(is_unfilled, :);
   audit = audit(~is_unfilled, :);
   times = filled.Properties.RowTimes;
   rows = {};
   for r = 1:height(provisional)
      channel = string(provisional.channel{r});
      if ~ismember(channel, string(filled.Properties.VariableNames))
         continue
      end
      % duration_hours is the cadence-independent, end-exclusive span.
      % Using it instead of the final timestamp lets a source-hour audit row
      % reconcile correctly after restoration to four delivered quarters.
      span_end = provisional.start_time(r) ...
         + hours(provisional.duration_hours(r));
      residual = times >= provisional.start_time(r) ...
         & times < span_end ...
         & ~isfinite(filled.(channel));
      detail = string(provisional.detail{r});
      % Append the union of last-tier denial notes over the residual
      % samples (usually one note per outage; the union keeps mixed
      % spans honest without per-sample row inflation).
      if isfield(denials, channel) && any(residual)
         notes = unique(denials.(channel)(residual));
         notes = notes(notes ~= "");
         if ~isempty(notes)
            detail = detail + "; final tier: " + strjoin(notes, "; ");
         end
      end
      segment_rows = icemodel.forcing.reconstruct.auditSegments( ...
         times, residual, channel, "unfilled", detail, ...
         context_id=string(provisional.context_id{r}));
      rows = [rows; segment_rows]; %#ok<AGROW>
   end
   if ~isempty(rows)
      audit = [audit; cell2table(vertcat(rows{:}), ...
         'VariableNames', audit.Properties.VariableNames)];
   end
end

function [series, location, winter_mask, filename, native_provenance, ...
      flat_run_findings] = loadStationMet(met_dir, site, opts, family)
   %LOADSTATIONMET Load one staged native met timetable and its point.
   % The PROMICE native builder historically filled every albedo gap and
   % selected
   % missing deep-night shortwave as physical zero. Where the source NetCDF
   % is recorded, replay the builder selection to distinguish observations
   % from both classes of derived values; only observations may fit or seed
   % reconstruction. Builder-derived SWD zeros re-enter the darkness rule,
   % while builder-derived SWU zeros re-enter the final albedo-times-SWD rule.
   % The returned winter mask lets the driver restore the approved dry-snow
   % constant honestly where no better method fills. That replay machinery
   % is promice-only; other families load through the same grammar with the
   % family token in place of promice (bead icemodel-g1n.49).
   hits = dir(fullfile(met_dir, ...
      sprintf('met_%s_%s_*_15m.mat', site, family)));
   % Match the native product grammar exactly; <family>_filled artifacts may
   % coexist in a caller-selected directory but can never become native input.
   names = string({hits.name});
   pattern = sprintf( ...
      '^met_%s_%s_[0-9]{8}_[0-9]{8}_15m\\.mat$', site, family);
   keep = arrayfun(@(name) ~isempty(regexp(char(name), pattern, 'once')), ...
      names);
   hits = hits(keep);
   if isempty(hits)
      error('icemodel:reconstruct:fillPromiceStation:missingNativeMet', ...
         'no staged native met for %s under %s', site, met_dir);
   end
   [series, filename] = ...
      icemodel.forcing.reconstruct.loadWidestTimetable(hits);
   if isempty(series)
      error('icemodel:reconstruct:fillPromiceStation:invalidNativeMet', ...
         'staged native met files for %s contain no timetable', site);
   end
   % Source names die at the read boundary (POLICY A16/D-24): staged
   % artifacts predating the rename carry the pypromice name usr, so accept
   % either spelling here and carry only the canonical swu forward. Both
   % names in one artifact is ambiguous corruption, never a merge.
   staged_names = string(series.Properties.VariableNames);
   if all(ismember(["usr", "swu"], staged_names))
      error(['icemodel:reconstruct:fillPromiceStation:' ...
         'ambiguousShortwaveChannel'], ...
         'staged native met for %s carries both usr and swu', site);
   end
   if ismember("usr", staged_names)
      series = renamevars(series, "usr", "swu");
   end
   if ismember("usr_provenance", staged_names)
      series = renamevars(series, "usr_provenance", "swu_provenance");
   end
    ud = series.Properties.UserData;
    % Identity may ride ud.site (the PROMICE staging grammar) or ud.station
    % (the IMAU staging grammar); either token must normalize to the
    % requested site (bead icemodel-g1n.49).
    has_site = isstruct(ud) && isfield(ud, 'site') ...
       && (ischar(ud.site) || (isstring(ud.site) && isscalar(ud.site)));
    has_station = isstruct(ud) && isfield(ud, 'station') ...
       && (ischar(ud.station) ...
       || (isstring(ud.station) && isscalar(ud.station)));
    identity = "";
    if has_site
       identity = string(ud.site);
    elseif has_station
       identity = string(ud.station);
    end
    if (~has_site && ~has_station) || ...
          icemodel.forcing.helpers.normalizedFileToken(identity) ~= ...
          icemodel.forcing.helpers.normalizedFileToken(site)
      found_site = "<missing>";
      if has_site || has_station
         found_site = identity;
      end
      error('icemodel:reconstruct:fillPromiceStation:nativeIdentityMismatch', ...
         ['native met metadata must identify the requested station %s; ' ...
         'found %s'], site, found_site);
   end
    if isstruct(ud) && isfield(ud, 'gapfill_product')
       error('icemodel:reconstruct:fillPromiceStation:nativeIdentityMismatch', ...
          'native met input for %s identifies a reconstructed product', site);
    end
    % PROMICE met metadata names the point lat/lon; accept the top-level
    % wgs84 spelling other families use, and the site_location struct the
    % IMAU staging writes — its top-level metadata carries lat_wgs84 and
    % lon_wgs84 but no elev_m, so the complete point lives only in the
    % struct form (bead icemodel-g1n.49) — so donors and non-promice
    % targets load through the same path.
   if isfield(ud, 'lat')
      location = struct('lat_wgs84', ud.lat, 'lon_wgs84', ud.lon, ...
         'elev_m', ud.elev);
   elseif all(isfield(ud, ["lat_wgs84", "lon_wgs84", "elev_m"]))
      location = struct('lat_wgs84', ud.lat_wgs84, ...
         'lon_wgs84', ud.lon_wgs84, 'elev_m', ud.elev_m);
   elseif isfield(ud, 'site_location') && isstruct(ud.site_location) ...
         && isscalar(ud.site_location) && all(isfield(ud.site_location, ...
         ["lat_wgs84", "lon_wgs84", "elev_m"]))
      location = struct('lat_wgs84', ud.site_location.lat_wgs84, ...
         'lon_wgs84', ud.site_location.lon_wgs84, ...
         'elev_m', ud.site_location.elev_m);
   else
      % No usable location form: keep the legacy failure shape (the same
      % top-level field access the pre-family loader made).
      location = struct('lat_wgs84', ud.lat_wgs84, ...
         'lon_wgs84', ud.lon_wgs84, 'elev_m', ud.elev_m);
   end

   % Opt-in builder estimates are derived candidates, never native
   % observations that may train or validate reconstruction. The flag and
   % policy must agree so a stale flag cannot relabel an empirical estimate.
   has_lwd_flag = isstruct(ud) && isfield(ud, 'lwd_estimated');
   lwd_estimated = false;
   if has_lwd_flag
      flag = ud.lwd_estimated;
      valid_flag = isscalar(flag) && (islogical(flag) ...
         || (isnumeric(flag) && isfinite(flag) && ismember(flag, [0, 1])));
      if ~valid_flag
         error(['icemodel:reconstruct:fillPromiceStation:' ...
            'inconsistentLwdProvenance'], ...
            'native met input for %s has an invalid lwd_estimated flag', site);
      end
      lwd_estimated = logical(flag);
   end
   has_lwd_policy = isstruct(ud) && isfield(ud, 'lwd_policy');
   policy_estimated = false;
   if has_lwd_policy
      lwd_policy = string(ud.lwd_policy);
      if ~isscalar(lwd_policy)
         error(['icemodel:reconstruct:fillPromiceStation:' ...
            'inconsistentLwdProvenance'], ...
            'native met input for %s has a non-scalar lwd policy', site);
      end
      lower_policy = lower(lwd_policy);
      policy_estimated = contains(lower_policy, "estimated") ...
         || contains(lower_policy, "empirical") ...
         || contains(lower_policy, "fill_lwd=true");
   end
   if (has_lwd_flag || policy_estimated) ...
         && lwd_estimated ~= policy_estimated
      error(['icemodel:reconstruct:fillPromiceStation:' ...
         'inconsistentLwdProvenance'], ...
         'native met input for %s has inconsistent estimated-LWD metadata', ...
         site);
   end
   if ismember("lwd", string(series.Properties.VariableNames)) ...
         && lwd_estimated
      series.lwd(:) = NaN;
   end

    % PROMICE-only builder machinery (bead icemodel-g1n.49): the raw-source
    % replay, the raw-fallback shortwave QC, and the winter-albedo stamp
    % recovery all reconstruct decisions the PROMICE builder is known to
    % have made (legacy albedo fills, source-selected shortwave). Other
    % families' staging pipelines carry no such builder state to replay.
    if family == "promice"
       [observed, raw_albedo, resolved, swd_darkness, native_provenance] = ...
          sourcePromiceMasks( ...
          series, site, ...
          icemodel.forcing.reconstruct.selectedDataRoot(met_dir), location);
       if any(swd_darkness)
          series.swd(swd_darkness) = NaN;
       end
       codes = icemodel.forcing.reconstruct.provenanceCodes();
       if isfield(native_provenance, 'swu')
          swu_darkness = native_provenance.swu == codes.darkness;
          series.swu(swu_darkness) = NaN;
       end

       % A15 hard limits apply at EVERY tier, including the builder's
       % raw-pyranometer fallback: the fallback exists to prefer raw over
       % nothing (A7), never to ship impossible incident energy, yet the
       % GC-Net-legacy cohort carries ~41k code-13 samples above the swd
       % ceiling with an evening-shifted diurnal signature (e.g. 945-962
       % W/m2 at ~1.6x TOA). Gate the raw-fallback support through the same
       % candidate validity rule the fill tiers obey; failing samples turn
       % MISSING so later tiers may fill them honestly. The staged native
       % artifact is untouched — only this run's working series drops them.
       % swd gates first so the swu pairing test sees the post-gate swd: a
       % raw swd proven impossible cannot vouch for its paired raw swu.
       for name = ["swd", "swu"]
          if ~isfield(native_provenance, name) ...
                || ~ismember(name, string(series.Properties.VariableNames))
             continue
          end
          is_raw = native_provenance.(name) == codes.raw_shortwave;
          if ~any(is_raw)
             continue
          end
           if name == "swd"
              valid = icemodel.forcing.reconstruct.physicalValidity("swd", ...
                 series.swd, series.Properties.RowTimes, ...
                 latitude=location.lat_wgs84, longitude=location.lon_wgs84, ...
                 interval=median(diff(series.Properties.RowTimes)));
          else
             valid = icemodel.forcing.reconstruct.physicalValidity("swu", ...
                series.swu, series.Properties.RowTimes, swd=series.swd);
          end
          drop = is_raw & ~valid;
          series.(name)(drop) = NaN;
          native_provenance.(name)(drop) = codes.missing;
       end

       winter_mask = false(height(series), 1);
       if ismember("albedo", string(series.Properties.VariableNames))
          m = month(series.Properties.RowTimes);
          stamped = ismember(m, opts.native_winter_months) ...
             & series.albedo == opts.native_winter_albedo;
          if resolved
             source_match = observed & series.albedo == raw_albedo;
             legacy_fill = isfinite(series.albedo) ...
                & (~observed | (stamped & ~source_match));
             winter_mask = stamped & legacy_fill;
             series.albedo(legacy_fill) = NaN;
          elseif any(isfinite(series.albedo))
             error( ...
                'icemodel:reconstruct:fillPromiceStation:missingAlbedoProvenance', ...
                ['staged %s has finite albedo without verifiable raw-source ' ...
                'provenance'], site)
          end
       end
    else
       % Non-promice albedo-provenance policy (bead icemodel-g1n.49): the
       % promice fail-closed missingAlbedoProvenance rule exists because
       % the PROMICE builder is known to have injected legacy winter fills
       % that must be proven absent by raw-source replay. IMAU and
       % K-transect staging has no such fill history, so finite native
       % albedo ships as observed, no raw replay runs, and the winter
       % stamp mask stays all-false.
       native_provenance = struct();
       winter_mask = false(height(series), 1);
    end
    % Screen PROMICE sensor-burial/rime signatures before this working
    % copy can train a method, set a seam scale, or act as a donor. The
    % staged artifact remains byte-identical; implicated samples become
    % reconstructable missing values with a persisted evidence table.
    flat_run_findings = table();
    if family == "promice"
       [series, native_provenance, flat_run_findings] = ...
          maskFlatRunFindings(series, native_provenance, location);
    end
    verifyNativeMetIdentity(filename, ...
       icemodel.forcing.reconstruct.selectedDataRoot(met_dir), site, family);
 end

function [series, native_provenance, findings] = ...
      maskFlatRunFindings(series, native_provenance, location)
   %MASKFLATRUNFINDINGS Exclude implicated native channels from working data.
   [~, findings] = icemodel.forcing.reconstruct.flatRunScreen( ...
      series, location.lat_wgs84, location.lon_wgs84);
   if isempty(findings)
      return
   end

   % Each finding names only the channels supported by its corroborating
   % evidence. Mask exactly those channels over the reported run.
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   names = string(series.Properties.VariableNames);
   times = series.Properties.RowTimes;
   for row = 1:height(findings)
      in_run = times >= findings.start_time(row) ...
         & times <= findings.end_time(row);
      channels = split(string(findings.channels(row)), ",");
      for channel = reshape(channels, 1, [])
         if ~ismember(channel, names)
            continue
         end
         series.(channel)(in_run) = NaN;
         if isfield(native_provenance, channel)
            native_provenance.(channel)(in_run) = codes.missing;
         end
      end
   end
end

function [observed, raw_albedo, resolved, swd_darkness, native_provenance] = ...
      sourcePromiceMasks(series, site, data_root, location)
   %SOURCEPROMICEMASKS Recover raw support for builder-derived native values.
   observed = false(height(series), 1);
   raw_albedo = nan(height(series), 1);
   resolved = false;
   swd_darkness = false(height(series), 1);
   native_provenance = struct();
   metadata = series.Properties.UserData;
   legacy_albedo = isstruct(metadata) ...
      && isfield(metadata, 'albedo_policy') ...
      && contains(string(metadata.albedo_policy), "fillPromiceAlbedo");
    mask_names = ["swd_raw_fallback", "swd_negative_clamped", ...
       "swd_darkness_fill", "swu_raw_fallback", "swu_negative_clamped", ...
       "swu_darkness_fill"];
    count_fields = ["swd_raw_fallback_count", ...
       "swd_negative_clamped_count", "swd_darkness_zero_filled_count", ...
       "swu_raw_fallback_count", "swu_negative_clamped_count", ...
       "swu_darkness_zero_filled_count"];
   counts = zeros(size(count_fields));
   for k = 1:numel(count_fields)
      if isstruct(metadata) && isfield(metadata, count_fields(k))
         counts(k) = double(metadata.(count_fields(k)));
      end
      if ~isscalar(counts(k)) || ~isfinite(counts(k)) ...
            || counts(k) < 0 || fix(counts(k)) ~= counts(k)
         error(['icemodel:reconstruct:fillPromiceStation:' ...
            'invalidShortwaveProvenance'], ...
            'staged %s has an invalid %s', site, count_fields(k));
      end
   end
   if ~legacy_albedo && ~any(counts)
      return
   end
   source_file = "";
   if isfield(metadata, 'source_file')
      source_file = string(metadata.source_file);
   end
   source_file = resolveSourceFile(source_file, data_root);
   if source_file == "" || ~isfile(source_file)
      if legacy_albedo
         error('icemodel:reconstruct:fillPromiceStation:missingAlbedoSource', ...
            ['staged %s albedo records legacy filling but its raw source ' ...
            'file is unavailable: %s'], site, source_file)
      end
      error('icemodel:reconstruct:fillPromiceStation:missingShortwaveSource', ...
         ['staged %s shortwave records source selection but its raw source ' ...
         'file is unavailable: %s'], site, source_file)
   end
   if ~icemodel.internal.isPathInside(source_file, data_root)
      if legacy_albedo
         error( ...
            'icemodel:reconstruct:fillPromiceStation:albedoSourceOutsideRoot', ...
            ['staged %s raw albedo source is outside the selected data ' ...
            'root %s: %s'], site, data_root, source_file)
      end
      error(['icemodel:reconstruct:fillPromiceStation:' ...
         'shortwaveSourceOutsideRoot'], ...
         ['staged %s raw shortwave source is outside the selected data ' ...
         'root %s: %s'], site, data_root, source_file)
   end
   source_info = dir(source_file);
   has_source_pin = isfield(metadata, 'source_size_bytes') ...
      && isscalar(metadata.source_size_bytes) ...
      && isfinite(metadata.source_size_bytes) ...
      && isfield(metadata, 'source_sha256') ...
      && isscalar(string(metadata.source_sha256)) ...
      && strlength(string(metadata.source_sha256)) == 64;
    % Legacy artifacts predate raw-byte pins. Their staged bytes must still
    % pass producer-manifest verification before loadStationMet returns, while
    % the support replay below binds the current raw file to every relevant
    % derived mask. Newer artifacts carry pins, which remain strict when present.
    if has_source_pin ...
          && (source_info.bytes ~= double(metadata.source_size_bytes) ...
          || icemodel.verification.setup.fileSha256(source_file) ...
          ~= lower(string(metadata.source_sha256)))
       error(['icemodel:reconstruct:fillPromiceStation:' ...
          'rawSourceIdentityMismatch'], ...
         'staged %s raw-source bytes do not match builder metadata: %s', ...
         site, source_file);
   end

   source_dir = string(fileparts(source_file));
   raw = icemodel.forcing.readPromiceAws(site, source_dir=source_dir);
   if legacy_albedo
      if ~ismember("albedo", string(raw.Properties.VariableNames))
         error('icemodel:reconstruct:fillPromiceStation:missingRawAlbedo', ...
            'raw PROMICE source for %s has no albedo channel', site);
      end

      % Reuse the canonical forward-support resampler so a finite hourly
      % observation owns exactly its four quarter-hour samples.
      raw_mask = timetable(raw.Properties.RowTimes, ...
          double(icemodel.forcing.promiceAlbedoSourceValid(raw.albedo)), ...
          raw.albedo, ...
          'VariableNames', {'observed', 'albedo'});
      support = icemodel.forcing.helpers.resampleMetTimestep(raw_mask, "15m");
      [tf, loc] = ismember(series.Properties.RowTimes, ...
         support.Properties.RowTimes);
      observed(tf) = support.observed(loc(tf)) == 1;
      raw_albedo(tf) = support.albedo(loc(tf));
      resolved = true;
   end

   if ~any(counts)
      return
   end
   swd_file_support = NaN;
   swu_file_support = NaN;
   if isfield(metadata, 'swd_source_file_observations_present')
      swd_file_support = double( ...
         metadata.swd_source_file_observations_present);
   end
   if isfield(metadata, 'swu_source_file_observations_present')
      swu_file_support = double( ...
         metadata.swu_source_file_observations_present);
   end
   [~, ~, ~, masks] = icemodel.forcing.helpers.promiceShortwave( ...
      raw, fill_darkness=true, latitude=location.lat_wgs84, ...
      longitude=location.lon_wgs84, ...
      swd_source_file_observations_present=swd_file_support, ...
      swu_source_file_observations_present=swu_file_support);
   in_window = ismember(raw.Properties.RowTimes, ...
      series.Properties.RowTimes);
   for k = 1:numel(mask_names)
      if nnz(masks.(mask_names(k)) & in_window) ~= counts(k)
         error(['icemodel:reconstruct:fillPromiceStation:' ...
            'shortwaveProvenanceMismatch'], ...
            ['staged %s shortwave source-selection counts cannot be ' ...
            'reproduced from its raw source'], site);
      end
   end

    % Carry source-selection provenance over the builder's forward-held
    % quarter-hour support. A clamp overrides raw or darkness when masks
    % overlap. The staged series already carries canonical names (the read
    % boundary renamed any legacy usr), so mask and channel names coincide.
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   for name = ["swd", "swu"]
      if ~ismember(name, string(series.Properties.VariableNames))
         continue
      end
       raw_mask = timetable(raw.Properties.RowTimes, ...
          double(masks.(name + "_raw_fallback")), ...
          double(masks.(name + "_negative_clamped")), ...
          double(masks.(name + "_darkness_fill")), ...
          'VariableNames', ...
          {'raw_fallback', 'negative_clamped', 'darkness_fill'});
      support = icemodel.forcing.helpers.resampleMetTimestep(raw_mask, "15m");
      [tf, loc] = ismember(series.Properties.RowTimes, ...
         support.Properties.RowTimes);
      code = repmat(codes.observed, height(series), 1);
       raw_fallback = false(height(series), 1);
       negative_clamped = false(height(series), 1);
       darkness_fill = false(height(series), 1);
       raw_fallback(tf) = support.raw_fallback(loc(tf)) == 1;
       negative_clamped(tf) = support.negative_clamped(loc(tf)) == 1;
       darkness_fill(tf) = support.darkness_fill(loc(tf)) == 1;
       finite = isfinite(series.(name));
       code(raw_fallback & finite) = codes.raw_shortwave;
       code(darkness_fill & finite) = codes.darkness;
       code(negative_clamped & finite) = codes.clamped_shortwave;
      native_provenance.(name) = code;
   end

   raw_mask = timetable(raw.Properties.RowTimes, ...
      double(masks.swd_darkness_fill), ...
      'VariableNames', {'builder_darkness'});
   support = icemodel.forcing.helpers.resampleMetTimestep(raw_mask, "15m");
   [tf, loc] = ismember(series.Properties.RowTimes, ...
      support.Properties.RowTimes);
   swd_darkness(tf) = support.builder_darkness(loc(tf)) == 1 ...
      & series.swd(tf) == 0;
end

function source_file = resolveSourceFile(source_file, data_root)
   %RESOLVESOURCEFILE Rebind recorded raw provenance after root relocation.
   if source_file == ""
      return
   end
   if isfile(source_file) && icemodel.internal.isPathInside(source_file, data_root)
      return
   end

   % A relative record resolves directly; historical absolute records fall
   % back to one unambiguous same-named source inside the selected data root.
   recorded = java.io.File(char(source_file));
   if ~recorded.isAbsolute()
      candidate = string(fullfile(data_root, source_file));
      if isfile(candidate) && icemodel.internal.isPathInside(candidate, data_root)
         source_file = candidate;
         return
      end
   end
   [~, name, ext] = fileparts(source_file);
   hits = dir(fullfile(data_root, '**', name + ext));
   hits = hits(~[hits.isdir]);
   if isscalar(hits)
      source_file = string(fullfile(hits.folder, hits.name));
   end
end

function verifyNativeMetIdentity(filename, data_root, site, family)
   %VERIFYNATIVEMETIDENTITY Match staged file identity to the producer manifest.
   % The family token selects the producer manifest data/eval/<family> and
   % its colocation.<family> staged-met leg (bead icemodel-g1n.49); the
   % identity rules below are family-neutral.
   manifest_file = fullfile(data_root, 'eval', char(family), 'manifest.json');
   % Family label used by every message here; promice keeps its historical
   % uppercase spelling because upper() reproduces it exactly.
   family_label = upper(family);
   if ~isfile(manifest_file)
      error(['icemodel:reconstruct:fillPromiceStation:' ...
         'missingNativeArtifactIdentity'], ...
         '%s producer manifest is unavailable for staged %s: %s', ...
         family_label, site, manifest_file);
   end

   manifest = jsondecode(fileread(manifest_file));
   if ~isstruct(manifest) || ~isfield(manifest, 'cases') ...
         || ~isstruct(manifest.cases)
      error(['icemodel:reconstruct:fillPromiceStation:' ...
         'missingNativeArtifactIdentity'], ...
         '%s producer manifest has no case records: %s', ...
         family_label, manifest_file);
   end
   cases = manifest.cases;
   case_ids = strings(numel(cases), 1);
   for n = 1:numel(cases)
      if isfield(cases(n), 'case_id')
         case_ids(n) = string(cases(n).case_id);
      elseif isfield(cases(n), 'site_id')
         case_ids(n) = string(cases(n).site_id);
      end
   end
   match_case = icemodel.forcing.helpers.normalizedFileToken(case_ids) ...
      == icemodel.forcing.helpers.normalizedFileToken(site);
   if nnz(match_case) ~= 1
      error(['icemodel:reconstruct:fillPromiceStation:' ...
         'missingNativeArtifactIdentity'], ...
         '%s producer manifest has no unique case for staged %s', ...
         family_label, site);
   end

   entry = cases(match_case);
    valid_leg = isfield(entry, 'colocation') ...
       && isstruct(entry.colocation) ...
       && isfield(entry.colocation, char(family)) ...
       && isstruct(entry.colocation.(char(family)));
    if ~valid_leg
       error(['icemodel:reconstruct:fillPromiceStation:' ...
          'missingNativeArtifactIdentity'], ...
          '%s producer manifest has no staged met leg for %s', ...
          family_label, site);
    end

    [~, stem, ext] = fileparts(filename);
    basename = string(stem) + string(ext);
    leg = entry.colocation.(char(family));
    has_identities = isfield(leg, 'met_file_identities') ...
       && isstruct(leg.met_file_identities) ...
       && ~isempty(leg.met_file_identities);
    if ~has_identities
       if family == "promice"
          error(['icemodel:reconstruct:fillPromiceStation:' ...
             'missingNativeArtifactIdentity'], ...
             ['PROMICE producer manifest has no size+SHA-256 staged met ' ...
             'identity for %s; refresh the manifest before reconstruction'], ...
             site);
       end
       % Legacy manifests pin the canonical staged path but predate per-file
       % hashes. Preserve that identity contract for non-PROMICE families;
       % A1/D-22 requires every PROMICE target and donor to take the
       % byte-strict branch below.
       if ~isfield(leg, 'met_files')
          error(['icemodel:reconstruct:fillPromiceStation:' ...
             'missingNativeArtifactIdentity'], ...
             '%s producer manifest has no staged met identity for %s', ...
             family_label, site);
       end
       declared = replace(string(leg.met_files), "\", "/");
       match_file = declared == basename | endsWith(declared, "/" + basename);
       if nnz(match_file) ~= 1
          error(['icemodel:reconstruct:fillPromiceStation:' ...
             'missingNativeArtifactIdentity'], ...
             '%s producer manifest has no unique staged path for %s', ...
             family_label, filename);
       end
       return
    end

    identities = leg.met_file_identities;
    declared = replace(string({identities.file}), "\", "/");
   match_file = declared == basename | endsWith(declared, "/" + basename);
   if nnz(match_file) ~= 1
      error(['icemodel:reconstruct:fillPromiceStation:' ...
         'missingNativeArtifactIdentity'], ...
         '%s producer manifest has no unique identity for %s', ...
         family_label, filename);
   end

   identity = identities(match_file);
   valid_identity = isfield(identity, 'size_bytes') ...
      && isscalar(identity.size_bytes) && isfinite(identity.size_bytes) ...
      && isfield(identity, 'sha256') ...
      && isscalar(string(identity.sha256)) ...
      && strlength(string(identity.sha256)) == 64;
   info = dir(filename);
   if ~valid_identity || info.bytes ~= double(identity.size_bytes) ...
         || icemodel.verification.setup.fileSha256(filename) ...
         ~= lower(string(identity.sha256))
      error(['icemodel:reconstruct:fillPromiceStation:' ...
         'nativeArtifactIdentityMismatch'], ...
         'staged native met bytes do not match the %s producer manifest: %s', ...
         family_label, filename);
   end
end

function donors = assembleDonors(site, met_dir, kwargs)
   %ASSEMBLEDONORS Build the role-typed donor pool for one target station.
   donor_cells = {};
   % Canonical staged met lives at <data_root>/input/met/<family>; every
   % donor family must stay inside that same selected root. Fixture and
   % caller-selected roots may omit the canonical input/ layer.
    [data_root, met_root] = ...
       icemodel.forcing.reconstruct.selectedDataRoot(met_dir);

   % Station-donor family set (bead icemodel-g1n.49): PROMICE is the only
   % staged station-donor registry today, so every target family — promice
   % itself and non-promice targets such as imau — draws station donors
   % from the promice network. A future donor family joins by extending
   % this set together with the per-neighbor family records the site
   % inventory carries.
   donor_families = "promice";
   % A promice target keeps its caller-selected met_dir (fixture roots may
   % relocate it); a non-promice target reads its promice donors from the
   % canonical sibling staged directory under the same met root.
   if kwargs.family == "promice"
      donor_met_dir = met_dir;
   else
      donor_met_dir = string(fullfile(met_root, 'promice'));
   end

   % Other staged PROMICE stations (nearest-first selection happens in the
   % plan's geometry gate, so loading order does not matter).
   donor_sites = kwargs.donor_sites;
   if isscalar(donor_sites) && donor_sites == "auto"
      % Preselect through the precomputed proximity inventory so only
      % geometry-plausible neighbors are loaded, not the whole network.
      donor_sites = inventoryDonors(data_root, site, kwargs.opts, ...
         donor_families);
      if isempty(donor_sites)
         if kwargs.family == "promice"
            hits = dir(fullfile(met_dir, 'met_*_promice_*_15m.mat'));
            tokens = extractBetween(string({hits.name}), "met_", "_promice");
            donor_sites = setdiff(tokens, site);
         else
            % A target family absent from the donor registry degrades to
            % no station donors rather than loading every staged promice
            % station unscreened (bead icemodel-g1n.49); K-transect and
            % GC-Net donors below still apply under their own toggles.
            warning(['icemodel:reconstruct:fillPromiceStation:' ...
               'noDonorInventory'], ...
               ['no donor inventory entry for %s (family %s); ' ...
               'proceeding with no station donors'], site, kwargs.family);
         end
      end
   end
   promice_cells = cell(numel(donor_sites), 1);
   for k = 1:numel(donor_sites)
      try
          [d_staged, d_location, d_winter_mask, ~, d_provenance, ~] = ...
             loadStationMet(donor_met_dir, donor_sites(k), kwargs.opts, ...
             "promice");
          [d_series, ~, ~, ~] = reconstructionAxis( ...
             d_staged, d_winter_mask, d_provenance);
      catch err
         if ismember(string(err.identifier), [ ...
                "icemodel:reconstruct:fillPromiceStation:missingAlbedoSource"
                "icemodel:reconstruct:fillPromiceStation:missingRawAlbedo"
                 "icemodel:reconstruct:fillPromiceStation:rawSourceIdentityMismatch"
                 "icemodel:reconstruct:fillPromiceStation:missingNativeArtifactIdentity"
                 "icemodel:reconstruct:fillPromiceStation:nativeArtifactIdentityMismatch"
                 "icemodel:reconstruct:fillPromiceStation:unverifiedNativeCadence"])
            rethrow(err)
         end
         continue
      end
      promice_cells{k} = struct('series', d_series, ...
         'station', donor_sites(k), 'family', "promice", ...
         'location', d_location, 'observed_mask', []);
   end
   donor_cells = [donor_cells; promice_cells(~cellfun(@isempty, ...
      promice_cells))];

   % Staged K-transect cases: native 30-minute observations align onto the
   % 15-minute axis by exact timestamps inside the plan.
   if kwargs.use_ktransect
      kt_root = fullfile(data_root, 'eval', 'ktransect');
      kt_manifest = fullfile(kt_root, 'manifest.json');
      if isfile(kt_manifest)
          kt = jsondecode(fileread(kt_manifest));
          kt_cells = cell(numel(kt.cases), 1);
          for k = 1:numel(kt.cases)
             entry = kt.cases(k);
             leg = entry.colocation.ktransect;
             evaluation_file = fullfile(kt_root, ...
                leg.evaluation_file);
             if ~icemodel.internal.isPathInside(evaluation_file, kt_root)
                error(['icemodel:reconstruct:fillPromiceStation:' ...
                   'ktransectPathOutsideRoot'], ...
                   'K-transect donor path escapes the selected root: %s', ...
                   evaluation_file);
             end
              required = ["doi", "bundle_doi", "license", "children"];
              info = dir(evaluation_file);
              has_size = isfield(leg, 'evaluation_size_bytes');
              has_hash = isfield(leg, 'evaluation_sha256');
              has_pin = has_size && has_hash;
              invalid_pin = xor(has_size, has_hash);
              if ~isempty(info) && has_pin
                 invalid_pin = ...
                    info.bytes ~= double(leg.evaluation_size_bytes) ...
                    || icemodel.verification.setup.fileSha256(evaluation_file) ...
                    ~= lower(string(leg.evaluation_sha256));
              end
              if isempty(info) || ~all(isfield(leg, required)) || invalid_pin
                 error(['icemodel:reconstruct:fillPromiceStation:' ...
                    'ktransectArtifactIdentityMismatch'], ...
                    'K-transect donor identity does not match its manifest: %s', ...
                    evaluation_file);
              end
             obs = load(evaluation_file);
             if ~isfield(obs, 'targets') || ~isstruct(obs.targets) ...
                   || ~isfield(obs.targets, 'data') ...
                   || ~istimetable(obs.targets.data) ...
                   || ~isfield(obs.targets, 'metadata')
                error(['icemodel:reconstruct:fillPromiceStation:' ...
                   'ktransectArtifactIdentityMismatch'], ...
                   'K-transect donor artifact has an invalid schema: %s', ...
                   evaluation_file);
             end
             metadata = obs.targets.data.Properties.UserData;
             if ~ktransectSourceIdentity( ...
                   obs.targets.metadata, leg, entry.site_id) ...
                   || ~ktransectSourceIdentity(metadata, leg, entry.site_id)
                error(['icemodel:reconstruct:fillPromiceStation:' ...
                   'ktransectIdentityMismatch'], ...
                   ['K-transect donor source identity does not match ' ...
                   'manifest site %s.'], string(entry.site_id));
             end
             metadata_station = "";
             if isstruct(metadata) && isfield(metadata, 'station') ...
                   && (ischar(metadata.station) ...
                   || (isstring(metadata.station) && isscalar(metadata.station)))
                metadata_station = string(metadata.station);
             end
              if ~strcmpi(metadata_station, string(entry.site_id))
                error(['icemodel:reconstruct:fillPromiceStation:' ...
                   'ktransectIdentityMismatch'], ...
                   ['K-transect donor metadata station %s does not match ' ...
                   'manifest site %s.'], metadata_station, ...
                    string(entry.site_id));
              end
              location = ktransectLocation(entry.site_location, ...
                 obs.targets.metadata, metadata, entry.site_id);
              has_heights = isstruct(metadata) ...
               && isfield(metadata, 'sensor_heights') ...
               && isstruct(metadata.sensor_heights) ...
               && isscalar(metadata.sensor_heights) ...
               && isfield(metadata.sensor_heights, 'present') ...
               && isscalar(metadata.sensor_heights.present) ...
               && logical(metadata.sensor_heights.present) ...
               && isfield(metadata.sensor_heights, 'records') ...
               && ~isempty(metadata.sensor_heights.records);
            if ~has_heights
               % POLICY A3: height provenance annotates and downgrades a
               % donor record; it never blocks donor use, and one
               % heights-less artifact must never abort the whole fill.
               % The admission gates still judge the donor's transfers on
               % held-out evidence, so quality protection is preserved.
               warning(['icemodel:reconstruct:fillPromiceStation:' ...
                  'missingKtransectHeightProvenance'], ...
                  ['K-transect donor %s lacks measured sensor-height ' ...
                  'provenance; admitted without height annotation ' ...
                  '(POLICY A3)'], string(entry.site_id));
            end
             % Half-hourly K-transect donors aggregate to hourly BEFORE
             % transfer (DesignSpec Resolution 4): all donor channels here
             % are states or flux densities, so the hourly mean is the
             % correct per-variable-class rule.
             kt_series = hourlyStateMean(obs.targets.data);
             % K-transect albedo is the screened ratio of the station's own
             % coincident measured shortwave components. A ratio of two
             % simultaneous measurements is itself a measurement product
             % (all albedo is), so it remains donor-eligible under the
             % measurement-only donor rule (POLICY A8); only
             % reconstructed/gap-filled values are barred.
              kt_cells{k} = struct('series', kt_series, ...
                 'station', string(entry.site_id), 'family', "ktransect", ...
                 'location', location, 'observed_mask', []);
         end
         donor_cells = [donor_cells; kt_cells];
      end
   end

   % GC-Net Vandecrux stations: ONLY origin-observed samples are donors
   % (reconstructed samples never launder into new fills).
   if kwargs.use_gcnet
      gc_dir = fullfile(data_root, 'verification', 'gcnet');
      gc_hits = dir(fullfile(gc_dir, '*_surface.nc'));
      gc_cells = cell(numel(gc_hits), 1);
      for k = 1:numel(gc_hits)
         gc_cells{k} = icemodel.forcing.helpers.readGcnetDonor( ...
            fullfile(gc_hits(k).folder, gc_hits(k).name));
      end
      donor_cells = [donor_cells; gc_cells(~cellfun(@isempty, gc_cells))];
   end

    if isempty(donor_cells)
       donors = struct('series', {}, 'station', {}, 'family', {}, ...
          'location', {}, 'observed_mask', {});
    else
       donors = vertcat(donor_cells{:});
       donor_tokens = icemodel.forcing.helpers.normalizedFileToken( ...
          icemodel.forcing.helpers.gcnetVandecruxStation( ...
          string({donors.station})));
       target_token = icemodel.forcing.helpers.normalizedFileToken( ...
          icemodel.forcing.helpers.gcnetVandecruxStation(site));
       donors = donors(donor_tokens ~= target_token);
    end
end

function tf = ktransectSourceIdentity(metadata, leg, site_id)
   %KTRANSECTSOURCEIDENTITY Match staged donor metadata to manifest pins.
   required = ["station", "doi", "bundle_doi", "license", "children"];
   tf = isstruct(metadata) && isscalar(metadata) ...
      && all(isfield(metadata, required));
   if ~tf
      return
   end
   children = metadata.children;
   pinned = leg.children;
   station = string(metadata.station);
   doi = string(metadata.doi);
   bundle_doi = string(metadata.bundle_doi);
   license = string(metadata.license);
   pinned_doi = string(leg.doi);
   pinned_bundle = string(leg.bundle_doi);
   pinned_license = string(leg.license);
   site_id = string(site_id);
   scalar_text = all([isscalar(station), isscalar(doi), ...
      isscalar(bundle_doi), isscalar(license), isscalar(pinned_doi), ...
      isscalar(pinned_bundle), isscalar(pinned_license), isscalar(site_id)]);
   tf = scalar_text && isstruct(children) && ~isempty(children) ...
      && isstruct(pinned) && ~isempty(pinned) ...
      && all(isfield(children, ["year", "doi"])) ...
      && all(isfield(pinned, ["year", "doi"])) ...
      && strcmpi(station, site_id) ...
      && doi == pinned_doi && bundle_doi == pinned_bundle ...
      && license == pinned_license ...
      && isequal([children.year], [pinned.year]) ...
      && isequal(string({children.doi}), string({pinned.doi}));
end

function location = ktransectLocation(manifest_location, ...
      target_metadata, data_metadata, site_id)
   %KTRANSECTLOCATION Authenticate donor geometry against pinned artifact data.
   records = {manifest_location, struct(), struct()};
   if isstruct(target_metadata) && isscalar(target_metadata) ...
         && isfield(target_metadata, 'site_location')
      records{2} = target_metadata.site_location;
   end
   if isstruct(data_metadata) && isscalar(data_metadata) ...
         && isfield(data_metadata, 'site_location')
      records{3} = data_metadata.site_location;
   end

   fields = ["lat_wgs84", "lon_wgs84", "elev_m"];
   points = nan(3, numel(fields));
   for n = 1:numel(records)
      record = records{n};
      if ~isstruct(record) || ~isscalar(record) ...
            || ~all(isfield(record, fields))
         error(['icemodel:reconstruct:fillPromiceStation:' ...
            'ktransectIdentityMismatch'], ...
            'K-transect donor %s has no authenticated site location.', ...
            string(site_id));
      end
      points(n, :) = [double(record.lat_wgs84), ...
         double(record.lon_wgs84), double(record.elev_m)];
   end
   if any(~isfinite(points), 'all') ...
         || ~isequal(points(1, :), points(2, :)) ...
         || ~isequal(points(1, :), points(3, :))
      error(['icemodel:reconstruct:fillPromiceStation:' ...
         'ktransectIdentityMismatch'], ...
         ['K-transect donor %s manifest geometry does not match its ' ...
         'pinned artifact.'], string(site_id));
   end
   location = records{3};
end

function channel_methods = plannedMethods(plan)
   %PLANNEDMETHODS Convert plan channels into reconstructSeries input.
   n = numel(plan.channels);
   slots = cell(n, 1);
   for c = 1:n
      slots{c} = struct('channel', plan.channels(c).channel, ...
         'methods', plan.channels(c).methods);
   end
   channel_methods = vertcat(slots{:});
end

function proxies = loadStagedProxies(site, location, catalog, selected_files)
   %LOADSTAGEDPROXIES Load the canonical staged per-site RCM proxy met.
   % ACCEPTANCEWINDOW has already validated and pinned the per-source file
   % set: the widest anchor plus any span extenders (POLICY A6). Consume
   % only those pinned paths so discovery cannot drift from the producer
   % manifest. The anchor owns every sample inside its own span; extender
   % files contribute rows strictly outside it, restricted to the
   % anchor's channels so the source series keeps one schema.
   slots = cell(numel(catalog), 1);
   for k = 1:numel(catalog)
      prefix = "met_" + site + "_" + catalog(k).storage + "_";
      selected = false(numel(selected_files), 1);
      for h = 1:numel(selected_files)
         [~, name, ext] = fileparts(selected_files(h));
         selected(h) = startsWith(string(name), prefix) && ext == ".mat";
      end
      if ~any(selected)
         continue
      end
      source_files = selected_files(selected);
      loaded = cell(numel(source_files), 1);
      spans = NaT(numel(source_files), 2, 'TimeZone', 'UTC');
      for f = 1:numel(source_files)
         proxy_file = source_files(f);
         hit = dir(proxy_file);
         [proxy_series, loaded_file] = ...
            icemodel.forcing.reconstruct.loadWidestTimetable(hit);
         if isempty(proxy_series) || loaded_file ~= proxy_file
            error(['icemodel:reconstruct:fillPromiceStation:' ...
               'proxySelectionMismatch'], ...
               'selected proxy file is unavailable for %s: %s', ...
               site, proxy_file);
         end
         metadata = proxy_series.Properties.UserData;
         if ~icemodel.forcing.reconstruct.proxyArtifactIdentity( ...
               metadata, site, location, catalog(k).storage)
            error(['icemodel:reconstruct:fillPromiceStation:' ...
               'proxyIdentityMismatch'], ...
               ['staged proxy metadata does not identify target %s and ' ...
               'source %s: %s'], site, catalog(k).label, proxy_file);
         end
         loaded{f} = proxy_series;
         proxy_times = proxy_series.Properties.RowTimes;
         spans(f, :) = [proxy_times(1), proxy_times(end)];
      end
      % acceptanceWindow pins each source's anchor FIRST (its selection
      % order is the contract), so the first prefix match anchors here
      % without re-deriving a width metric that could rank near-ties
      % differently. Every other pinned file appends only samples not
      % already covered by anything merged so far — masking against the
      % anchor span alone would duplicate rows when two extenders
      % overlap each other outside it.
      anchor = 1;
      combined = loaded{anchor};
      anchor_names = string(combined.Properties.VariableNames);
      for f = 1:numel(source_files)
         if f == anchor
            continue
         end
         extender = loaded{f};
         outside = extender.Properties.RowTimes < spans(anchor, 1) ...
            | extender.Properties.RowTimes > spans(anchor, 2);
         outside = outside & ~ismember( ...
            extender.Properties.RowTimes, combined.Properties.RowTimes);
         if ~any(outside)
            continue
         end
         extender = extender(outside, :);
         % Channels the anchor lacks are dropped; channels the extender
         % lacks arrive missing, keeping one honest source schema.
         extra = timetable(extender.Properties.RowTimes);
         for name = anchor_names
            if ismember(name, ...
                  string(extender.Properties.VariableNames))
               extra.(name) = extender.(name);
            else
               extra.(name) = nan(height(extender), 1);
            end
         end
         combined = sortrows([combined; extra]);
      end
      % The merged source must carry exactly one row per posting; a
      % duplicate means the dedup contract above regressed, and silent
      % duplicates would let file order decide which value wins.
      if numel(unique(combined.Properties.RowTimes)) ~= height(combined)
         error(['icemodel:reconstruct:fillPromiceStation:' ...
            'duplicateProxyPostings'], ...
            'merged %s proxy series repeats postings for %s', ...
            catalog(k).label, site);
      end
      slots{k} = struct('series', combined, ...
         'name', catalog(k).code_name, ...
         'code_name', catalog(k).code_name);
   end
   present = ~cellfun(@isempty, slots);
   proxies = struct('series', {}, 'name', {}, 'code_name', {});
   if any(present)
      proxies = vertcat(slots{present});
   end
end

function [filled, provenance, audit] = adoptModisAlbedo( ...
      filled, provenance, audit, site, modis_dir, native, codes, opts)
    %ADOPTMODISALBEDO Fill residual albedo from staged MODIS daily userdata.
    % POLICY A11/B12 (activated by D-15): GEUS C6 daily albedo is an
    % albedo-only observational source ranking ahead of the RCM proxies in
    % the last-resort order. The staged artifact attaches through
    % icemodel.forcing.modisToMetCadence — the single daily->met-cadence
    % conversion path — so no second interpolation rule can drift. A site
    % with no staged artifact (or none with usable retrievals) simply
    % leaves the gap for the RCM tier; absence is not an error because
    % bedrock sites legitimately stage no_source_coverage artifacts.
    if ~ismember("albedo", string(filled.Properties.VariableNames))
       return
    end
    hits = dir(fullfile(modis_dir, sprintf('%s_modis_*_86400s.mat', site)));
    if isempty(hits)
       return
    end
    % The staging layer writes one artifact per site; the newest wins if a
    % regeneration ever leaves two.
    [~, newest] = max([hits.datenum]);
    S = load(fullfile(hits(newest).folder, hits(newest).name));
    if ~isfield(S, 'Data') || ~istimetable(S.Data) ...
          || ~ismember("albedo", string(S.Data.Properties.VariableNames))
       return
    end

    times = filled.Properties.RowTimes;
    [candidate, support] = icemodel.forcing.modisToMetCadence( ...
       S.Data.albedo, S.Data.Properties.RowTimes, times);
    needs = ~isfinite(filled.albedo);
    valid = support & icemodel.forcing.reconstruct.physicalValidity( ...
       "albedo", candidate, times);
    adopt = needs & valid;
    if ~any(adopt)
       return
    end

    % Seam-blend the adoption against the native record like every other
    % fallback, then keep only post-blend-valid samples.
    x = filled.albedo;
    full_candidate = x;
    full_candidate(adopt) = candidate(adopt);
    native_albedo = x;
    if ismember("albedo", string(native.Properties.VariableNames))
       native_albedo = native.albedo;
    end
    [blended, seam_note] = ...
       icemodel.forcing.reconstruct.blendFallbackSeams( ...
       times, native_albedo, x, full_candidate, adopt, ...
       jump_factor=opts.jump_factor, blend_hours=opts.blend_hours);
    target = find(adopt);
    post_valid = icemodel.forcing.reconstruct.physicalValidity( ...
       "albedo", blended(target), times(target));
    target = target(post_valid);
    if isempty(target)
       return
    end
    adopted = false(numel(times), 1);
    adopted(target) = true;
    x(target) = blended(target);
    filled.albedo = x;
    code = provenance.albedo;
    code(target) = codes.modis;
    provenance.albedo = code;
    rows = icemodel.forcing.reconstruct.auditSegments( ...
       times, adopted, "albedo", "modis:daily_albedo", ...
       "staged GEUS C6 daily albedo, linear to cadence" ...
       + string(seam_note));
    if ~isempty(rows)
       audit = [audit; cell2table(vertcat(rows{:}), ...
          'VariableNames', audit.Properties.VariableNames)];
    end
end

function [filled, provenance, audit] = adoptPrecip(filled, provenance, ...
      audit, proxies, codes, opts, native)
    %ADOPTPRECIP Adopt proxy total precipitation and its source phase split.
    % Proxy order is the policy order (MAR, then MERRA-2). NO partitioning
    % happens at reconstruction (POLICY A10/D-18): finite native components
    % are preserved; a missing complement derives by exact arithmetic from
    % the total; both-missing phases adopt the PROXY'S OWN split scaled to
    % the tapered total when the source provides one, and otherwise stay
    % missing for the runtime phase option to resolve.
    times = filled.Properties.RowTimes;
    required = icemodel.forcing.helpers.precipitationVariables();
    if ~all(ismember(required, string(filled.Properties.VariableNames)))
       return
    end

    ppt = filled.ppt;
    rain = filled.rainf;
    snow = filled.snowf;
    ppt_code = initialProvenance(provenance, "ppt", ppt, codes);
    rain_code = initialProvenance(provenance, "rainf", rain, codes);
    snow_code = initialProvenance(provenance, "snowf", snow, codes);

    % Normalize pre-existing precipitation before any adoption. Native rain
    % is the protected observation (A10), including a finite invalid value:
    % retain it so the publication boundary refuses the source defect rather
    % than silently replacing an observation. Invalid totals and snow phases
    % re-enter as missing so a valid proxy split can replace them.
    invalid_total = isfinite(ppt) ...
       & ~icemodel.forcing.reconstruct.scalarValidity("ppt", ppt);
    invalid_snow = isfinite(snow) & snow < 0;
    ppt(invalid_total) = NaN;
    snow(invalid_snow) = NaN;
    ppt_code(invalid_total) = codes.missing;
    snow_code(invalid_snow) = codes.missing;

    finite_total = icemodel.forcing.reconstruct.scalarValidity("ppt", ppt);
    rain_exceeds_total = finite_total & isfinite(rain) & rain > ppt;
    ppt(rain_exceeds_total) = NaN;
    ppt_code(rain_exceeds_total) = codes.missing;
    snow(rain_exceeds_total) = NaN;
    snow_code(rain_exceeds_total) = codes.missing;
    finite_total = icemodel.forcing.reconstruct.scalarValidity("ppt", ppt);
    snow_exceeds_total = finite_total & isfinite(snow) & snow > ppt;
    snow(snow_exceeds_total) = NaN;
    snow_code(snow_exceeds_total) = codes.missing;
    full_split = finite_total & isfinite(rain) & isfinite(snow);
    inconsistent = full_split ...
       & ~icemodel.forcing.helpers.precipitationConsistency( ...
       ppt, rain, snow);
    snow(inconsistent) = NaN;
    snow_code(inconsistent) = codes.missing;
    no_total_pair = ~isfinite(ppt) & isfinite(rain) & isfinite(snow);
    snow(no_total_pair) = NaN;
    snow_code(no_total_pair) = codes.missing;

    % A finite total plus ONE finite component determines the other by
    % exact arithmetic — that complement is bookkeeping, not partitioning,
    % so it is the only phase derivation reconstruction performs here
    % (POLICY A10/D-18). Both-missing phases wait for a proxy source split
    % or the runtime phase option.
    finite_total = icemodel.forcing.reconstruct.scalarValidity("ppt", ppt);
    rain_missing = ~isfinite(rain);
    snow_missing = ~isfinite(snow);
    candidate_rain = rain;
    candidate_snow = snow;
    only_rain = ~rain_missing & snow_missing;
    only_snow = rain_missing & ~snow_missing;
    candidate_snow(only_rain) = ppt(only_rain) - rain(only_rain);
    candidate_rain(only_snow) = ppt(only_snow) - snow(only_snow);
    complement = (only_rain | only_snow) & finite_total ...
       & candidate_rain >= 0 & candidate_snow >= 0;
    rain(complement & rain_missing) = candidate_rain(complement & ...
       rain_missing);
    snow(complement & snow_missing) = candidate_snow(complement & ...
       snow_missing);
    rain_code(complement & rain_missing) = ppt_code(complement & ...
       rain_missing);
    snow_code(complement & snow_missing) = ppt_code(complement & ...
       snow_missing);
    rain_rows = icemodel.forcing.reconstruct.auditSegments( ...
       times, complement & rain_missing, "rainf", ...
       "complement:total_minus_snow", "derived rain complement");
    snow_rows = icemodel.forcing.reconstruct.auditSegments( ...
       times, complement & snow_missing, "snowf", ...
       "complement:total_minus_rain", "derived snow complement");
    rows = [rain_rows; snow_rows];

    edge = diff([false; ~isfinite(ppt); false]);
    starts = find(edge == 1);
    stops = find(edge == -1) - 1;
    for g = 1:numel(starts)
       idx = (starts(g):stops(g)).';
       chosen = 0;
       fallback = 0;
       chosen_source = nan(numel(idx), 1);
       chosen_mask = false(numel(idx), 1);
       fallback_source = nan(numel(idx), 1);
       fallback_mask = false(numel(idx), 1);
       chosen_rain_proxy = nan(numel(idx), 1);
       chosen_snow_proxy = nan(numel(idx), 1);
       fallback_rain_proxy = nan(numel(idx), 1);
       fallback_snow_proxy = nan(numel(idx), 1);
       for p = 1:numel(proxies)
          [source, proxy_rain, proxy_snow, adopt] = ...
             precipProxyPlan(proxies(p), times, idx, rain, snow);
          if any(adopt) && fallback == 0
             fallback = p;
             fallback_source = source;
             fallback_mask = adopt;
             fallback_rain_proxy = proxy_rain;
             fallback_snow_proxy = proxy_snow;
          end
          if all(adopt)
             chosen = p;
             chosen_source = source;
             chosen_mask = adopt;
             chosen_rain_proxy = proxy_rain;
             chosen_snow_proxy = proxy_snow;
             break
          end
       end
       if chosen == 0
          chosen = fallback;
          chosen_source = fallback_source;
          chosen_mask = fallback_mask;
          chosen_rain_proxy = fallback_rain_proxy;
          chosen_snow_proxy = fallback_snow_proxy;
       end
       if chosen == 0
          continue
       end

       % One source owns the whole outage; later proxies never fill leftovers.
       % Taper its total before resolving the source split so the phase
       % identity remains exact.
       adopted = false(numel(times), 1);
       adopted(idx(chosen_mask)) = true;
       candidate = ppt;
       candidate(adopted) = chosen_source(chosen_mask);
       [candidate, seam_note] = ...
          icemodel.forcing.reconstruct.blendFallbackSeams( ...
          times, native.ppt, ppt, candidate, adopted, ...
          jump_factor=opts.jump_factor, blend_hours=opts.blend_hours);
       chosen_source(chosen_mask) = candidate(adopted);
       post_blend_conflict = chosen_mask ...
          & precipitationPhaseConflict( ...
          chosen_source, rain(idx), snow(idx));
       chosen_mask(post_blend_conflict) = false;
       [chosen_rain, chosen_snow, phase_known] = sourceSplit( ...
          chosen_source, chosen_rain_proxy, chosen_snow_proxy, ...
          rain(idx), snow(idx));
       target = idx(chosen_mask);
       rain_missing = ~isfinite(rain(idx));
       snow_missing = ~isfinite(snow(idx));
       % Phases fill only where the source split (or a complement) exists;
       % the total still adopts on its own so the runtime phase option can
       % resolve any phase-unknown samples (POLICY A10/D-18).
       derived_rain = chosen_mask & rain_missing & phase_known;
       derived_snow = chosen_mask & snow_missing & phase_known;
       source_code = codes.(char(proxies(chosen).code_name));
       ppt(target) = chosen_source(chosen_mask);
       ppt_code(target) = source_code;
       rain(idx(derived_rain)) = chosen_rain(derived_rain);
       snow(idx(derived_snow)) = chosen_snow(derived_snow);
       rain_code(idx(derived_rain)) = source_code;
       snow_code(idx(derived_snow)) = source_code;

       % Every filled precipitation channel gets contiguous audit rows.
       for item = {["ppt", "total"], ["rainf", "derived rain"], ...
             ["snowf", "derived snow"]}
          pair = item{1};
          local_mask = chosen_mask;
          if pair(1) == "rainf"
             local_mask = derived_rain;
          elseif pair(1) == "snowf"
             local_mask = derived_snow;
          end
          mask = false(numel(times), 1);
          mask(idx(local_mask)) = true;
           segment_rows = icemodel.forcing.reconstruct.auditSegments( ...
              times, mask, pair(1), ...
              "proxy:" + proxies(chosen).name + ":precip_adoption", ...
              pair(2) + string(seam_note));
          rows = [rows; segment_rows]; %#ok<AGROW>
       end
    end

    filled.ppt = ppt;
    filled.rainf = rain;
    filled.snowf = snow;
    provenance.ppt = ppt_code;
    provenance.rainf = rain_code;
    provenance.snowf = snow_code;
    % Optional phase gaps are legitimate at runtime, but they are still
    % shipped missing values and the final report must explain every one.
    % Record the post-adoption state directly because precipitation never
    % enters the statistically planned provisional-unfilled ledger.
    for item = { ...
          ["ppt", "no valid staged proxy total precipitation value"], ...
          ["rainf", "source rain phase unavailable; runtime phase option required"], ...
          ["snowf", "source snow phase unavailable; runtime phase option required"]}
       pair = item{1};
       name = pair(1);
       values = filled.(name);
       missing = ~isfinite(values);
       segment_rows = icemodel.forcing.reconstruct.auditSegments( ...
          times, missing, name, "unfilled", ...
          "final tier: " + pair(2), context_id="precip_adoption");
       rows = [rows; segment_rows]; %#ok<AGROW>
    end
    if ~isempty(rows)
       audit = [audit; cell2table(vertcat(rows{:}), ...
          'VariableNames', audit.Properties.VariableNames)];
    end
end

function [source, proxy_rain, proxy_snow, adopt] = ...
      precipProxyPlan(proxy, times, idx, rain, snow)
    %PRECIPPROXYPLAN Test one proxy as the sole source for one ppt outage.
    % Returns the aligned total plus the proxy's own phase split (NaN when
    % the source ships no split; a rain channel absent beside ppt+snowf
    % derives as the exact complement per D-31) so adoption can rescale
    % the source fraction instead of partitioning (POLICY A10/D-18).
    aligned = alignProxyChannel(proxy, "ppt", times);
    source = aligned(idx);
    aligned_rain = alignProxyChannel(proxy, "rainf", times);
    proxy_rain = aligned_rain(idx);
    aligned_snow = alignProxyChannel(proxy, "snowf", times);
    proxy_snow = aligned_snow(idx);
    % Staged MERRA-2 ships the total (PRECTOTCORR -> ppt) and snowfall
    % (PRECSNOCORR -> snowf) but no rain channel — that is the staged
    % channel inventory, not missing physics. The source's own split is
    % still fully determined, so derive the missing component as the
    % exact complement ppt - snowf, floored at zero against sub-epsilon
    % source inconsistency, so adoptions ship a full split instead of
    % phase-unknown samples (POLICY D-31). Subtraction propagates NaN,
    % so samples without both operands stay phase-unknown.
    proxy_names = string(proxy.series.Properties.VariableNames);
    if ~ismember("rainf", proxy_names) ...
          && all(ismember(["ppt", "snowf"], proxy_names))
       complement = source - proxy_snow;
       complement(complement < 0) = 0;
       proxy_rain = complement;
    end
    valid = icemodel.forcing.reconstruct.physicalValidity( ...
       "ppt", source, times(idx));
    % A finite native component larger than the candidate total is a
    % conflict with the observation. Two finite native phases must also
    % sum to the candidate total; otherwise preserving them beside that
    % total would publish an impossible phase identity. The observations
    % win and that sample refuses adoption (the C24 native-preserving veto).
    rain_native = rain(idx);
    snow_native = snow(idx);
    conflict = precipitationPhaseConflict( ...
       source, rain_native, snow_native);
    adopt = valid & ~conflict;
end

function [candidate_rain, candidate_snow, phase_known] = ...
      sourceSplit(total, proxy_rain, proxy_snow, rain_native, snow_native)
    %SOURCESPLIT Resolve phase components without any partitioning.
    % Order of truth: finite native components; then the exact complement
    % from the adopted total; then the proxy's OWN split rescaled so its
    % phase FRACTION survives the seam taper of the total. Samples with
    % none of those stay phase-unknown for the runtime option
    % (POLICY A10/D-18).
    candidate_rain = rain_native;
    candidate_snow = snow_native;
    rain_missing = ~isfinite(rain_native);
    snow_missing = ~isfinite(snow_native);
    only_rain = ~rain_missing & snow_missing;
    only_snow = rain_missing & ~snow_missing;
    both_missing = rain_missing & snow_missing;
    candidate_snow(only_rain) = total(only_rain) - rain_native(only_rain);
    candidate_rain(only_snow) = total(only_snow) - snow_native(only_snow);

    % Rescale the source split onto the tapered total: the source's phase
    % fraction is physical information; its absolute magnitudes are not
    % once the seam taper adjusts the total.
    proxy_total = proxy_rain + proxy_snow;
    scalable = both_missing & isfinite(proxy_total) & proxy_total > 0 ...
       & isfinite(total);
    frac = proxy_rain(scalable) ./ proxy_total(scalable);
    candidate_rain(scalable) = frac .* total(scalable);
    candidate_snow(scalable) = (1 - frac) .* total(scalable);
    % A zero source total with a finite adopted total carries no fraction
    % information; a zero adopted total splits to exact zeros.
    zero_total = both_missing & isfinite(total) & total == 0;
    candidate_rain(zero_total) = 0;
    candidate_snow(zero_total) = 0;

    phase_known = isfinite(candidate_rain) & isfinite(candidate_snow) ...
       & candidate_rain >= 0 & candidate_snow >= 0;
end

function conflict = precipitationPhaseConflict(total, rain, snow)
    %PRECIPITATIONPHASECONFLICT Protect finite native phase observations.
    both = isfinite(rain) & isfinite(snow);
    consistent = icemodel.forcing.helpers.precipitationConsistency( ...
       total, rain, snow);
    conflict = (isfinite(rain) & rain > total) ...
       | (isfinite(snow) & snow > total) ...
       | (both & ~consistent);
end

function code = initialProvenance(provenance, channel, values, codes)
    %INITIALPROVENANCE Preserve existing codes or initialize native/missing.
    if ismember(channel, string(provenance.Properties.VariableNames))
       code = provenance.(channel);
       return
    end
    code = repmat(codes.missing, numel(values), 1);
    code(isfinite(values)) = codes.observed;
end

function source = alignProxyChannel(proxy, channel, times)
    %ALIGNPROXYCHANNEL Align one staged proxy channel to the target axis.
    source = nan(numel(times), 1);
    if ~ismember(channel, string(proxy.series.Properties.VariableNames))
       return
    end
    [tf, loc] = ismember(times, proxy.series.Properties.RowTimes);
    source(tf) = proxy.series.(channel)(loc(tf));
end

function [filled, provenance, audit] = fillBoomHeight( ...
      filled, provenance, audit, codes)
   %FILLBOOMHEIGHT Preserve boom-height gaps without complete visit boundaries.
   % A native product with no boom channel at all still ships one: the
   % runtime loader requires the channel's presence for provenance
   % integrity, and an all-missing column lands every sample on the A3
   % nominal rung instead of blocking the load (POLICY A3).
   if ~ismember("boom_height", ...
         string(filled.Properties.VariableNames))
      filled.boom_height = nan(height(filled), 1);
   end

   % Station-transition metadata records handovers, not every maintenance
   % visit. Without the latter, no gap can be proven to stay within one sensor
   % mounting interval, so fail closed and stamp only native observations.
   % Preserving the gaps never blocks a run: the runtime POLICY A3 chain in
   % icemodel.loadmet resolves them (measured -> interpolated -> nominal
   % 2.6 m) with the outcome recorded in opts.boom_height_source.
   native = filled.boom_height;
   provenance.boom_height = ...
      initialProvenance(provenance, "boom_height", native, codes);
end

function readiness = readinessLedger(site, native, filled, plan, opts, location)
   %READINESSLEDGER Per-year forcing verdicts and confidence advisories.
   % Two absolute verdicts per station-year (POLICY A5/D-16):
   % ready_icemodel grades the icemodel required set; ready_snowmodel
   % additionally requires snowfall input (finite total ppt OR snowf per
   % sample). Verdicts grade COMPLETENESS plus the SCALAR registry
   % bounds only (finite and inside physicalBounds); RELATIONAL
   % exceedances — swd above the TOA ceiling, swu above swd — are
   % reported in the worst_relational_invalid and relational_diagnostic
   % columns and NEVER flip a verdict (POLICY A15/D-28: twilight diffuse
   % light is real incident energy the TOA model cannot represent, and
   % native sensor noise must not brand usable years not-ready; the
   % relational rules stay hard gates for fill CANDIDATES in
   % physicalValidity). A required channel absent from the product is
   % wholly missing, not a channel to skip. Boom height never grades a
   % verdict (POLICY A3): the runtime fallback chain guarantees usable
   % heights, so geometry is reported nowhere here. Phase consistency is
   % graded only where all three precipitation channels are finite
   % (partial phases are legitimate; runtime resolves them).
   required = ...
      icemodel.forcing.reconstruct.icemodelRequiredChannels();
   present = intersect(required, ...
      string(filled.Properties.VariableNames), 'stable');
   times = filled.Properties.RowTimes;
   names = string(filled.Properties.VariableNames);
   precip_names = icemodel.forcing.helpers.precipitationVariables();

   % Snowfall-input availability per sample, and phase-identity validity
   % over the samples where a full split exists.
   snow_input = false(size(times));
   phase_conflict = false(size(times));
   % Only physically valid totals or snowfall count as snowfall input; both
   % share ppt's scalar registry because they are accumulation rates.
   if ismember("ppt", names)
      snow_input = snow_input ...
         | icemodel.forcing.reconstruct.scalarValidity("ppt", filled.ppt);
   end
   if ismember("snowf", names)
      snow_input = snow_input ...
         | icemodel.forcing.reconstruct.scalarValidity("ppt", filled.snowf);
   end
   if all(ismember(precip_names, names))
      full_split = isfinite(filled.ppt) & isfinite(filled.rainf) ...
         & isfinite(filled.snowf);
      phase_conflict = full_split ...
         & ~icemodel.forcing.helpers.precipitationConsistency( ...
         filled.ppt, filled.rainf, filled.snowf);
   end

   record_years = unique(year(times));
   rows = cell(numel(record_years), 1);
   for y = 1:numel(record_years)
      in_year = year(times) == record_years(y);
      core = intersect(opts.core_channels, ...
         string(native.Properties.VariableNames), 'stable');
      native_core = 0;
      if numel(core) == numel(opts.core_channels)
         native_core = mean(all(isfinite(native{in_year, ...
            cellstr(core)}), 2));
      end

      % Grade the icemodel set once; the snowmodel verdict reuses it and
      % appends the snowfall-input and phase-identity requirements.
      % Grading is deliberately SCALAR-only (A15/D-28): completeness plus
      % the physicalBounds registry, never the relational rules.
      invalid_after = strings(numel(required) + 2, 1);
      n_invalid = 0;
      worst_residual = 0;
      for name = required
         if ismember(name, present)
            residual = mean(scalarInvalid(name, filled.(name)(in_year)));
         else
            % Absent required channel = 100% missing for every year.
            residual = 1;
         end
         worst_residual = max(worst_residual, residual);
         if residual > 0
            n_invalid = n_invalid + 1;
            invalid_after(n_invalid) = ...
               sprintf("%s %.1f%% invalid", name, 100 * residual);
         end
      end
      ice_reasons = invalid_after(1:n_invalid);

      % Relational DIAGNOSTICS (A15/D-28): the fraction of scalar-valid
      % samples that exceed the paired rule, reported beside the
      % verdicts but never allowed to flip one. The swd check reuses
      % physicalValidity so the diagnostic and the candidate gate can
      % never disagree about the ceiling.
      relational_notes = strings(2, 1);
      n_relational = 0;
      worst_relational = 0;
      if ismember("swd", required) && ismember("swd", present)
         values = filled.swd(in_year);
          ceiling_ok = icemodel.forcing.reconstruct.physicalValidity( ...
             "swd", values, times(in_year), ...
             latitude=location.lat_wgs84, longitude=location.lon_wgs84, ...
             interval=median(diff(times)));
         frac = mean(~scalarInvalid("swd", values) & ~ceiling_ok);
         worst_relational = max(worst_relational, frac);
         if frac > 0
            n_relational = n_relational + 1;
            relational_notes(n_relational) = sprintf( ...
               "swd %.1f%% above TOA ceiling", 100 * frac);
         end
      end
      % swu is derived and never required (A5/B10), so the diagnostic
      % gates on PRODUCT presence — the required-set gate left this
      % branch dead in every policy-default run (review pass 9).
      if ismember("swu", names)
         values = filled.swu(in_year);
         swd_reference = nan(size(values));
         if ismember("swd", names)
            swd_reference = filled.swd(in_year);
         end
         % A finite swu without a finite swd reference (or above it)
         % counts as a relational exceedance: the pairing cannot be
         % proven, which is exactly what the diagnostic reports.
         exceed = ~scalarInvalid("swu", values) ...
            & ~(isfinite(swd_reference) & values <= swd_reference);
         frac = mean(exceed);
         worst_relational = max(worst_relational, frac);
         if frac > 0
            n_relational = n_relational + 1;
            relational_notes(n_relational) = sprintf( ...
               "swu %.1f%% above swd", 100 * frac);
         end
      end
      % Guard the empty case explicitly: a clean year reports "" rather
      % than relying on strjoin's empty-input behavior.
      relational_diagnostic = "";
      if n_relational > 0
         relational_diagnostic = strjoin( ...
            relational_notes(1:n_relational), "; ");
      end

      % The snowmodel verdict reuses the icemodel reasons and appends at
      % most two of its own, so the extras preallocate to that maximum.
      snow_extra = strings(2, 1);
      n_extra = 0;
      n_year_samples = nnz(in_year);
      n_snow_missing = nnz(~snow_input(in_year));
      snow_missing = n_snow_missing / n_year_samples;
      if snow_missing > 0
         n_extra = n_extra + 1;
         snow_extra(n_extra) = sprintf( ...
            "snowfall input %d/%d samples (%.4g%%) missing", ...
            n_snow_missing, n_year_samples, 100 * snow_missing);
      end
      phase_bad = mean(phase_conflict(in_year));
      if phase_bad > 0
         n_extra = n_extra + 1;
         snow_extra(n_extra) = sprintf( ...
            "precipitation phase %.1f%% inconsistent", 100 * phase_bad);
      end
      snow_reasons = [ice_reasons; snow_extra(1:n_extra)];
      snow_worst = max([worst_residual, snow_missing, phase_bad]);

      % Readiness grades the shipped forcing, never how much native data
      % happened to constrain it. Sparse native support remains a visible
      % scientific-confidence advisory, but it cannot excuse an unfilled
      % in-window product or demote a complete one (POLICY A13/D-34).
      verdict_ice = "ready";
      reason_ice = "";
      if ~isempty(ice_reasons)
         verdict_ice = "not_forcing_ready";
         reason_ice = strjoin(ice_reasons, "; ");
      end
      verdict_snow = "ready";
      reason_snow = "";
      if ~isempty(snow_reasons)
         verdict_snow = "not_forcing_ready";
         reason_snow = strjoin(snow_reasons, "; ");
      end
      worth_filling = native_core >= opts.min_native_core_coverage;
      triage_note = "";
      if ~worth_filling
         triage_note = sprintf( ...
            "native core coverage %.0f%% below %.0f%% advisory", ...
            100 * native_core, 100 * opts.min_native_core_coverage);
      end
      rows{y} = {char(site), record_years(y), native_core, ...
         worst_residual, char(verdict_ice), char(reason_ice), ...
         snow_worst, char(verdict_snow), char(reason_snow), ...
         worst_relational, char(relational_diagnostic), worth_filling, ...
         char(triage_note)};
   end
   readiness = cell2table(vertcat(rows{:}), 'VariableNames', ...
      {'site', 'year', 'native_core_coverage', ...
      'worst_residual_invalid', 'verdict_icemodel', 'reason_icemodel', ...
      'worst_residual_snowmodel', 'verdict_snowmodel', ...
      'reason_snowmodel', 'worst_relational_invalid', ...
      'relational_diagnostic', 'worth_filling', 'triage_note'});
   readiness.n_admitted_methods = repmat(sum(arrayfun(@(c) ...
      numel(c.methods), plan.channels)), height(readiness), 1);
end

function invalid = scalarInvalid(channel, values)
   %SCALARINVALID Finite-and-in-scalar-bounds failure mask for one channel.
   % Verdict grading is deliberately scalar-only (POLICY A15/D-28), and
   % the shared registry keeps the bound values single-sourced.
   invalid = ~icemodel.forcing.reconstruct.scalarValidity(channel, values);
end

function met_file = writeArtifacts(site, filled, provenance, audit, ...
      plan, readiness, seam_quality, flat_run_findings, out_dir, qa_dir, ...
      codes, native_file, ...
      acceptance_window, proxy_window_files, split_manifest, split_file, ...
      donor_sites, family)
   %WRITEARTIFACTS Transactionally persist met, provenance, plan, and ledger.
   % The family token names the shipped product <family>_filled — the
   % writemet source token, the gapfill_product identity stamp, and the
   % superseded-window glob all derive from it (bead icemodel-g1n.49).
   icemodel.helpers.ensureDirExists(out_dir);
   icemodel.helpers.ensureDirExists(fullfile(qa_dir, 'plans'));
   icemodel.helpers.ensureDirExists(fullfile(qa_dir, 'ledger'));
   icemodel.helpers.ensureDirExists(fullfile(qa_dir, 'splits'));
   stage_root = string(tempname(fileparts(char(out_dir))));
   mkdir(stage_root);
   stage_cleanup = onCleanup(@() removeArtifactTree(stage_root));
   stage_met_root = fullfile(stage_root, 'met');
   stage_qa_dir = fullfile(stage_root, 'qa');
   icemodel.helpers.ensureDirExists(stage_met_root);
   icemodel.helpers.ensureDirExists(fullfile(stage_qa_dir, 'plans'));
   icemodel.helpers.ensureDirExists(fullfile(stage_qa_dir, 'ledger'));
   icemodel.helpers.ensureDirExists(fullfile(stage_qa_dir, 'splits'));

   % Inline provenance channels ride the met artifact (registered names),
   % and the registry/policy identity lands in UserData beside the existing
   % artifact metadata. Attachment is union-driven (POLICY A7 provenance
   % honesty): every provenance column whose channel ships in the met
   % timetable ships beside it, so channels stamped outside plan.channels
   % (adopted precipitation, boom_height, and the driver-created swu that
   % deriveUpwardShortwave fills) always carry their codes instead of
   % depending on per-list attachment that once skipped them.
   met = filled;
   for name = string(provenance.Properties.VariableNames)
      if ismember(name, string(met.Properties.VariableNames))
         met.(name + "_provenance") = provenance.(name);
      end
   end
   ud = met.Properties.UserData;
     % Runtime identity uses the canonical compact token, not the
     % separator-bearing display name inherited from native metadata.
     ud.site = site;
     ud.gapfill_registry = codes;
     ud.gapfill_seed = plan.split.seed;
     ud.gapfill_product = char(family + "_filled");
     ud.gapfill_channels = string({plan.channels.channel});
     ud.gapfill_engine_version = string(icemodel.internal.version());
    policy_file = fullfile(fileparts(mfilename('fullpath')), 'POLICY.md');
    ud.gapfill_policy_sha256 = ...
       icemodel.verification.setup.fileSha256(policy_file);
    ud.gapfill_donors = donor_sites(:).';
    met.Properties.UserData = ud;

   % Stage through the canonical writer so metadata and window naming match
   % every other met artifact without exposing a partial final product.
   met = icemodel.forcing.helpers.stampMetadata(met);
   met_files = icemodel.forcing.helpers.writemet(met, site, ...
      family + "_filled", outdir=stage_met_root, naming="window", ...
      dt_out="", overwrite=true);
   staged_met_file = string(met_files(1));
   [~, met_name, met_ext] = fileparts(staged_met_file);
   met_file = fullfile(out_dir, string(met_name) + string(met_ext));

   % Plan summary (JSON; fitted parameter structs and estimates stay in the
   % MAT sidecar to keep the JSON reviewable) plus audit and ledger rows.
   plan_record = stripEstimates(plan);
   audit_record = audit;
   seam_quality_record = seam_quality;
   flat_run_findings_record = flat_run_findings;
   plan_file = fullfile(qa_dir, 'plans', sprintf('%s-plan.mat', site));
   staged_plan_file = fullfile(stage_qa_dir, 'plans', ...
      sprintf('%s-plan.mat', site));
   save(staged_plan_file, 'plan_record', 'audit_record', ...
      'seam_quality_record', 'flat_run_findings_record');
   summary = planSummary(plan);
   summary_file = fullfile(qa_dir, 'plans', ...
      sprintf('%s-plan-summary.json', site));
   staged_summary_file = fullfile(stage_qa_dir, 'plans', ...
      sprintf('%s-plan-summary.json', site));
   fid = fopen(staged_summary_file, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fprintf(fid, '%s', jsonencode(summary));
   clear cleaner
   readiness_file = fullfile(qa_dir, 'ledger', ...
      sprintf('%s-readiness.csv', site));
   staged_readiness_file = fullfile(stage_qa_dir, 'ledger', ...
      sprintf('%s-readiness.csv', site));
   writetable(readiness, staged_readiness_file);
   staged_split_file = fullfile(stage_qa_dir, 'splits', ...
      sprintf('%s-split.json', site));
   if ~isfile(split_manifest)
      error('icemodel:reconstruct:fillPromiceStation:missingSplitManifest', ...
         'station plan did not produce its staged split manifest: %s', ...
         split_manifest);
   end
   copyfile(split_manifest, staged_split_file);

   % Freeze every direct report input plus each proxy filename that defines
   % the acceptance-window verdict. Hash staged bytes under their final paths.
   n_proxy = numel(proxy_window_files);
   roles = ["native"; "filled"; "plan"; "readiness"; ...
      repmat("proxy_window", n_proxy, 1)];
    final_files = [string(native_file); string(met_file); ...
       string(plan_file); string(readiness_file); proxy_window_files(:)];
    hash_files = [string(native_file); staged_met_file; ...
       string(staged_plan_file); string(staged_readiness_file); ...
       proxy_window_files(:)];
    [data_root, ~] = ...
       icemodel.forcing.reconstruct.selectedDataRoot(string(out_dir));
    for k = 1:numel(final_files)
       if ~icemodel.internal.isPathInside(final_files(k), data_root)
          error('icemodel:reconstruct:fillPromiceStation:artifactOutsideRoot', ...
             'report input must stay inside selected data root %s: %s', ...
             data_root, final_files(k));
       end
    end
    relative_files = icemodel.verification.setup.relpaths( ...
       final_files, data_root);
    artifacts = repmat(struct('role', "", 'path', "", 'bytes', 0, ...
       'sha256', ""), numel(final_files), 1);
    for k = 1:numel(final_files)
       info = dir(hash_files(k));
       artifacts(k) = struct('role', roles(k), 'path', relative_files(k), ...
          'bytes', info.bytes, 'sha256', ...
          icemodel.verification.setup.fileSha256(hash_files(k)));
   end
   window = struct('start', ...
      icemodel.verification.setup.formatManifestTime( ...
      acceptance_window(1)), 'end', ...
      icemodel.verification.setup.formatManifestTime( ...
      acceptance_window(2)));
    manifest = struct('site', site, ...
       'path_base', "selected_data_root", 'artifacts', artifacts, ...
       'acceptance_window', window);
   manifest_file = fullfile(qa_dir, 'plans', ...
      sprintf('%s-report-inputs.json', site));
   staged_manifest_file = fullfile(stage_qa_dir, 'plans', ...
      sprintf('%s-report-inputs.json', site));
   fid = fopen(staged_manifest_file, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fprintf(fid, '%s', jsonencode(manifest));
   clear cleaner

   % Back up every destination, install the complete staged set, and restore
   % all prior artifacts if any move fails. Superseded met windows participate
   % in the same transaction rather than being pruned early by writemet.
    sources = [staged_met_file; string(staged_plan_file); ...
       string(staged_summary_file); string(staged_readiness_file); ...
       string(staged_manifest_file); string(staged_split_file)];
    destinations = [string(met_file); string(plan_file); ...
       string(summary_file); string(readiness_file); string(manifest_file); ...
       string(split_file)];
   prior = dir(fullfile(out_dir, sprintf( ...
      'met_%s_%s_filled_*.mat', site, family)));
   stale = string(fullfile({prior.folder}, {prior.name})).';
   stale = setdiff(stale, string(met_file), 'stable');
   publishArtifacts(sources, destinations, stale, fileparts(out_dir));
   clear stage_cleanup
end

function publishArtifacts(sources, destinations, stale, backup_parent)
   %PUBLISHARTIFACTS Install one complete product set with rollback.
   all_destinations = [destinations(:); stale(:)];
   backup_root = string(tempname(backup_parent));
   mkdir(backup_root);
   backed = false(numel(all_destinations), 1);
   backups = strings(numel(all_destinations), 1);
   installed = false(numel(sources), 1);
   try
      % Preserve every prior file before installing any replacement.
      for k = 1:numel(all_destinations)
         if isfile(all_destinations(k)) || isfolder(all_destinations(k))
            backups(k) = fullfile(backup_root, sprintf('%06d', k));
            [ok, message] = movefile( ...
               all_destinations(k), backups(k));
            if ~ok
               error(['icemodel:reconstruct:fillPromiceStation:' ...
                  'backupFailed'], 'could not back up %s: %s', ...
                  all_destinations(k), message);
            end
            backed(k) = true;
         end
      end
      % Install only after the complete prior set is recoverable.
      for k = 1:numel(sources)
         [ok, message] = movefile(sources(k), destinations(k));
         if ~ok
            error(['icemodel:reconstruct:fillPromiceStation:' ...
               'installFailed'], 'could not install %s: %s', ...
               destinations(k), message);
         end
         installed(k) = true;
      end
   catch err
      % Remove only artifacts installed by this attempt.
      for k = numel(installed):-1:1
         if installed(k)
            removeArtifactPath(destinations(k));
         end
      end
      % Restore the complete prior set; retain backups if restoration fails.
      restored = true;
      restore_message = "";
      for k = numel(backed):-1:1
         if backed(k)
            [ok, message] = movefile( ...
               backups(k), all_destinations(k));
            restored = restored && ok;
            if ~ok
               restore_message = restore_message + newline + message;
            end
         end
      end
      if restored
         removeArtifactTree(backup_root);
         rethrow(err)
      end
      error(['icemodel:reconstruct:fillPromiceStation:' ...
         'rollbackFailed'], ...
         'artifact publish failed; recovery files remain under %s:%s', ...
         backup_root, restore_message);
   end
   removeArtifactTree(backup_root);
end

function retirePublishedArtifacts(site, out_dir, qa_dir, family)
   %RETIREPUBLISHEDARTIFACTS Remove one invalidated station product set.
   met = dir(fullfile(out_dir, sprintf( ...
      'met_%s_%s_filled_*.mat', site, family)));
   candidates = [ ...
      string(fullfile({met.folder}, {met.name})).'; ...
      string(fullfile(qa_dir, 'plans', sprintf('%s-plan.mat', site))); ...
      string(fullfile(qa_dir, 'plans', ...
      sprintf('%s-plan-summary.json', site))); ...
      string(fullfile(qa_dir, 'plans', ...
      sprintf('%s-report-inputs.json', site))); ...
      string(fullfile(qa_dir, 'ledger', ...
      sprintf('%s-readiness.csv', site))); ...
      string(fullfile(qa_dir, 'splits', sprintf('%s-split.json', site)))];
   existing = unique(candidates(arrayfun(@isfile, candidates)), 'stable');
   if isempty(existing)
      return
   end
   publishArtifacts(strings(0, 1), strings(0, 1), existing, ...
      fileparts(out_dir));
end

function removeArtifactPath(pathname)
   %REMOVEARTIFACTPATH Remove one explicitly enumerated installed artifact.
   if isfolder(pathname)
      rmdir(pathname, 's');
   elseif isfile(pathname)
      delete(pathname);
   end
end

function removeArtifactTree(pathname)
   %REMOVEARTIFACTTREE Remove one isolated staging or backup directory.
   if isfolder(pathname)
      rmdir(pathname, 's');
   end
end

function s = stripEstimates(s)
   %STRIPESTIMATES Drop bulky estimate vectors from the persisted plan.
   for c = 1:numel(s.channels)
      for m = 1:numel(s.channels(c).methods)
         s.channels(c).methods(m).estimate = [];
      end
   end
end

function summary = planSummary(plan)
   %PLANSUMMARY Reviewable JSON summary of the admitted methods.
   channels = cell(numel(plan.channels), 1);
   for c = 1:numel(plan.channels)
      methods = plan.channels(c).methods;
      rows = cell(numel(methods), 1);
      for m = 1:numel(methods)
         rows{m} = struct('name', methods(m).name, ...
            'code', methods(m).code, 'buckets', methods(m).buckets, ...
            'seasons', methods(m).seasons, ...
            'audit_context_id', methods(m).audit_context_id, ...
            'max_validated_hours', methods(m).max_validated_hours, ...
            'selection_rmse', methods(m).selection_rmse);
      end
      channels{c} = struct('channel', plan.channels(c).channel, ...
         'methods', {rows});
   end
   summary = struct('station', plan.station, ...
      'seed', plan.split.seed, ...
      'years_selection', plan.split.years_selection, ...
      'years_evaluation', plan.split.years_evaluation, ...
      'fixed_methods', plan.fixed_methods, ...
      'audit_contexts', table2struct(plan.audit_contexts), ...
      'channels', {channels});
end

function records = auditContextRecords(plan, audit)
   %AUDITCONTEXTRECORDS Bind every emitted audit context to the saved plan.
   context_id = unique(string(audit.context_id), 'stable');
   method = strings(numel(context_id), 1);
   record_type = repmat("runtime_policy", numel(context_id), 1);
   channels = strings(numel(context_id), 1);
   n_segments = zeros(numel(context_id), 1);
   candidate_ids = strings(0, 1);
   for c = 1:numel(plan.channels)
      candidate_ids = [candidate_ids; string( ...
         {plan.channels(c).methods.audit_context_id}).']; %#ok<AGROW>
   end
   fixed_ids = string({plan.fixed_methods.audit_context_id}).';
   for k = 1:numel(context_id)
      rows = string(audit.context_id) == context_id(k);
      methods = unique(string(audit.method(rows)));
      if strlength(context_id(k)) == 0 || numel(methods) ~= 1
         error( ...
            'icemodel:reconstruct:fillPromiceStation:auditContextConflict', ...
            'audit context must be nonblank and identify one method: %s', ...
            context_id(k));
      end
      method(k) = methods;
      channels(k) = strjoin(unique(string(audit.channel(rows)), ...
         'stable'), ",");
      n_segments(k) = nnz(rows);
      if ismember(context_id(k), candidate_ids)
         record_type(k) = "fitted_method";
      elseif ismember(context_id(k), fixed_ids)
         record_type(k) = "fixed_policy";
      end
   end
   records = table(context_id, method, record_type, channels, n_segments);
end

function donor_sites = inventoryDonors(data_root, site, opts, donor_families)
   %INVENTORYDONORS Geometry-plausible donor-family neighbors from the
   % inventory. The registry records each neighbor's family; only members
   % of the caller's donor family set are admitted (bead icemodel-g1n.49).
   donor_sites = strings(1, 0);
   inventory_file = fullfile(data_root, 'preview', 'qa', 'gapfill', ...
      'site_inventory.json');
   if ~isfile(inventory_file)
      return
   end
   inventory = jsondecode(fileread(inventory_file));
   for k = 1:numel(inventory)
      if string(inventory(k).site) ~= site
         continue
      end
      neighbors = inventory(k).nearest_donors;
      keep = false(numel(neighbors), 1);
      for n = 1:numel(neighbors)
         % Apply the policy geometry gate at preselection so only loadable
         % candidates are read; the plan re-applies the same gate.
         keep(n) = ismember(string(neighbors(n).family), donor_families) && ...
            neighbors(n).km <= opts.max_donor_distance_km && ...
            abs(neighbors(n).delev_m) <= opts.max_donor_elev_diff_m;
      end
      neighbors = neighbors(keep);
      donor_sites = reshape(string({neighbors.donor}), 1, []);
      return
   end
end

function hourly = hourlyStateMean(series)
   %HOURLYSTATEMEAN Aggregate a sub-hourly donor to hourly state means.
   % Non-numeric or flag channels (aws_type) drop out of the donor role;
   % the mean of the remaining state/flux-density channels is the policy's
   % per-variable-class rule for this donor family.
   names = string(series.Properties.VariableNames);
   metadata = series.Properties.UserData;
   if ~isstruct(metadata)
      metadata = struct();
   end
   if ismember("aws_type", names)
      metadata.aws_types = unique(series.aws_type(isfinite(series.aws_type))).';
   end
   numeric = varfun(@isnumeric, series, 'OutputFormat', 'uniform');
   donor_channels = setdiff(names(numeric), "aws_type", 'stable');
   hourly = retime(series(:, cellstr(donor_channels)), 'hourly', 'mean');
   hourly.Properties.UserData = metadata;
end
