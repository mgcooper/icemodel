function report = buildGapFillReport(kwargs)
   %BUILDGAPFILLREPORT Build the gap-fill before/after Quarto report inputs.
   %
   %  report = icemodel.verification.report.buildGapFillReport()
   %  report = icemodel.verification.report.buildGapFillReport( ...
   %     sites="kanm", render=false)
   %
   % Role
   %  The `.5` report generator (DesignSpec report contract): renders from
   %  SAVED artifacts only — the filled met products, their native
   %  counterparts, the per-site plan summaries, readiness ledgers, and
   %  audit sidecars — never rerunning a model or the engine. It writes
   %  each station's appendix figure set — one 8x1 full-period overview
   %  (POLICY A14/D-19) into the dedicated overview/ subfolder, and the
   %  a bounded method-diverse set of windowed before/after detail figures
   %  into the sibling detail/ subfolder (POLICY A14/D-31), and a small
   %  cohort-level scientific-interpretation set into interpretation/ —
   %  all beneath the
   %  data/preview/figures/gapfill/ namespace (outside the accepted
   %  seasonal/firn ledgers), a figure ledger that must reconcile exactly,
   %  compact summary CSVs, and two QMDs. The main report's fixed
   %  structure is Executive
   %  Summary, Background, Methods, Results, Summary, Appendices. The
   %  appendix embeds only station overview figures. A companion detail
   %  report embeds the method/gap figures and their per-station tables;
   %  both reports have clickable station links. Results carries the
   %  cohort verdict table,
   %  the provenance-derived fill volume per method family, the
   %  admission/held-out skill aggregate, and the stations without any
   %  filled product, all from saved artifacts. Results also publishes an
   %  exhaustive residual-gap table, every non-ready station-year with its
   %  reason, and independently verified native/proxy spans for stations
   %  without products; unexplained product absence blocks publication.
   %  The Results interpretation
   %  catalog states each declared category's reproducible selector,
   %  status, mechanism, policy basis, and figure (or explicit absence).
   %  Detail figures plot ONLY
   %  the filled period plus context: each side's pad equals the filled
   %  period, floored at min_context_days and capped at max_context_days
   %  — never the full record; the station overview is the one deliberate
   %  full-period exception. Each detail panel accents ONLY its own
   %  method; fills by any other method in the context window render
   %  muted grey, keyed by the per-sample provenance registry and the
   %  plan audit (POLICY D-31).
   %
   % Name-value
   %  sites : site tokens to include (default "all" = every site with a
   %     filled product on disk).
   %  max_segments_per_method : representative gap panels per method in
   %     each station/channel/method figure (default 6, spanning the
   %     observed gap-duration range).
   %  max_detail_figures_per_site : maximum station detail figures
   %     (default 6). Selection interleaves duration ranks within each
   %     channel x method-family group, then takes the longest candidates.
   %  min_context_days : context floor per side (default 2), so every
   %     figure shows at least two diurnal cycles of native record on
   %     each side of the fill.
   %  max_context_days : context cap per side (default 366), so a
   %     year-plus fill still shows one full annual cycle of native
   %     record without unbounded windows; between the bounds each side
   %     gets the filled period itself, keeping the fill and its context
   %     visually comparable.
   %  render : run Quarto PDF and HTML rendering after generation
   %     (default true).
   %
   % Returns
   %  report : struct — qmd_file, pdf_file, html_file, detail_qmd_file,
   %     detail_pdf_file, and detail_html_file (rendered paths are ""
   %     when render=false),
   %     figure_ledger (table), summary (per site x channel table),
   %     method_diagnostics (per candidate stratum table),
   %     interpretation_catalog (cohort scientific evidence table), and
   %     fill_by_family (per-station percent of channel-samples per
   %     provenance family), residual_gaps, readiness_blockers, and
   %     absent_products (the exhaustive gap/readiness explanation tables).
   %
   % See also: icemodel.forcing.reconstruct.fillPromiceStation

   arguments
      kwargs.sites (1, :) string = "all"
      kwargs.max_segments_per_method (1, 1) double ...
         {mustBePositive, mustBeInteger} = 6
      kwargs.max_detail_figures_per_site (1, 1) double ...
         {mustBePositive, mustBeInteger} = 6
      kwargs.min_context_days (1, 1) double {mustBePositive} = 2
      kwargs.max_context_days (1, 1) double {mustBePositive} = 366
      kwargs.render (1, 1) logical = true
      kwargs.filled_dir (1, 1) string = ""
      kwargs.qa_dir (1, 1) string = ""
      kwargs.fig_dir (1, 1) string = ""
      kwargs.report_dir (1, 1) string = ""
   end

   % Directory overrides keep the generator testable against fixture trees;
   % production calls use the canonical preview locations.
   repo = icemodel.internal.fullpath;
   filled_dir = defaultDir(kwargs.filled_dir, repo, ...
      fullfile('data', 'input', 'met', 'promice_filled'));
   qa_dir = defaultDir(kwargs.qa_dir, repo, ...
      fullfile('data', 'preview', 'qa', 'gapfill'));
   fig_dir = defaultDir(kwargs.fig_dir, repo, ...
      fullfile('data', 'preview', 'figures', 'gapfill'));
   report_dir = defaultDir(kwargs.report_dir, repo, ...
      fullfile('data', 'preview', 'report'));

   sites = kwargs.sites;
   full_cohort = isscalar(sites) && sites == "all";
   if full_cohort
      hits = dir(fullfile(filled_dir, 'met_*_promice_filled_*_15m.mat'));
      sites = unique(extractBetween(string({hits.name}), "met_", ...
         "_promice_filled"));
   end
   if isempty(sites)
      error('icemodel:report:buildGapFillReport:noFilledProducts', ...
         'no filled met products under %s', filled_dir);
   end
   icemodel.forcing.reconstruct.mustBeStationToken(sites);

   % Validate every consumed artifact against the producer's byte manifest
   % before creating or changing any report output.
   inputs = cell(numel(sites), 1);
   for s = 1:numel(sites)
      inputs{s} = verifySiteInputs(sites(s), qa_dir, filled_dir);
   end
   report_options = pinnedReconstructionOptions(inputs);
   [staging_root, stage_qa_dir, stage_fig_dir, stage_report_dir] = ...
      stageReportDirectories(qa_dir, fig_dir, report_dir);
   staging_cleanup = onCleanup(@() removeStagingDir(staging_root));

   % Build every report artifact in an isolated tree with the same relative
   % layout as production, so Quarto validates the exact links before publish.
   figure_rows = cell(numel(sites), 1);
   summary_rows = cell(numel(sites), 1);
   diagnostic_rows = cell(numel(sites), 1);
   family_rows = cell(numel(sites), 1);
   for s = 1:numel(sites)
      [figure_rows{s}, summary_rows{s}, diagnostic_rows{s}, ...
         family_rows{s}] = siteFigures(sites(s), inputs{s}, ...
         stage_fig_dir, ...
         days([kwargs.min_context_days, kwargs.max_context_days]), ...
         kwargs.max_segments_per_method, ...
         kwargs.max_detail_figures_per_site);
   end
   % Sites can legitimately contribute zero example figures; keep the
   % ledger a schema-stable table so writetable and the QMD always work.
   kept_sites = figure_rows(~cellfun(@isempty, figure_rows));
   figure_ledger = vertcat(kept_sites{:});
   if isempty(figure_ledger)
      figure_ledger = table('Size', [0 10], 'VariableTypes', {'string', ...
         'string', 'string', 'datetime', 'datetime', 'double', ...
         'double', 'logical', 'logical', 'string'}, 'VariableNames', {'site', ...
         'channel', 'method', 'gap_start', 'gap_end', 'duration_hours', ...
         'segments_shown', 'has_before_context', 'has_after_context', ...
         'file'});
   end
   summary = vertcat(summary_rows{:});
   method_diagnostics = vertcat(diagnostic_rows{:});
   fill_by_family = vertcat(family_rows{:});
   [interpretation_catalog, interpretation_rows] = ...
      scientificInterpretations(sites, inputs, summary, ...
      method_diagnostics, stage_fig_dir);
   if ~isempty(interpretation_rows)
      figure_ledger = [figure_ledger; interpretation_rows];
   end

   % Validate the staged ledger before replacing any in-scope PNG. Once the
   % staged set is exact, replace only this report cohort and leave unrelated
   % subset figures untouched. Ledger file values carry their overview/ or
   % detail/ subfolder so reconciliation covers the split layout exactly.
   if numel(unique(figure_ledger.file)) ~= height(figure_ledger)
      error('icemodel:report:buildGapFillReport:duplicateLedgerFiles', ...
         'figure ledger contains duplicate filenames');
   end
   staged_names = stagedFigureNames(stage_fig_dir);
   missing = setdiff(figure_ledger.file, staged_names);
   extra = setdiff(staged_names, figure_ledger.file);
   if ~isempty(missing) || ~isempty(extra)
      error('icemodel:report:buildGapFillReport:ledgerMismatch', ...
         'staged figure ledger mismatch: %d missing, %d unledgered', ...
         numel(missing), numel(extra));
   end
   % Keep all tabular outputs staged with the figures until document
   % generation and optional rendering have both succeeded.
   writetable(figure_ledger, fullfile(stage_qa_dir, ...
      'gapfill_figure_ledger.csv'));
   writetable(summary, fullfile(stage_qa_dir, 'gapfill_summary.csv'));
   writetable(method_diagnostics, fullfile(stage_qa_dir, ...
      'gapfill_method_diagnostics.csv'));
   writetable(fill_by_family, fullfile(stage_qa_dir, ...
      'gapfill_fill_by_family.csv'));
   writetable(interpretation_catalog, fullfile(stage_qa_dir, ...
      'gapfill_interpretation_catalog.csv'));

   % Assemble the combined readiness ledger the report and acceptance use.
   ledgers = cell(numel(inputs), 1);
   for k = 1:numel(inputs)
      % The delimiter is pinned because semicolon-joined reason strings
      % can outnumber commas on heavily-invalid station-years, flipping
      % readtable's delimiter auto-detection (LYN_T/THU_L2 in production).
      ledgers{k} = readtable(inputs{k}.readiness, 'Delimiter', ',', ...
         'TextType', 'string');
   end
   readiness = vertcat(ledgers{:});
   readiness = readiness(ismember(string(readiness.site), sites), :);
   % Policy view over the absolute verdicts uses each producer-pinned
   % acceptance window. Live proxy-directory changes cannot alter a report
   % while all recorded report inputs remain byte-identical.
   readiness = applyAcceptanceWindow(readiness, sites, inputs);
   writetable(readiness, fullfile(stage_qa_dir, ...
      'promice_filled_readiness.csv'));

   % Explain every residual forcing gap from the reconciled producer audit,
   % and refuse a report whose ledger does not exactly match the shipped
   % missing runs. Finite-but-invalid readiness blockers remain a separate
   % table because they are not gaps and need their own exact reason.
   residual_gaps = residualGapTable(sites, inputs);
   readiness_blockers = readinessBlockerTable(readiness);
   writetable(residual_gaps, fullfile(stage_qa_dir, ...
      'gapfill_residual_gaps.csv'));
   writetable(readiness_blockers, fullfile(stage_qa_dir, ...
      'gapfill_readiness_blockers.csv'));

   % Stations with staged native met but no filled product anywhere on
   % disk are a key cohort result. Independently verify their staged proxy
   % support; filename absence alone is never allowed to invent a cause.
   absent_products = sitesWithoutProducts( ...
      inputs, filled_dir, report_options, sites, full_cohort);
   writetable(absent_products, fullfile(stage_qa_dir, ...
      'gapfill_absent_products.csv'));

   stage_qmd = writeQmd(stage_report_dir, stage_qa_dir, stage_fig_dir, ...
      fill_by_family, readiness, method_diagnostics, interpretation_catalog, ...
      sites, residual_gaps, readiness_blockers, absent_products, ...
      report_options, [kwargs.min_context_days, kwargs.max_context_days], ...
      kwargs.max_segments_per_method, ...
      kwargs.max_detail_figures_per_site);
   stage_detail_qmd = writeDetailQmd(stage_report_dir, stage_qa_dir, ...
      stage_fig_dir, summary, method_diagnostics, sites, ...
      kwargs.max_detail_figures_per_site);
   stage_pdf = "";
   stage_detail_pdf = "";
   if kwargs.render
      [stage_pdf, ~] = renderQmd(stage_qmd);
      [stage_detail_pdf, ~] = renderQmd(stage_detail_qmd);
   end
   % Refused stations have no newly staged figures, but their prior report
   % figures are stale evidence and belong in the same rollback-safe cleanup.
   figure_sites = union(sites, string(absent_products.site), 'stable');
   installReportArtifacts(stage_qa_dir, stage_fig_dir, stage_report_dir, ...
      qa_dir, fig_dir, report_dir, figure_sites);

   qmd_file = fullfile(report_dir, 'gapfill-report.qmd');
   detail_qmd_file = fullfile(report_dir, 'gapfill-detail-report.qmd');
   pdf_file = "";
   html_file = "";
   detail_pdf_file = "";
   detail_html_file = "";
   if stage_pdf ~= ""
      pdf_file = fullfile(report_dir, 'gapfill-report.pdf');
      html_file = fullfile(report_dir, 'gapfill-report.html');
   end
   if stage_detail_pdf ~= ""
      detail_pdf_file = fullfile(report_dir, ...
         'gapfill-detail-report.pdf');
      detail_html_file = fullfile(report_dir, ...
         'gapfill-detail-report.html');
   end

   report = struct('qmd_file', string(qmd_file), ...
      'pdf_file', string(pdf_file), ...
      'html_file', string(html_file), ...
      'detail_qmd_file', string(detail_qmd_file), ...
      'detail_pdf_file', string(detail_pdf_file), ...
      'detail_html_file', string(detail_html_file), ...
      'figure_ledger', figure_ledger, 'summary', summary, ...
      'method_diagnostics', method_diagnostics, ...
      'fill_by_family', fill_by_family, ...
      'interpretation_catalog', interpretation_catalog, ...
      'residual_gaps', residual_gaps, ...
      'readiness_blockers', readiness_blockers, ...
      'absent_products', absent_products);
end

function subs = figureSubfolders()
   %FIGURESUBFOLDERS Canonical figure subfolder names for the split layout.
   % POLICY D-31: overview and detail figures live in separate folders and
   % report sections. Scientific Results figures use a third folder so they
   % do not inflate either station appendix set. One named source keeps
   % staging, export, ledger paths, transaction, and QMD links coherent.
   subs = struct('overview', "overview", 'detail', "detail", ...
      'interpretation', "interpretation");
end

function staged_names = stagedFigureNames(stage_fig_dir)
   %STAGEDFIGURENAMES List staged PNGs as subfolder-relative ledger paths.
   subs = figureSubfolders();
   folders = string(struct2cell(subs)).';
   name_cells = cell(numel(folders), 1);
   for f = 1:numel(folders)
      hits = dir(fullfile(stage_fig_dir, char(folders(f)), '*.png'));
      % Forward slashes match the ledger file column on every platform.
      name_cells{f} = folders(f) + "/" + string({hits.name}).';
   end
   staged_names = vertcat(name_cells{:});
end

function [staging_root, stage_qa, stage_fig, stage_report] = ...
      stageReportDirectories(qa_dir, fig_dir, report_dir)
   %STAGEREPORTDIRECTORIES Mirror the output layout in an isolated tree.
   destinations = string([qa_dir, fig_dir, report_dir]);
   output_root = commonDirectory(destinations);
   if output_root == string(filesep)
      error('icemodel:report:buildGapFillReport:outputLayout', ...
         'report outputs must share a writable parent below the filesystem root');
   end
   staging_root = string(tempname(fileparts(char(output_root))));
   mkdir(staging_root);
   stage_qa = fullfile(staging_root, relativeOutput(output_root, qa_dir));
   stage_fig = fullfile(staging_root, relativeOutput(output_root, fig_dir));
   stage_report = fullfile(staging_root, ...
      relativeOutput(output_root, report_dir));
   mkdir(stage_qa);
   mkdir(stage_fig);
   mkdir(stage_report);
   % The split figure layout (POLICY D-31) exists from the start of the
   % staged tree so exports never create ad hoc directories.
   subs = figureSubfolders();
   for folder = reshape(string(struct2cell(subs)), 1, [])
      mkdir(fullfile(stage_fig, char(folder)));
   end
end

function root = commonDirectory(paths)
   %COMMONDIRECTORY Find the shared directory boundary of output paths.
   root = paths(1);
   while ~all(paths == root | startsWith(paths, root + string(filesep)))
      parent = string(fileparts(char(root)));
      if parent == root || parent == ""
         root = string(filesep);
         return
      end
      root = parent;
   end
end

function relative = relativeOutput(root, pathname)
   %RELATIVEOUTPUT Return a path below the shared output root.
   if pathname == root
      relative = "";
   else
      relative = extractAfter(pathname, root + string(filesep));
   end
end

function installReportArtifacts(stage_qa, stage_fig, stage_report, ...
      qa_dir, fig_dir, report_dir, sites)
   %INSTALLREPORTARTIFACTS Publish the complete report with rollback.
   icemodel.helpers.ensureDirExists(qa_dir);
   icemodel.helpers.ensureDirExists(fig_dir);
   icemodel.helpers.ensureDirExists(report_dir);
   [sources, destinations] = artifactPairs(stage_qa, qa_dir);
   % Figures install per split subfolder (POLICY D-31) so the transaction
   % replaces individual PNGs and never moves a whole subfolder over
   % unrelated stations' figures.
   subs = figureSubfolders();
   folders = string(struct2cell(subs)).';
   fig_sources = cell(numel(folders), 1);
   fig_destinations = cell(numel(folders), 1);
   for f = 1:numel(folders)
      sub = char(folders(f));
      icemodel.helpers.ensureDirExists(fullfile(fig_dir, sub));
      [fig_sources{f}, fig_destinations{f}] = artifactPairs( ...
         fullfile(stage_fig, sub), fullfile(fig_dir, sub));
   end
   sources = [sources; vertcat(fig_sources{:})];
   destinations = [destinations; vertcat(fig_destinations{:})];
   [src, dst] = artifactPairs(stage_report, report_dir);
   sources = [sources; src];
   destinations = [destinations; dst];
   staged_report_names = string({dir(stage_report).name});
   managed_report = ["gapfill-report.pdf"; "gapfill-report.html"; ...
      "gapfill-report_files"; "gapfill-detail-report.pdf"; ...
      "gapfill-detail-report.html"; "gapfill-detail-report_files"];
   stale_report = setdiff(managed_report, staged_report_names, 'stable');
   destinations = [destinations; ...
      fullfile(string(report_dir), stale_report)];

   % Stale selected-site figures are part of the same transaction even
   % though they have no staged replacement. Each subfolder reconciles
   % against its own staged set; any selected-site PNG still sitting flat
   % in fig_dir comes from the pre-split layout and is stale by definition.
   obsolete_cells = cell(numel(folders) + 1, 1);
   for f = 1:numel(folders)
      sub = char(folders(f));
      prior_names = pngNames(fullfile(fig_dir, sub));
      staged_sub = pngNames(fullfile(stage_fig, sub));
      selected = selectedSiteMask(prior_names, sites);
      % reshape pins the column orientation: an all-replaced subfolder
      % makes setdiff/fullfile return row-shaped empties that would
      % break the destination concatenation.
      obsolete_cells{f} = reshape(fullfile(string(fig_dir), sub, ...
         setdiff(prior_names(selected), staged_sub, 'stable')), [], 1);
   end
   flat_names = pngNames(fig_dir);
   obsolete_cells{end} = reshape(fullfile(string(fig_dir), ...
      flat_names(selectedSiteMask(flat_names, sites))), [], 1);
   destinations = [destinations; vertcat(obsolete_cells{:})];

   % Move every current destination aside before installing anything. On
   % any failure, remove only newly installed paths and restore every backup.
   output_root = commonDirectory(string([qa_dir, fig_dir, report_dir]));
   backup_root = string(tempname(fileparts(char(output_root))));
   mkdir(backup_root);
   backed = false(numel(destinations), 1);
   backups = strings(numel(destinations), 1);
   installed = false(numel(sources), 1);
   try
      for k = 1:numel(destinations)
         if isfile(destinations(k)) || isfolder(destinations(k))
            backups(k) = fullfile(backup_root, sprintf('%06d', k));
            [ok, message] = movefile(destinations(k), backups(k), 'f');
            if ~ok
               error('icemodel:report:buildGapFillReport:backupFailed', ...
                  'could not back up %s: %s', destinations(k), message);
            end
            backed(k) = true;
         end
      end
      for k = 1:numel(sources)
         [ok, message] = movefile(sources(k), destinations(k), 'f');
         if ~ok
            error('icemodel:report:buildGapFillReport:installFailed', ...
               'could not install %s: %s', destinations(k), message);
         end
         installed(k) = true;
      end
   catch err
      for k = numel(installed):-1:1
         if installed(k)
            removePath(destinations(k));
         end
      end
      restored = true;
      restore_message = "";
      for k = numel(backed):-1:1
         if backed(k)
            [ok, message] = movefile(backups(k), destinations(k), 'f');
            restored = restored && ok;
            if ~ok
               restore_message = restore_message + newline + message;
            end
         end
      end
      if restored
         removeStagingDir(backup_root);
         rethrow(err)
      end
      error('icemodel:report:buildGapFillReport:rollbackFailed', ...
         'report publish failed and rollback is incomplete under %s:%s', ...
         backup_root, restore_message);
   end
   removeStagingDir(backup_root);
end

function names = pngNames(dirpath)
   %PNGNAMES List PNG basenames in one directory as a string column.
   hits = dir(fullfile(dirpath, '*.png'));
   names = string({hits.name}).';
end

function selected = selectedSiteMask(names, sites)
   %SELECTEDSITEMASK Mark figure basenames belonging to selected sites.
   selected = false(size(names));
   for site = reshape(sites, 1, [])
      selected = selected | startsWith(names, site + "_");
   end
end

function [sources, destinations] = artifactPairs(stage_dir, final_dir)
   %ARTIFACTPAIRS Pair immediate staged files/directories with destinations.
   entries = dir(stage_dir);
   entries = entries(~ismember(string({entries.name}), [".", ".."]));
   sources = strings(numel(entries), 1);
   destinations = strings(numel(entries), 1);
   for k = 1:numel(entries)
      sources(k) = fullfile(entries(k).folder, entries(k).name);
      destinations(k) = fullfile(final_dir, entries(k).name);
   end
end

function removePath(pathname)
   %REMOVEPATH Remove one explicitly enumerated transaction artifact.
   if isfolder(pathname)
      rmdir(pathname, 's');
   elseif isfile(pathname)
      delete(pathname);
   end
end

function removeStagingDir(pathname)
   %REMOVESTAGINGDIR Remove an isolated report staging directory.
   if isfolder(pathname)
      rmdir(pathname, 's');
   end
end

function readiness = applyAcceptanceWindow(readiness, sites, inputs)
   %APPLYACCEPTANCEWINDOW Add producer-pinned policy columns to the ledger.
   % One window lookup per site; a year overlaps the window when any part
   % of the calendar year falls inside it. Absolute verdicts pass
   % through; only not_forcing_ready years wholly outside the window
   % become out_of_policy_window. The policy view grades run readiness,
   % so it derives from the icemodel verdict; the snowmodel verdict
   % (POLICY A5) rides along untouched in its own column.
   window_start = NaT(height(readiness), 1, 'TimeZone', 'UTC');
   window_end = NaT(height(readiness), 1, 'TimeZone', 'UTC');
   policy = string(readiness.verdict_icemodel);
   for s = 1:numel(sites)
      in_site = readiness.site == sites(s);
      window = inputs{s}.acceptance_window;
      window_start(in_site) = window(1);
      window_end(in_site) = window(2);
      year_start = datetime(readiness.year(in_site), 1, 1, ...
         'TimeZone', 'UTC');
      year_end_exclusive = datetime(readiness.year(in_site) + 1, 1, 1, ...
         'TimeZone', 'UTC');
      site_policy = policy(in_site);
      if ~isnat(window(1))
         in_window = year_end_exclusive > window(1) ...
            & year_start <= window(2);
         site_policy(~in_window ...
            & site_policy == "not_forcing_ready") = ...
            "out_of_policy_window";
      end
      policy(in_site) = site_policy;
   end
   readiness.window_start = window_start;
   readiness.window_end = window_end;
   readiness.policy_verdict = policy;
end

function [rows, summary_rows, diagnostics, family_row] = siteFigures( ...
      site, inputs, fig_dir, context_bounds, max_segments_per_method, ...
      max_detail_figures_per_site)
   %SITEFIGURES Render the overview plus a bounded grouped-detail set.
   % The overview leads (POLICY A14/D-19) and lands in the overview/
   % subfolder. Candidate detail figures remain one channel and method
   % apiece (POLICY D-31), but POLICY A14 bounds the per-station set and
   % requires method diversity rather than one figure for every candidate.
   filled = loadMet(inputs.filled, site, "promice_filled");
   native = loadMet(inputs.native, site, "promice");
   A = load(inputs.plan);
   audit = A.audit_record;
   diagnostics = methodDiagnostics(site, A.plan_record);
   % The per-family fill volume derives from the shipped provenance
   % channels, the registry the engine stamps, so the Results table can
   % never disagree with the product (POLICY D-31 Results contract).
   family_row = siteFillFamilies(site, filled);

   % Composition refusal rows (method 'unfilled') explain residual gaps;
   % they are not fills, so example figures and fill totals must exclude
   % them — the summary instead counts them per channel. Darkness
   % zero-fill rows ARE fills but summarize a record-spanning set of
   % night samples, so plotting one as an "example gap" would violate
   % the report's windowed-figure contract; they stay in the totals and
   % out of the figure pool.
   method_names = string(audit.method);
   is_fill = ~ismember(method_names, ["unfilled", "native_context"]);
   % Record-spanning summary rows (darkness zero-fill, the winter-albedo
   % seasonal bridge) are fills but not windowable example gaps.
   is_example = is_fill & ~ismember(method_names, ...
      ["darkness_zero", "winter_albedo_bridge"]);

   % Every shipped science channel has a paired provenance column. Derive
   % the exhaustive summary scope from that contract, not from audit rows:
   % a fully native or wholly unresolved channel may have no fill segment.
   channels = provenanceChannels(filled);
   if isempty(channels)
      channels = unique(string(audit.channel), 'stable');
   end
   summary_cells = cell(numel(channels), 1);
   candidate_template = struct('channel', "", 'method', "", ...
      'duration_hours', NaN, 'segments', table(), ...
      'method_audit', table());
   candidates = candidate_template([], 1);
   n_candidates = 0;
   for c = 1:numel(channels)
      channel = channels(c);
      in_channel = string(audit.channel) == channel;
      n_unfilled = countMissingRuns(~isfinite(filled.(channel)));
      fills = audit(in_channel & is_fill, :);
      segment_vars = {'start_time', 'end_time', 'duration_hours', 'method'};
      seg = sortrows(audit(in_channel & is_example, segment_vars), ...
         'duration_hours', 'descend');
      seg = representativeSegments(seg, max_segments_per_method);
      methods = unique(string(seg.method), 'stable');
      % Different methods generally describe different gaps. Preserve one
      % method per candidate figure so a panel set never implies a same-gap
      % comparison; curation happens only after all station candidates exist.
      for m = 1:numel(methods)
         n_candidates = n_candidates + 1;
         method_seg = seg(string(seg.method) == methods(m), :);
         candidates(n_candidates, 1).channel = channel;
         candidates(n_candidates, 1).method = methods(m);
         candidates(n_candidates, 1).duration_hours = ...
            max(method_seg.duration_hours);
         candidates(n_candidates, 1).segments = method_seg;
         % ALL audited fill segments of this method — not only displayed
         % representatives — attribute samples to the accent layer.
         candidates(n_candidates, 1).method_audit = ...
            fills(string(fills.method) == methods(m), :);
      end
      summary_cells{c} = table(site, channel, height(fills), ...
         sum(fills.duration_hours), 0, n_unfilled, ...
         100 * mean(~isfinite(filled.(channel))), ...
         'VariableNames', {'site', 'channel', 'segments_filled', ...
         'hours_filled', 'figures', 'unfilled_segments', ...
         'residual_missing_pct'});
   end

   % Interleave duration ranks within channel x method-family groups, then
   % cut at the station budget. This keeps manual review and render time
   % bounded without mixing methods inside a figure.
   selected = curateSiteDetails(candidates, max_detail_figures_per_site);
   method_rows = cell(numel(selected), 1);
   for k = 1:numel(selected)
      candidate = candidates(selected(k));
      method_rows{k} = gapFigure(site, candidate.channel, ...
         candidate.segments, candidate.method_audit, native, filled, ...
         fig_dir, context_bounds);
   end
   rows = table();
   if ~isempty(method_rows)
      rows = vertcat(method_rows{:});
   end
   selected_channels = strings(0, 1);
   if ~isempty(selected)
      selected_channels = string({candidates(selected).channel}).';
   end
   for c = 1:numel(channels)
      summary_cells{c}.figures = nnz(selected_channels == channels(c));
   end

   % The station overview leads the appendix (POLICY A14/D-19), so its
   % ledger row goes first for this site.
   overview_row = overviewFigure(site, filled, fig_dir);
   if isempty(rows)
      rows = overview_row;
   else
      rows = [overview_row; rows];
   end
   summary_rows = vertcat(summary_cells{:});
end

function channels = provenanceChannels(filled)
   %PROVENANCECHANNELS Science channels shipping a paired provenance column.
   % One derivation shared by the figure scope and the fill-by-family
   % summary so the two can never diverge.
   names = string(filled.Properties.VariableNames);
   provenance_names = names(endsWith(names, "_provenance"));
   channels = erase(provenance_names, "_provenance");
   channels = channels(ismember(channels, names));
end

function families = fillFamilyNames()
   %FILLFAMILYNAMES Reconstruction family names from the provenance registry.
   % The provenance code registry is the single source of family names
   % (STYLE.local SSOT rule); observed-equivalent codes (observed plus the
   % raw/clamped shortwave measurements, POLICY A7) and the missing
   % sentinel are report layers, not fill families.
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   names = reshape(string(fieldnames(codes)), 1, []);
   families = names(~ismember(names, ["observed", "raw_shortwave", ...
      "clamped_shortwave", "missing"]));
end

function row = siteFillFamilies(site, filled)
   %SITEFILLFAMILIES Percent of channel-samples per provenance family.
   % The denominator is every sample of every provenance-paired science
   % channel, so observed, each family, and missing sum to 100 within
   % rounding; the counts come from the stamped codes alone, never from
   % re-derived audit spans.
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   families = fillFamilyNames();
   channels = provenanceChannels(filled);
   n_total = max(1, numel(channels) * height(filled));
   counts = zeros(1, numel(families));
   n_observed = 0;
   n_missing = 0;
   for c = 1:numel(channels)
      provenance = double(filled.(channels(c) + "_provenance"));
      n_observed = n_observed + nnz(ismember(provenance, double([ ...
         codes.observed, codes.raw_shortwave, codes.clamped_shortwave])));
      n_missing = n_missing + nnz(provenance == double(codes.missing));
      for f = 1:numel(families)
         counts(f) = counts(f) ...
            + nnz(provenance == double(codes.(families(f))));
      end
   end
   values = 100 * [n_observed, counts, n_missing] / n_total;
   row = array2table(round(values, 3), 'VariableNames', ...
      cellstr(["observed_pct", families + "_pct", "missing_pct"]));
   row = addvars(row, string(site), 'Before', 1, ...
      'NewVariableNames', 'site');
end

function selected = representativeSegments(seg, max_per_method)
   %REPRESENTATIVESEGMENTS Select bounded duration-spanning panels per method.
   selected = seg([], :);
   methods = unique(string(seg.method), 'stable');
   selected_cells = cell(numel(methods), 1);
   for m = 1:numel(methods)
      candidates = sortrows(seg(string(seg.method) == methods(m), :), ...
         'duration_hours', 'descend');
      n_pick = min(max_per_method, height(candidates));
      % Evenly spaced duration ranks expose scale sensitivity while keeping
      % manual review bounded independently of the number of gaps.
      pick = unique(round(linspace(1, height(candidates), n_pick)), ...
         'stable');
      selected_cells{m} = candidates(pick, :);
   end
   if ~isempty(selected_cells)
      selected = vertcat(selected_cells{:});
   end
end

function selected = curateSiteDetails(candidates, max_figures)
   %CURATESITEDETAILS Bound one station's channel/method review surface.
   if isempty(candidates)
      selected = zeros(0, 1);
      return
   end
   channel = string({candidates.channel}).';
   method = string({candidates.method}).';
   family = extractBefore(method + ":", ":");
   duration_hours = [candidates.duration_hours].';
   index = (1:numel(candidates)).';
   [groups, ~, group_of] = unique(table(channel, family), 'rows');
   ranked = cell(height(groups), 1);
   for g = 1:height(groups)
      in_group = table(index(group_of == g), ...
         duration_hours(group_of == g), channel(group_of == g), ...
         method(group_of == g), 'VariableNames', ...
         {'index', 'duration_hours', 'channel', 'method'});
      in_group = sortrows(in_group, ...
         {'duration_hours', 'channel', 'method'}, ...
         {'descend', 'ascend', 'ascend'});
      in_group.rank = (1:height(in_group)).';
      ranked{g} = in_group;
   end
   pooled = sortrows(vertcat(ranked{:}), ...
      {'rank', 'duration_hours', 'channel', 'method'}, ...
      {'ascend', 'descend', 'ascend', 'ascend'});
   selected = pooled.index(1:min(max_figures, height(pooled)));
end

function n_runs = countMissingRuns(missing)
   %COUNTMISSINGRUNS Count contiguous residual runs from final product truth.
   edges = diff([false; missing(:); false]);
   n_runs = nnz(edges == 1);
end

function diagnostics = methodDiagnostics(site, plan)
   %METHODDIAGNOSTICS Export selection/evaluation evidence per stratum.
    metric_fields = ["n", "coverage", "bias", "rmse", "correlation", ...
       "variability_ratio", "within_gap_observed_spread", ...
       "bound_violations", ...
      "boundary_jump_rate", "sigma1_coverage", "sigma2_coverage", ...
      "provenance_accounting"];
    variable_names = ["site", "channel", "candidate", "season", ...
       "bucket", "decision", "max_validated_hours", ...
       "selection_baseline", "selection_baseline_rmse", ...
       "uncertainty_status", ...
      "selection_" + metric_fields, "evaluation_" + metric_fields, ...
      "denial_reason"];
   rows = cell(0, 1);
   n_rows = 0;
   for c = 1:numel(plan.channels)
      channel = string(plan.channels(c).channel);
      methods = plan.channels(c).methods;
      denials = plan.channels(c).denials;
      for m = 1:numel(methods)
         season = string(methods(m).seasons);
         for b = 1:numel(methods(m).buckets)
            bucket = methods(m).buckets(b);
             selection = metricStratum(methods(m).selection, ...
                bucket, season);
             evaluation = evaluationStratum(methods(m).evaluation, ...
                bucket, season);
             baseline = "";
             baseline_rmse = NaN;
             if ismember('baseline', selection.Properties.VariableNames)
                baseline = string(selection.baseline);
                baseline_rmse = selection.baseline_rmse;
             end
             n_rows = n_rows + 1;
             rows{n_rows, 1} = [{char(site), char(channel), ...
                char(methods(m).name), char(season), bucket, ...
                'admitted', methods(m).max_validated_hours(b), ...
                char(baseline), baseline_rmse, ...
                char(methods(m).uncertainty)}, ...
               metricCells(selection, metric_fields), ...
               metricCells(evaluation, metric_fields), {''}];
         end
      end
      for d = 1:height(denials)
         season = "all";
         if ismember('season', denials.Properties.VariableNames)
            season = string(denials.season{d});
         end
          n_rows = n_rows + 1;
          rows{n_rows, 1} = [{char(site), char(channel), ...
             denials.candidate{d}, char(season), denials.bucket(d), ...
             'denied', NaN, '', NaN, ''}, ...
            num2cell(nan(1, 2 * numel(metric_fields))), ...
            denials.reasons(d)];
      end
   end
   if n_rows == 0
       variable_types = [repmat("string", 1, 4), "double", "string", ...
          "double", "string", "double", "string", ...
         repmat("double", 1, 2 * numel(metric_fields)), ...
         "string"];
      diagnostics = table('Size', [0 numel(variable_names)], ...
         'VariableTypes', cellstr(variable_types), ...
         'VariableNames', cellstr(variable_names));
   else
      diagnostics = cell2table(vertcat(rows{:}), 'VariableNames', ...
         cellstr(variable_names));
   end
end

function values = metricCells(metrics, fields)
   %METRICCELLS Return schema-ordered values or explicit NaNs.
   values = num2cell(nan(1, numel(fields)));
   if isempty(metrics)
      return
   end
   for k = 1:numel(fields)
      values{k} = metrics.(fields(k));
   end
end

function metrics = metricStratum(table_rows, bucket, season)
   %METRICSTRATUM Select one bucket-season metric row.
   metrics = table();
   if isempty(table_rows)
      return
   end
   keep = table_rows.bucket == bucket & table_rows.season == season;
   if any(keep)
      metrics = table_rows(keep, :);
   end
end

function metrics = evaluationStratum(grades, bucket, season)
   %EVALUATIONSTRATUM Select one evaluation row from bucket grades.
   metrics = table();
   if isempty(grades)
      return
   end
   for k = 1:numel(grades)
      if grades(k).bucket == bucket
         metrics = metricStratum(grades(k).by_stratum, bucket, season);
         return
      end
   end
end

function row = gapFigure(site, channel, seg, method_audit, native, ...
      filled, fig_dir, context_bounds)
   %GAPFIGURE Render grouped windowed panels and one figure-ledger row.
   n_segments = height(seg);
   n_cols = min(2, n_segments);
   n_rows = ceil(n_segments / n_cols);
   method_names = unique(string(seg.method), 'stable');
   if numel(method_names) ~= 1
      error('icemodel:report:buildGapFillReport:mixedMethodFigure', ...
         'a grouped figure must contain exactly one method');
   end
   method = method_names(1);
   native = observedNative(native, filled, channel);
   has_before = false(n_segments, 1);
   has_after = false(n_segments, 1);

   % POLICY D-31 accent contract: the filled series splits into the
   % panel's own method (accent), other methods' fills (muted context),
   % and observed samples, keyed by the shipped provenance channel plus
   % the audit spans of this method. The provenance column is a product
   % contract; a detail figure without it must fail loudly.
   provenance_name = channel + "_provenance";
   if ~ismember(provenance_name, string(filled.Properties.VariableNames))
      error('icemodel:report:buildGapFillReport:missingProvenance', ...
         'detail figure requires %s in the filled product', ...
         provenance_name);
   end
   own_mask = methodSampleMask(filled.Properties.RowTimes, method_audit);
   layers = icemodel.verification.report.methodFillLayers( ...
      filled.(channel), double(filled.(provenance_name)), own_mask);
   % Layer views subset the filled product so variable metadata (units)
   % rides along and the shared plotting layer converts all three series
   % identically.
   own_view = filled(:, channel);
   own_view.(channel) = layers.own_fill;
   other_view = filled(:, channel);
   other_view.(channel) = layers.other_fill;

   % One bounded multi-panel figure per station/channel/method makes the
   % manual review surface independent of the total number of filled gaps.
   fig = icemodel.plot.newFigure(width=1200, ...
      height=max(360, 280 * n_rows));
   cleaner = onCleanup(@() close(fig));
   layout = tiledlayout(fig, n_rows, n_cols, ...
      'TileSpacing', 'compact', 'Padding', 'compact');
   for k = 1:n_segments
      gap_start = seg.start_time(k);
      % Audit end_time is the final interval-start posting. Plot spans use
      % audited duration after delivery support returns to quarter-hourly.
      gap_end = gap_start + hours(seg.duration_hours(k));
      gap_len = gap_end - gap_start;
      pad = min(max(gap_len, context_bounds(1)), context_bounds(2));
      [window_start, window_end, has_before(k), has_after(k)] = ...
         observedContextWindow(native, channel, gap_start, gap_end, pad, ...
         context_bounds);

      % Each panel retains the original gap-scale context contract while
      % multiple duration regimes for one method share one reviewable PNG.
      ax = nexttile(layout);
      hold(ax, 'on');
      out = icemodel.plot.compareTimeseries({own_view, other_view, ...
         native}, channel, axes=ax, ...
         names=[method, "other methods", "promice"], ...
         startdate=window_start, enddate=window_end);
      % compareTimeseries keys colors by source name; re-key the three
      % layers to the accent contract (own method accent, foreign fills
      % muted grey, observed dark). Markers keep a single-posting fill
      % visible between its NaN-masked neighbors; the plotted map keeps
      % the layer-to-handle pairing exact when a layer is empty in this
      % window and was skipped.
      layer_colors = {layers.own_color, layers.other_color, ...
         layers.observed_color};
      plotted = find(out.plotted);
      for j = 1:numel(plotted)
         set(out.lines(j), 'Color', layer_colors{plotted(j)});
         if plotted(j) < 3
            set(out.lines(j), 'Marker', '.', 'MarkerSize', 6);
         end
      end
      icemodel.plot.markTimeSpan(ax, gap_start, gap_end);
      title(ax, icemodel.plot.formatDuration(seg.duration_hours(k)), ...
         'Interpreter', 'none', 'FontSize', 10, ...
         'FontWeight', 'normal', 'Color', [0.2 0.2 0.2]);
      grid(ax, 'on');
   end
   sgtitle(layout, sprintf('%s %s — %s (%d representative gaps)', ...
      upper(site), channel, method, n_segments), ...
      'Interpreter', 'none', 'FontSize', 12, ...
      'FontWeight', 'normal', 'Color', [0.2 0.2 0.2]);

   % Stable station/channel/method names ensure one file per grouped
   % figure; detail figures land in the split detail/ subfolder and the
   % ledger records the subfolder-relative path (POLICY D-31).
   subs = figureSubfolders();
   method_token = matlab.lang.makeValidName(method);
   filename = sprintf('%s_%s_%s_gapfill.png', site, channel, method_token);
   exportgraphics(fig, fullfile(fig_dir, char(subs.detail), filename), ...
      'Resolution', 130);
   clear cleaner
   [~, primary] = max(seg.duration_hours);
   primary_start = seg.start_time(primary);
   primary_end = primary_start + hours(seg.duration_hours(primary));
   row = table(site, channel, method, primary_start, ...
      primary_end, seg.duration_hours(primary), n_segments, ...
      all(has_before), all(has_after), subs.detail + "/" + filename, ...
      'VariableNames', {'site', 'channel', 'method', 'gap_start', ...
      'gap_end', 'duration_hours', 'segments_shown', ...
      'has_before_context', 'has_after_context', 'file'});
end

function mask = methodSampleMask(times, method_audit)
   %METHODSAMPLEMASK Mark samples inside a method's audited fill segments.
   % Audit end_time is the final interval-start POSTING on the hourly
   % reconstruction axis, but the delivered product repeats each posting
   % over its quarter-hour support — closed-bounds containment would
   % strand the last posting's +15/+30/+45 rows in the "other methods"
   % layer (review pass 9). Span containment mirrors gapFigure's
   % duration-based window instead.
   mask = false(numel(times), 1);
   for k = 1:height(method_audit)
      mask = mask | (times >= method_audit.start_time(k) ...
         & times < method_audit.start_time(k) ...
         + hours(method_audit.duration_hours(k)));
   end
end

function row = overviewFigure(site, filled, fig_dir)
   %OVERVIEWFIGURE Render the eight-channel full-period station overview.
   % POLICY A14/D-19: each station's appendix leads with ONE full-period
   % overview — the science channels stacked 8x1 — so a reviewer sees
   % every fill in the context of the whole record before the windowed
   % detail figures. Observed samples draw dark; filled samples overlay
   % in the shared accent color. Hourly-or-coarser decimation keeps
   % multi-decade 15-minute records renderable without changing the
   % full-period visual story; detail figures retain local resolution.
   % The D-19 eight-channel layout: the A5 seven-channel SSOT plus the
   % derived swu (B10 — planned nowhere, shipped everywhere), filtered to
   % the channels the product actually carries with provenance.
   channels = [icemodel.forcing.reconstruct.icemodelRequiredChannels(), ...
      "swu"];
   names = string(filled.Properties.VariableNames);
   channels = channels(ismember(channels, names) ...
      & ismember(channels + "_provenance", names));
   times = filled.Properties.RowTimes;
   style = icemodel.verification.report.gapfillFigureStyle();
   % The point budget prevents exportgraphics from silently returning a
   % white canvas for long records with eight densely populated axes.
   hourly_step = round(hours(1) / median(diff(times)));
   budget_step = ceil(height(filled) / style.max_overview_points);
   step = max([1, hourly_step, budget_step]);
   pick = (1:step:height(filled)).';
   codes = icemodel.forcing.reconstruct.provenanceCodes();

   fig = icemodel.plot.newFigure(width=1400, ...
      height=max(900, 170 * numel(channels)));
   cleaner = onCleanup(@() close(fig));
   layout = tiledlayout(fig, numel(channels), 1, ...
      'TileSpacing', 'compact', 'Padding', 'compact');
   for c = 1:numel(channels)
      channel = channels(c);
      values = filled.(channel)(pick);
      provenance = filled.(channel + "_provenance")(pick);
      % Raw-fallback and clamped shortwave (codes 13/14) are builder-
      % selected raw MEASUREMENTS (A7), so they render as observed here
      % exactly as observedNative treats them in the detail figures.
      is_observed = ismember(provenance, [codes.observed, ...
         codes.raw_shortwave, codes.clamped_shortwave]);
      observed = values;
      observed(~is_observed) = NaN;
      reconstructed = values;
      reconstructed(is_observed) = NaN;
      ax = nexttile(layout);
      hold(ax, 'on');
      % Fills draw first so the observed record reads as the foreground
      % truth wherever both exist at plot resolution.
      plot(ax, times(pick), reconstructed, '-', ...
         'Color', style.accent, 'LineWidth', 0.4);
      plot(ax, times(pick), observed, '-', ...
         'Color', style.observed, 'LineWidth', 0.4);
      ylabel(ax, channel, 'Interpreter', 'none');
      grid(ax, 'on');
      if c < numel(channels)
         ax.XTickLabel = {};
      end
   end
   sgtitle(layout, sprintf( ...
      '%s — full-period overview (observed dark, filled orange)', ...
      upper(site)), 'Interpreter', 'none', 'FontSize', 12, ...
      'FontWeight', 'normal', 'Color', [0.2 0.2 0.2]);

   % Overviews land in their split overview/ subfolder and the ledger
   % records the subfolder-relative path (POLICY D-31).
   subs = figureSubfolders();
   filename = sprintf('%s_overview.png', site);
   % Force the large multi-axis canvas to finish painting before export.
   % Without this synchronization, long real records can produce a valid
   % white PNG while short fixture records happen to render in time.
   drawnow
   exportgraphics(fig, fullfile(fig_dir, char(subs.overview), filename), ...
      'Resolution', 130);
   clear cleaner
   % The ledger row mirrors the gapFigure schema so one table holds every
   % appendix figure; the "gap" span is the whole product period.
   span_hours = hours(times(end) - times(1));
   row = table(string(site), "all", "overview", times(1), times(end), ...
      span_hours, numel(channels), true, true, ...
      subs.overview + "/" + filename, ...
      'VariableNames', {'site', 'channel', 'method', 'gap_start', ...
      'gap_end', 'duration_hours', 'segments_shown', ...
      'has_before_context', 'has_after_context', 'file'});
end

function native = observedNative(native, filled, channel)
   %OBSERVEDNATIVE Mask source values the filled-product ledger did not observe.
   provenance_name = channel + "_provenance";
   if ~ismember(provenance_name, ...
         string(filled.Properties.VariableNames))
      return
   end
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   observed = false(height(native), 1);
   [tf, loc] = ismember(native.Properties.RowTimes, ...
      filled.Properties.RowTimes);
    source_code = filled.(provenance_name)(loc(tf));
    observed(tf) = source_code == codes.observed ...
       | source_code == codes.raw_shortwave ...
       | source_code == codes.clamped_shortwave;
   native.(channel)(~observed) = NaN;
end

function [window_start, window_end, has_before, has_after] = ...
      observedContextWindow(native, channel, gap_start, gap_end, pad, ...
      context_bounds)
   %OBSERVEDCONTEXTWINDOW Extend a plot to the nearest observed side context.
   times = native.Properties.RowTimes;
   finite = isfinite(native.(channel));
   before = find(finite & times < gap_start ...
      & times >= gap_start - context_bounds(2), 1, 'last');
   after = find(finite & times >= gap_end ...
      & times <= gap_end + context_bounds(2), 1, 'first');
   has_before = ~isempty(before);
   has_after = ~isempty(after);
   window_start = gap_start - pad;
   window_end = gap_end + pad;
   if has_before
      window_start = min(window_start, ...
         times(before) - context_bounds(1));
      window_start = max(window_start, gap_start - context_bounds(2));
   end
   if has_after
      window_end = max(window_end, times(after) + context_bounds(1));
      window_end = min(window_end, gap_end + context_bounds(2));
   end
end

function absent = sitesWithoutProducts(inputs, filled_dir, opts, ...
      selected_sites, full_cohort)
   %SITESWITHOUTPRODUCTS Verify staged stations lacking any filled product.
   % A missing filename proves only that no product was published. For each
   % absent station, validate the staged proxy inventory with the same
   % acceptance-window function as production and record the exact native
   % and proxy spans. An empty validated proxy inventory is an explicit
   % noProxyWindow refusal; any real overlap is an unexplained missing
   % product and blocks the report.
   native_dirs = strings(numel(inputs), 1);
   for k = 1:numel(inputs)
      native_dirs(k) = string(fileparts(inputs{k}.native));
   end
   native_dirs = unique(native_dirs);
   native_cells = cell(numel(native_dirs), 1);
   for k = 1:numel(native_dirs)
      hits = dir(fullfile(native_dirs(k), 'met_*_promice_*_15m.mat'));
      names = string({hits.name}).';
      % The native glob would also match filled artifacts if any were
      % ever co-located; exclude them so only true native stations count.
      names = names(~contains(names, "_promice_filled_"));
      native_cells{k} = extractBetween(names, "met_", "_promice_");
   end
   native_sites = unique(vertcat(native_cells{:}));
   filled_hits = dir(fullfile(filled_dir, ...
      'met_*_promice_filled_*_15m.mat'));
   filled_sites = unique(extractBetween(string({filled_hits.name}).', ...
      "met_", "_promice_filled"));
   absent_sites = setdiff(native_sites, filled_sites);
   % An explicit subset report owns only its requested stations. Other
   % staged native stations may be awaiting production and must neither
   % block nor be mislabeled by this narrower report. The default "all"
   % cohort retains the exhaustive repository-wide absence gate.
   if ~full_cohort
      absent_sites = intersect(absent_sites, selected_sites);
   end
   names = {'site', 'native_start', 'native_end', 'proxy_start', ...
      'proxy_end', 'reason', 'explanation'};
   if isempty(absent_sites)
      absent = table('Size', [0 numel(names)], 'VariableTypes', ...
         {'string', 'datetime', 'datetime', 'datetime', 'datetime', ...
         'string', 'string'}, 'VariableNames', names);
      return
   end

   % Derived tokens feed report text; validate them like every other public
   % station token so stray files fail loudly, never silently.
   icemodel.forcing.reconstruct.mustBeStationToken(absent_sites);
   rows = cell(numel(absent_sites), 1);
   for k = 1:numel(absent_sites)
      site = absent_sites(k);
      hits = [];
      native_dir = "";
      for d = reshape(native_dirs, 1, [])
         candidate = dir(fullfile(d, ...
            sprintf('met_%s_promice_*_15m.mat', site)));
         if ~isempty(candidate)
            if ~isempty(hits)
               error( ...
                  'icemodel:report:buildGapFillReport:ambiguousNativeProduct', ...
                  'native met for absent station %s appears in two roots', ...
                  site);
            end
            hits = candidate;
            native_dir = d;
         end
      end
      [native, native_file] = ...
         icemodel.forcing.reconstruct.loadWidestTimetable(hits);
      if isempty(native)
         error('icemodel:report:buildGapFillReport:missingNativeProduct', ...
            'cannot verify absent product because native met is empty: %s', ...
            native_file);
      end
      native_times = native.Properties.RowTimes;
      location = reportArtifactLocation(native.Properties.UserData, site);
      window = icemodel.forcing.reconstruct.acceptanceWindow( ...
         site, met_dir=native_dir, location=location, opts=opts);
      native_start = native_times(1);
      native_end = native_times(end);
      if any(isnat(window))
         reason = "noProxyWindow";
         explanation = ...
            "No validated staged 15-minute MAR/MERRA proxy exists, so " + ...
            "POLICY A6 defines no fillable product period.";
      else
         overlaps = any(native_times >= window(1) ...
            & native_times <= window(2));
         if overlaps
            error( ...
               ['icemodel:report:buildGapFillReport:' ...
               'unexplainedMissingProduct'], ...
               ['%s has native samples inside its validated proxy window ' ...
               'but no filled product'], site);
         end
         reason = "windowRecordDisjoint";
         explanation = ...
            "The native record and staged proxy window do not overlap.";
      end
      rows{k} = {site, native_start, native_end, window(1), window(2), ...
         reason, explanation};
   end
   absent = cell2table(vertcat(rows{:}), 'VariableNames', names);
end

function location = reportArtifactLocation(metadata, site)
   %REPORTARTIFACTLOCATION Recover the saved point needed for proxy checks.
   % PROMICE products use top-level lat/lon/elev; other canonical builders
   % may use WGS84 field names or a site_location record.
   if isfield(metadata, 'lat')
      location = struct('lat_wgs84', metadata.lat, ...
         'lon_wgs84', metadata.lon, 'elev_m', metadata.elev);
   elseif all(isfield(metadata, ["lat_wgs84", "lon_wgs84", "elev_m"]))
      location = struct('lat_wgs84', metadata.lat_wgs84, ...
         'lon_wgs84', metadata.lon_wgs84, 'elev_m', metadata.elev_m);
   elseif isfield(metadata, 'site_location') ...
         && all(isfield(metadata.site_location, ...
         ["lat_wgs84", "lon_wgs84", "elev_m"]))
      location = metadata.site_location;
   else
      error('icemodel:report:buildGapFillReport:missingNativeLocation', ...
         'native met for %s lacks a saved location identity', site);
   end
end

function residual = residualGapTable(sites, inputs)
   %RESIDUALGAPTABLE Reconcile every shipped forcing gap to its audit row.
   % The report may summarize reasons only after the producer's final
   % unfilled ledger exactly matches contiguous missing runs in the product.
   names = {'site', 'year', 'channel', 'start_time', 'end_time', ...
      'duration_hours', 'reason', 'context_id'};
   rows = cell(numel(sites), 1);
   required = icemodel.forcing.reconstruct.icemodelRequiredChannels();
   for s = 1:numel(sites)
      site = sites(s);
      filled = loadMet(inputs{s}.filled, site, "promice_filled");
      saved = load(inputs{s}.plan, 'audit_record');
      audit = saved.audit_record;
      product_names = string(filled.Properties.VariableNames);
      is_unfilled = string(audit.method) == "unfilled";
      % The seven IceModel channels are always in scope. Dedicated
      % precipitation processing can leave optional phase columns missing;
      % add only channels the producer itself reconciled as unfilled (for
      % example canonical ppt), while the separate readiness table explains
      % the ppt-OR-snowf snow-model contract.
      tracked = intersect(unique([required, ...
         string(audit.channel(is_unfilled)).'], 'stable'), ...
         product_names, 'stable');
      is_residual = is_unfilled ...
         & ismember(string(audit.channel), tracked);
      ledger = audit(is_residual, :);

      % Derive the expected contiguous runs from the shipped values, then
      % compare the segment identity rather than trusting either source
      % independently.
      actual_cells = cell(numel(tracked), 1);
      for c = 1:numel(tracked)
         channel = tracked(c);
         missing = ~isfinite(filled.(channel));
         segments = icemodel.forcing.reconstruct.auditSegments( ...
            filled.Properties.RowTimes, missing, channel, ...
            "unfilled", "");
         if ~isempty(segments)
            actual_cells{c} = cell2table(vertcat(segments{:}), ...
               'VariableNames', audit.Properties.VariableNames);
         end
      end
      actual_cells = actual_cells(~cellfun(@isempty, actual_cells));
      actual = audit([], :);
      if ~isempty(actual_cells)
         actual = vertcat(actual_cells{:});
      end
      if ~sameGapCoverage(actual, ledger, ...
            filled.Properties.RowTimes, tracked)
         product_span = gapTableSpan(actual);
         ledger_span = gapTableSpan(ledger);
         error('icemodel:report:buildGapFillReport:residualGapMismatch', ...
               ['shipped missing runs and reconciled unfilled audit rows ' ...
            'disagree for %s (product segments=%d, span=%s; ' ...
            'ledger segments=%d, span=%s)'], site, height(actual), ...
            product_span, height(ledger), ledger_span);
      end
      if isempty(ledger)
         rows{s} = table();
         continue
      end
      actual = sortrows(actual, {'channel', 'start_time'});
      ledger = sortrows(ledger, {'channel', 'start_time'});

      % Keep the producer's final-tier denial text verbatim as the reason;
      % publish the shipped 15-minute endpoints rather than the producer's
      % hourly posting endpoints. A shipped contiguous run can span several
      % producer strata, so retain every distinct overlapping denial.
      n = height(actual);
      reason = strings(n, 1);
      context_id = strings(n, 1);
      for r = 1:n
         run_end = actual.start_time(r) ...
            + hours(actual.duration_hours(r));
         ledger_end = ledger.start_time + hours(ledger.duration_hours);
         overlap = string(ledger.channel) == string(actual.channel(r)) ...
            & ledger.start_time < run_end ...
            & ledger_end > actual.start_time(r);
         reasons = unique(string(ledger.detail(overlap)), 'stable');
         reasons = reasons(strlength(reasons) > 0);
         if isempty(reasons)
            error( ...
               'icemodel:report:buildGapFillReport:residualGapMismatch', ...
               'missing final-tier denial for %s %s at %s', ...
               site, string(actual.channel(r)), ...
               string(actual.start_time(r)));
         end
         reason(r) = strjoin(reasons, " | ");
         contexts = unique(string(ledger.context_id(overlap)), 'stable');
         contexts = contexts(strlength(contexts) > 0);
         context_id(r) = strjoin(contexts, " | ");
      end
      rows{s} = table(repmat(site, n, 1), ...
         year(actual.start_time), string(actual.channel), ...
         actual.start_time, actual.end_time, actual.duration_hours, ...
         reason, context_id, ...
         'VariableNames', names);
   end
   present = ~cellfun(@isempty, rows);
   if any(present)
      residual = sortrows(vertcat(rows{present}), ...
         {'site', 'year', 'channel', 'start_time'});
   else
      residual = table('Size', [0 numel(names)], 'VariableTypes', ...
         {'string', 'double', 'string', 'datetime', 'datetime', ...
         'double', 'string', 'string'}, 'VariableNames', names);
   end
end

function span = gapTableSpan(gaps)
   %GAPTABLESPAN Compact diagnostic span for a failed coverage reconciliation.
   span = "<empty>";
   if ~isempty(gaps)
      first = min(gaps.start_time);
      exclusive_end = gaps.start_time + hours(gaps.duration_hours);
      span = string(first) + ".." + string(max(exclusive_end));
   end
end

function tf = sameGapCoverage(actual, ledger, times, channels)
   %SAMEGAPCOVERAGE Compare audit and product missingness on shipped support.
   % Producer rows can split one continuous gap at season/method strata and
   % use hourly posting endpoints, while the artifact ships quarter-hour
   % support. Exclusive intervals derived from start + duration are the
   % cadence-independent identity.
   tf = true;
   if isempty(times)
      tf = isempty(actual) && isempty(ledger);
      return
   end
   dt = median(diff(times));
   support_end = times(end) + dt;
   for source = {actual, ledger}
      rows = source{1};
      if any(rows.start_time < times(1) ...
            | rows.start_time + hours(rows.duration_hours) > support_end)
         tf = false;
         return
      end
   end
   for channel = reshape(channels, 1, [])
      actual_mask = false(numel(times), 1);
      ledger_mask = false(numel(times), 1);
      actual_rows = actual(string(actual.channel) == channel, :);
      ledger_rows = ledger(string(ledger.channel) == channel, :);
      for r = 1:height(actual_rows)
         actual_mask = actual_mask ...
            | (times >= actual_rows.start_time(r) ...
            & times < actual_rows.start_time(r) ...
            + hours(actual_rows.duration_hours(r)));
      end
      for r = 1:height(ledger_rows)
         ledger_mask = ledger_mask ...
            | (times >= ledger_rows.start_time(r) ...
            & times < ledger_rows.start_time(r) ...
            + hours(ledger_rows.duration_hours(r)));
      end
      if ~isequal(actual_mask, ledger_mask)
         tf = false;
         return
      end
   end
end

function blockers = readinessBlockerTable(readiness)
   %READINESSBLOCKERTABLE Explain every non-ready use without aggregation.
   keep = string(readiness.policy_verdict) ~= "ready" ...
      | string(readiness.verdict_snowmodel) ~= "ready";
   names = {'site', 'year', 'policy_verdict', 'reason_icemodel', ...
      'verdict_snowmodel', 'reason_snowmodel', 'window_start', ...
      'window_end', 'plain_explanation'};
   if ~any(keep)
      blockers = table('Size', [0 numel(names)], 'VariableTypes', ...
         {'string', 'double', 'string', 'string', 'string', 'string', ...
         'datetime', 'datetime', 'string'}, 'VariableNames', names);
      return
   end
   source = readiness(keep, :);
   explanation = strings(height(source), 1);
   for k = 1:height(source)
      ice_explanation = "";
      if string(source.policy_verdict(k)) == "out_of_policy_window"
         ice_explanation = ...
            "The calendar year lies outside the staged proxy window.";
      elseif string(source.verdict_icemodel(k)) ~= "ready"
         ice_explanation = "IceModel: " ...
            + string(source.reason_icemodel(k)) + ".";
      end
      snow_explanation = "";
      if string(source.verdict_snowmodel(k)) ~= "ready"
         snow_explanation = "Snow model: " ...
            + string(source.reason_snowmodel(k)) + ".";
      end
      explanation(k) = strtrim(ice_explanation + " " + snow_explanation);
   end
   blockers = table(string(source.site), source.year, ...
      string(source.policy_verdict), string(source.reason_icemodel), ...
      string(source.verdict_snowmodel), string(source.reason_snowmodel), ...
      source.window_start, source.window_end, explanation, ...
      'VariableNames', names);
end

function verdicts = verdictSummary(readiness)
   %VERDICTSUMMARY Count station-years and stations per policy verdict.
   % The counts derive from the same policy_verdict column the executive
   % summary sentence uses, so the table and the sentence can never
   % disagree.
   names = {'policy_verdict', 'station_years', 'stations'};
   if isempty(readiness)
      verdicts = table('Size', [0 3], 'VariableTypes', ...
         {'string', 'double', 'double'}, 'VariableNames', names);
      return
   end
   [groups, verdict] = findgroups(string(readiness.policy_verdict));
   rows = cell(numel(verdict), 1);
   for g = 1:numel(verdict)
      in_group = groups == g;
      rows{g} = table(verdict(g), nnz(in_group), ...
         numel(unique(string(readiness.site(in_group)))), ...
         'VariableNames', names);
   end
   verdicts = vertcat(rows{:});
end

function admission = admissionSummary(method_diagnostics)
   %ADMISSIONSUMMARY Aggregate plan diagnostics into a cohort skill table.
   % Compact per-channel, per-family view of the admission results
   % (POLICY D-31 Results contract). Held-out improvement is
   % 100 * (1 - selection_rmse / selection_baseline_rmse) — the same
   % baseline-relative criterion the admission gate enforces. The family
   % token is the audit/registry label prefix (the part before ':' in
   % names such as "donor:aws10"; plain labels are their own family), so
   % no alias mapping is restated here.
   names = {'channel', 'family', 'admitted_strata', ...
      'denied_candidates', 'median_improvement_pct'};
   if isempty(method_diagnostics)
      admission = table('Size', [0 5], 'VariableTypes', {'string', ...
         'string', 'double', 'double', 'double'}, 'VariableNames', names);
      return
   end
   channel = string(method_diagnostics.channel);
   family = extractBefore(string(method_diagnostics.candidate) + ":", ":");
   decision = string(method_diagnostics.decision);
   improvement = 100 * (1 - method_diagnostics.selection_rmse ...
      ./ method_diagnostics.selection_baseline_rmse);
   [groups, key_channel, key_family] = findgroups(channel, family);
   rows = cell(numel(key_channel), 1);
   for g = 1:numel(key_channel)
      in_group = groups == g;
      admitted = in_group & decision == "admitted";
      rows{g} = table(key_channel(g), key_family(g), nnz(admitted), ...
         nnz(in_group & decision == "denied"), ...
         round(median(improvement(admitted), 'omitnan'), 1), ...
         'VariableNames', names);
   end
   admission = sortrows(vertcat(rows{:}), {'channel', 'family'});
end

function [catalog, rows] = scientificInterpretations(sites, inputs, ...
      summary, diagnostics, fig_dir)
   %SCIENTIFICINTERPRETATIONS Select and render the declared evidence set.
   % The category list is fixed; selectors operate on producer-pinned
   % artifacts and absent categories remain explicit catalog rows.
   cases = [ ...
      newInterpretationCase("swd_fill_geometry", ...
      "expected_or_accepted", ...
      "Within KANL, choose the station-year with the most finite " + ...
      "positive reconstructed SWD samples.", ...
      "Filled values occur only at formerly missing timestamps, so " + ...
      "their plotted amplitude is conditional on gap timing and solar " + ...
      "geometry rather than a second full-record climatology.", ...
      "A7, A15, D-28, D-29"); ...
      newInterpretationCase("long_gap_limited_context", ...
      "expected_or_accepted", ...
      "Choose the longest audited fill among required forcing channels, " + ...
      "excluding record-spanning darkness and winter summary rows.", ...
      "Long outages can leave little or no native context inside a " + ...
      "one-year diagnostic window; the selected accepted method and " + ...
      "audit duration remain explicit.", ...
      "A6, A11, B1, B8"); ...
      newInterpretationCase("climatology_variability", ...
      "diagnostic_limitation", ...
      "Among admitted climatology strata used by a fill, choose the " + ...
      "finite selection variability ratio farthest from one.", ...
      "Day-of-year/hour climatology preserves seasonal and diurnal " + ...
      "structure but can compress event-scale variability.", ...
      "B1, B7, B8, D-29"); ...
      newInterpretationCase("boundary_behavior", ...
      "diagnostic_limitation", ...
      "Among admitted methods used by a fill, choose the largest " + ...
      "positive selection/evaluation boundary-jump rate.", ...
      "Seam blending limits discontinuities where reconstruction meets " + ...
      "native data; the jump metric remains diagnostic evidence.", ...
      "A15, B6, B8"); ...
      newInterpretationCase("hourly_disaggregation", ...
      "expected_or_accepted", ...
      "Choose the longest audited MAR or MERRA-2 proxy fill among " + ...
      "required forcing channels.", ...
      "Filled hourly proxy postings are expanded to 15 minutes with a " + ...
      "smooth mean-preserving curve; observed postings remain exact " + ...
      "held copies.", ...
      "A7, D-30"); ...
      newInterpretationCase("residual_missing", ...
      "unresolved", ...
      "Choose the station/channel with the largest final residual-missing " + ...
      "percentage, then its longest contiguous missing run.", ...
      "A value rejected by every valid tier remains missing and is " + ...
      "ledgered rather than fabricated.", ...
      "A5, A11, A15"); ...
      newInterpretationCase("native_sensor_flat_run", ...
      "diagnostic_limitation", ...
      "Choose the longest producer-recorded corroborated native " + ...
      "sensor-burial/rime run in the selected cohort.", ...
      "The staged native artifact remains unchanged, but implicated " + ...
      "channels are excluded from training and donor truth and are " + ...
      "reconstructed with explicit provenance.", ...
      "A1, A7, D-37")];

   audit = combinedAudit(sites, inputs);
   required = icemodel.forcing.reconstruct.icemodelRequiredChannels();
   audit = audit(ismember(string(audit.channel), required), :);
   methods = string(audit.method);
   excluded = ismember(methods, ["unfilled", "native_context", ...
      "darkness_zero", "winter_albedo_bridge"]);

   % KANL is the motivating SWD visual question, but its window is selected
   % by the same declared annual rule every time rather than hand-picked.
   kanl = find(lower(sites) == "kanl", 1);
   if ~isempty(kanl)
      cases(1).present = true;
      cases(1).site = sites(kanl);
      cases(1).channel = "swd";
      cases(1).method = "mixed_reconstruction";
   end

   eligible = find(~excluded);
   if ~isempty(eligible)
      [~, local] = max(audit.duration_hours(eligible));
      cases(2) = assignAuditCase(cases(2), audit(eligible(local), :), ...
         audit.duration_hours(eligible(local)));
   end

   % Climatology selection is based on held-out variability evidence, then
   % joined to an actually used audit segment before it can enter Results.
   admitted = string(diagnostics.decision) == "admitted";
   ratios = diagnostics.selection_variability_ratio;
   climate = admitted & string(diagnostics.candidate) == "climatology" ...
      & isfinite(ratios) & ratios > 0;
   climate_rows = find(climate);
   [~, order] = sort(abs(log(ratios(climate_rows))), 'descend');
   for i = reshape(climate_rows(order), 1, [])
      match = matchingAudit(audit, diagnostics(i, :));
      if ~isempty(match)
         [~, longest] = max(match.duration_hours);
         cases(3) = assignAuditCase(cases(3), match(longest, :), ...
            ratios(i));
         break
      end
   end

   % A positive diagnostic rate is materially present; an all-zero cohort
   % remains an explicit absent row rather than a contrived example.
   boundary = max([diagnostics.selection_boundary_jump_rate, ...
      diagnostics.evaluation_boundary_jump_rate], [], 2, 'omitnan');
   boundary_rows = find(admitted & isfinite(boundary) & boundary > 0);
   [~, order] = sort(boundary(boundary_rows), 'descend');
   for i = reshape(boundary_rows(order), 1, [])
      match = matchingAudit(audit, diagnostics(i, :));
      if ~isempty(match)
         [~, longest] = max(match.duration_hours);
         cases(4) = assignAuditCase(cases(4), match(longest, :), ...
            boundary(i));
         break
      end
   end

   proxy = find(~excluded & (startsWith(methods, "proxy:mar") ...
      | startsWith(methods, "proxy:merra2")));
   if ~isempty(proxy)
      [~, local] = max(audit.duration_hours(proxy));
      cases(5) = assignAuditCase(cases(5), audit(proxy(local), :), ...
         audit.duration_hours(proxy(local)));
   end

   science_summary = summary(ismember(string(summary.channel), required), :);
   if ~isempty(science_summary)
      [largest_missing, i] = max(science_summary.residual_missing_pct, [], ...
         'omitnan');
      if isfinite(largest_missing) && largest_missing > 0
         cases(6).present = true;
         cases(6).site = string(science_summary.site(i));
         cases(6).channel = string(science_summary.channel(i));
         cases(6).method = "unfilled";
         cases(6).score = largest_missing;
      end
   end

   % Producer-recorded native QC findings are systematic report evidence,
   % not a hand-picked station note. The longest finding represents the
   % category; an unaffected cohort keeps the fixed absent row.
   flat_runs = combinedFlatRunFindings(sites, inputs);
   if ~isempty(flat_runs)
      [longest, i] = max(flat_runs.n_days);
      channels = split(string(flat_runs.channels(i)), ",");
      cases(7).present = true;
      cases(7).site = string(flat_runs.site(i));
      cases(7).channel = channels(1);
      cases(7).method = "native_sensor_exclusion";
      cases(7).focus_start = flat_runs.start_time(i);
      cases(7).focus_end = flat_runs.end_time(i);
      cases(7).score = longest;
   end

   row_cells = cell(numel(cases), 1);
   for k = 1:numel(cases)
      if cases(k).present
         [cases(k), row_cells{k}] = renderInterpretationCase( ...
            cases(k), sites, inputs, fig_dir);
      else
         cases(k).interpretation = ...
            "No material example is present in the selected cohort.";
      end
   end
   kept = row_cells(~cellfun(@isempty, row_cells));
   rows = emptyFigureLedger();
   if ~isempty(kept)
      rows = vertcat(kept{:});
   end
   catalog = struct2table(cases);
   catalog = removevars(catalog, {'focus_start', 'focus_end'});
end

function value = newInterpretationCase(category, status, selector, ...
      mechanism, policy_basis)
   %NEWINTERPRETATIONCASE One schema-stable scientific evidence category.
   missing_time = NaT(1, 1, 'TimeZone', 'UTC');
   value = struct('category', string(category), 'present', false, ...
      'status', string(status), 'site', "", 'channel', "", ...
      'method', "", 'focus_start', missing_time, ...
      'focus_end', missing_time, 'window_start', missing_time, ...
      'window_end', missing_time, 'score', NaN, ...
      'selection_rule', string(selector), ...
      'mechanism', string(mechanism), ...
      'policy_basis', string(policy_basis), ...
      'interpretation', "", 'figure', "");
end

function value = assignAuditCase(value, audit_row, score)
   %ASSIGNAUDITCASE Attach one audited segment to a category.
   value.present = true;
   value.site = string(audit_row.site);
   value.channel = string(audit_row.channel);
   value.method = string(audit_row.method);
   value.focus_start = audit_row.start_time;
   value.focus_end = audit_row.start_time ...
      + hours(audit_row.duration_hours);
   value.score = score;
end

function audit = combinedAudit(sites, inputs)
   %COMBINEDAUDIT Join producer-pinned station audit tables.
   cells = cell(numel(sites), 1);
   for s = 1:numel(sites)
      saved = load(inputs{s}.plan, 'audit_record');
      cells{s} = addvars(saved.audit_record, ...
         repmat(sites(s), height(saved.audit_record), 1), ...
         'Before', 1, 'NewVariableNames', 'site');
   end
   audit = vertcat(cells{:});
end

function findings = combinedFlatRunFindings(sites, inputs)
   %COMBINEDFLATRUNFINDINGS Join producer-pinned native-QC evidence.
   cells = cell(numel(sites), 1);
   for s = 1:numel(sites)
      saved = load(inputs{s}.plan, 'flat_run_findings_record');
      if ~isfield(saved, 'flat_run_findings_record') ...
            || isempty(saved.flat_run_findings_record)
         cells{s} = table();
         continue
      end
      cells{s} = addvars(saved.flat_run_findings_record, ...
         repmat(sites(s), height(saved.flat_run_findings_record), 1), ...
         'Before', 1, 'NewVariableNames', 'site');
   end
   present = ~cellfun(@isempty, cells);
   findings = table();
   if any(present)
      findings = vertcat(cells{present});
   end
end

function match = matchingAudit(audit, diagnostic)
   %MATCHINGAUDIT Join one admitted candidate to its used audit spans.
   candidate = string(diagnostic.candidate);
   methods = string(audit.method);
   keep = string(audit.site) == string(diagnostic.site) ...
      & string(audit.channel) == string(diagnostic.channel) ...
      & (methods == candidate | startsWith(methods, candidate + ":"));
   match = audit(keep, :);
end

function [value, row] = renderInterpretationCase(value, sites, inputs, ...
      fig_dir)
   %RENDERINTERPRETATIONCASE Render one selected month-to-year diagnostic.
   input_index = find(sites == value.site, 1);
   filled = loadMet(inputs{input_index}.filled, value.site, ...
      "promice_filled");
   native = loadMet(inputs{input_index}.native, value.site, "promice");
   times = filled.Properties.RowTimes;
   cadence = median(diff(times));

   if value.category == "swd_fill_geometry"
      [value, fill_mask, observed_mask] = selectKanlSwdYear(value, filled);
   else
      [fill_mask, observed_mask] = reconstructionMasks(filled, ...
         value.channel);
   end
   if value.category == "residual_missing"
      missing = ~isfinite(filled.(value.channel));
      [value.focus_start, value.focus_end] = longestMissingRun( ...
         times, missing, cadence);
      if isnat(value.focus_start)
         value.present = false;
         value.interpretation = ...
            "No residual missing run remains in the selected product.";
         row = emptyFigureLedger();
         return
      end
   end

   [value.window_start, value.window_end] = interpretationWindow( ...
      times, value.focus_start, value.focus_end, value.category);
   is_sensor = value.category == "native_sensor_flat_run";
   if ~is_sensor
      native = observedNative(native, filled, value.channel);
   end
   reconstructed = filled(:, value.channel);
   reconstructed.(value.channel)(~fill_mask) = NaN;
   observed = native(:, value.channel);
   style = icemodel.verification.report.gapfillFigureStyle();

   is_cadence = value.category == "hourly_disaggregation";
   is_geometry = value.category == "swd_fill_geometry";
   is_residual = value.category == "residual_missing";
   n_tiles = 1 + is_cadence + 3 * is_geometry + is_residual;
   fig = icemodel.plot.newFigure(width=1200, height=340 * n_tiles);
   cleaner = onCleanup(@() close(fig));
   layout = tiledlayout(fig, n_tiles, 1, ...
      'TileSpacing', 'compact', 'Padding', 'compact');
   ax = nexttile(layout);
   plotInterpretationTimeseries(ax, reconstructed, observed, ...
      value.channel, value.window_start, value.window_end, style);
   markVisibleFocus(ax, value.focus_start, value.focus_end, ...
      value.window_start, value.window_end);
   title(ax, "Month-to-year context", 'FontSize', 10, ...
      'FontWeight', 'normal', 'Color', [0.2 0.2 0.2]);

   if is_cadence || is_residual
      defaults = interpretationFigureDefaults();
      center = value.focus_start ...
         + (value.focus_end - value.focus_start) / 2;
      zoom_start = center - days(defaults.cadence_zoom_days / 2);
      zoom_end = center + days(defaults.cadence_zoom_days / 2);
      ax = nexttile(layout);
      plotInterpretationTimeseries(ax, reconstructed, observed, ...
         value.channel, zoom_start, zoom_end, style);
      markVisibleFocus(ax, value.focus_start, value.focus_end, ...
         zoom_start, zoom_end);
      if is_cadence
         detail_title = "Three-day cadence detail";
      else
         detail_title = "Three-day residual detail";
      end
      title(ax, detail_title, 'FontSize', 10, ...
         'FontWeight', 'normal', 'Color', [0.2 0.2 0.2]);
   elseif is_geometry
      [zoom_start, zoom_end] = swdGeometryZoom( ...
         filled, fill_mask, value.window_start, value.window_end);
      ax = nexttile(layout);
      plotInterpretationTimeseries(ax, reconstructed, observed, ...
         value.channel, zoom_start, zoom_end, style);
      title(ax, "30-day highest-fill detail", 'FontSize', 10, ...
         'FontWeight', 'normal', 'Color', [0.2 0.2 0.2]);
      [boundary_start, boundary_end, boundary_time] = ...
         swdBoundaryZoom(filled, fill_mask, observed_mask, ...
         value.window_start, value.window_end);
      ax = nexttile(layout);
      plotInterpretationTimeseries(ax, reconstructed, observed, ...
         value.channel, boundary_start, boundary_end, style);
      xline(ax, boundary_time, ':', 'Color', [0.45 0.45 0.45], ...
         'HandleVisibility', 'off');
      title(ax, "Three-day largest fill/observed transition", ...
         'FontSize', 10, 'FontWeight', 'normal', ...
         'Color', [0.2 0.2 0.2]);
      ax = nexttile(layout);
      plotSwdGeometry(ax, filled, fill_mask, observed_mask, ...
         value.window_start, value.window_end, style);
   end
   sgtitle(layout, sprintf('%s %s — %s', upper(value.site), ...
      value.channel, strrep(value.category, "_", " ")), ...
      'Interpreter', 'none', 'FontSize', 12, ...
      'FontWeight', 'normal', 'Color', [0.2 0.2 0.2]);

   subs = figureSubfolders();
   filename = sprintf('%s_%s_%s.png', value.site, value.channel, ...
      matlab.lang.makeValidName(value.category));
   exportgraphics(fig, fullfile(fig_dir, char(subs.interpretation), ...
      filename), 'Resolution', 130);
   clear cleaner
   value.figure = subs.interpretation + "/" + filename;
   value.interpretation = interpretationText(value);
   row = table(value.site, value.channel, ...
      "interpretation:" + value.category, value.focus_start, ...
      value.focus_end, hours(value.focus_end - value.focus_start), 1, ...
      value.window_start < value.focus_start, ...
      value.window_end > value.focus_end, value.figure, ...
      'VariableNames', {'site', 'channel', 'method', 'gap_start', ...
      'gap_end', 'duration_hours', 'segments_shown', ...
      'has_before_context', 'has_after_context', 'file'});
end

function [zoom_start, zoom_end, boundary_time] = swdBoundaryZoom( ...
      filled, fill_mask, observed_mask, window_start, window_end)
   %SWDBOUNDARYZOOM Select the largest finite fill/observed transition.
   defaults = interpretationFigureDefaults();
   times = filled.Properties.RowTimes;
   values = filled.swd;
   transition = (fill_mask(1:end - 1) & observed_mask(2:end)) ...
      | (observed_mask(1:end - 1) & fill_mask(2:end));
   transition = transition & isfinite(values(1:end - 1)) ...
      & isfinite(values(2:end)) ...
      & times(2:end) >= window_start & times(2:end) < window_end;
   jumps = abs(diff(values));
   jumps(~transition) = -Inf;
   [~, selected] = max(jumps);
   boundary_time = times(selected + 1);
   zoom_start = max(window_start, boundary_time ...
      - days(defaults.cadence_zoom_days / 2));
   zoom_end = min(window_end, boundary_time ...
      + days(defaults.cadence_zoom_days / 2));
end

function [zoom_start, zoom_end] = swdGeometryZoom( ...
      filled, fill_mask, window_start, window_end)
   %SWDGEOMETRYZOOM Select the 30-day span with most positive SWD fills.
   defaults = interpretationFigureDefaults();
   times = filled.Properties.RowTimes;
   day_start = dateshift(window_start, 'start', 'day');
   day_end = dateshift(window_end, 'start', 'day');
   days_in_window = (day_start:days(1):day_end - days(1)).';
   counts = zeros(size(days_in_window));
   positive_fill = fill_mask & filled.swd > 0;
   for k = 1:numel(days_in_window)
      counts(k) = nnz(positive_fill ...
         & times >= days_in_window(k) ...
         & times < days_in_window(k) + days(1));
   end
   rolling = movsum(counts, [0, defaults.swd_zoom_days - 1]);
   [~, selected] = max(rolling);
   zoom_start = days_in_window(selected);
   zoom_end = min(zoom_start + days(defaults.swd_zoom_days), window_end);
end

function markVisibleFocus(ax, focus_start, focus_end, ...
      window_start, window_end)
   %MARKVISIBLEFOCUS Mark only the focus span visible in one plot window.
   visible_start = max(focus_start, window_start);
   visible_end = min(focus_end, window_end);
   if visible_end > visible_start
      icemodel.plot.markTimeSpan(ax, visible_start, visible_end);
   end
end

function [value, fill_mask, observed_mask] = selectKanlSwdYear( ...
      value, filled)
   %SELECTKANLSWDYEAR Apply the declared annual SWD evidence selector.
   [fill_mask, observed_mask] = reconstructionMasks(filled, "swd");
   times = filled.Properties.RowTimes;
   years_present = unique(year(times));
   counts = zeros(size(years_present));
   for k = 1:numel(years_present)
      counts(k) = nnz(year(times) == years_present(k) ...
         & fill_mask & filled.swd > 0);
   end
   [n_fill, selected] = max(counts);
   selected_year = years_present(selected);
   value.focus_start = datetime(selected_year, 1, 1, ...
      'TimeZone', 'UTC');
   value.focus_end = datetime(selected_year + 1, 1, 1, ...
      'TimeZone', 'UTC');
   value.score = n_fill / nnz(year(times) == selected_year);
end

function [fill_mask, observed_mask] = reconstructionMasks(filled, channel)
   %RECONSTRUCTIONMASKS Split provenance into reconstructed and observed.
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   provenance = double(filled.(channel + "_provenance"));
   observed_mask = ismember(provenance, double([codes.observed, ...
      codes.raw_shortwave, codes.clamped_shortwave]));
   fill_mask = isfinite(filled.(channel)) & ~observed_mask ...
      & provenance ~= double(codes.missing);
end

function [run_start, run_end] = longestMissingRun(times, missing, cadence)
   %LONGESTMISSINGRUN Return the longest true run as an end-exclusive span.
   edges = diff([false; missing(:); false]);
   starts = find(edges == 1);
   ends = find(edges == -1) - 1;
   if isempty(starts)
      run_start = NaT(1, 1, 'TimeZone', 'UTC');
      run_end = run_start;
      return
   end
   [~, longest] = max(ends - starts + 1);
   run_start = times(starts(longest));
   run_end = times(ends(longest)) + cadence;
end

function [window_start, window_end] = interpretationWindow(times, ...
      focus_start, focus_end, category)
   %INTERPRETATIONWINDOW Bound diagnostics to readable months or one year.
   defaults = interpretationFigureDefaults();
   if category == "swd_fill_geometry"
      target = focus_end - focus_start;
   elseif category == "boundary_behavior" ...
         || category == "hourly_disaggregation"
      target = days(defaults.min_days);
   else
      target = min(max(focus_end - focus_start, ...
         days(defaults.min_days)), days(defaults.max_days));
   end
   if category == "boundary_behavior"
      center = focus_start;
   else
      center = focus_start + (focus_end - focus_start) / 2;
   end
   window_start = center - target / 2;
   window_end = center + target / 2;
   cadence = median(diff(times));
   record_start = times(1);
   record_end = times(end) + cadence;
   if window_start < record_start
      window_end = min(record_end, ...
         window_end + (record_start - window_start));
      window_start = record_start;
   end
   if window_end > record_end
      window_start = max(record_start, ...
         window_start - (window_end - record_end));
      window_end = record_end;
   end
end

function defaults = interpretationFigureDefaults()
   %INTERPRETATIONFIGUREDEFAULTS Single source of diagnostic view spans.
   defaults = struct('min_days', 90, 'max_days', 366, ...
      'cadence_zoom_days', 3, 'swd_zoom_days', 30);
end

function plotInterpretationTimeseries(ax, reconstructed, observed, ...
      channel, start_time, end_time, style)
   %PLOTINTERPRETATIONTIMESERIES Use the shared comparison plotting layer.
   hold(ax, 'on');

   % Draw the complete reconstructed product as a thin neutral underlay.
   % The observed and filled overlays remain provenance-honest, while the
   % underlay makes the actual connection across their color boundary
   % visible instead of creating an artificial two-line break.
   combined = reconstructed;
   [present, observed_index] = ismember( ...
      combined.Properties.RowTimes, observed.Properties.RowTimes);
   observed_values = NaN(height(combined), 1);
   observed_values(present) = observed.(channel)(observed_index(present));
   use_observed = isfinite(observed_values);
   combined.(channel)(use_observed) = observed_values(use_observed);
   in_window = combined.Properties.RowTimes >= start_time ...
      & combined.Properties.RowTimes <= end_time;
   plot(ax, combined.Properties.RowTimes(in_window), ...
      combined.(channel)(in_window), 'Color', style.context, ...
      'LineWidth', 0.5, 'HandleVisibility', 'off');

   % Provenance-specific overlays retain the report's established colors
   % and legend labels.
   out = icemodel.plot.compareTimeseries({reconstructed, observed}, ...
      channel, axes=ax, names=["reconstructed", "promice"], ...
      startdate=start_time, enddate=end_time);
   colors = {style.accent, style.observed};
   plotted = find(out.plotted);
   for k = 1:numel(plotted)
      set(out.lines(k), 'Color', colors{plotted(k)});
   end
   grid(ax, 'on');
end

function plotSwdGeometry(ax, filled, fill_mask, observed_mask, ...
      start_time, end_time, style)
   %PLOTSWDGEOMETRY Relate the KANL annual fill amplitudes to solar angle.
   metadata = filled.Properties.UserData;
   in_window = filled.Properties.RowTimes >= start_time ...
      & filled.Properties.RowTimes < end_time;
   elevation = icemodel.forcing.helpers.solarElevation( ...
      filled.Properties.RowTimes(in_window), metadata.lat, metadata.lon);
   values = filled.swd(in_window);
   is_fill = fill_mask(in_window) & values > 0;
   is_observed = observed_mask(in_window) & values > 0;
   scatter(ax, elevation(is_observed), values(is_observed), 5, ...
      style.observed, 'filled', 'MarkerFaceAlpha', 0.15, ...
      'DisplayName', 'promice');
   hold(ax, 'on');
   scatter(ax, elevation(is_fill), values(is_fill), 9, ...
      style.accent, 'filled', 'MarkerFaceAlpha', 0.55, ...
      'DisplayName', 'reconstructed');
   xlabel(ax, 'solar elevation [degrees]');
   ylabel(ax, 'swd [W m-2]');
   title(ax, "Positive SWD versus solar elevation", 'FontSize', 10, ...
      'FontWeight', 'normal', 'Color', [0.2 0.2 0.2]);
   legend(ax, 'Location', 'best', 'Color', [1 1 1], ...
      'TextColor', [0.2 0.2 0.2], 'EdgeColor', [0.7 0.7 0.7]);
   grid(ax, 'on');
end

function text = interpretationText(value)
   %INTERPRETATIONTEXT Emit concise data-specific scientific interpretation.
   switch value.category
      case "swd_fill_geometry"
         text = sprintf([ ...
            'Selected year has %.2f%% finite positive reconstructed SWD ' ...
            'samples. Orange values occupy only missing timestamps and ' ...
            'cluster at low solar elevation because KANL outages are ' ...
            'concentrated in shoulder hours; they are not a random ' ...
            'sample of the full observed amplitude distribution. The ' ...
            '30-day and worst-transition panels expose local continuity, ' ...
            'while the solar-elevation panel shows the timestamp-' ...
            'selection geometry.'], ...
            100 * value.score);
      case "long_gap_limited_context"
         text = sprintf('Longest audited fill: %s using %s.', ...
            icemodel.plot.formatDuration( ...
            hours(value.focus_end - value.focus_start)), value.method);
      case "climatology_variability"
         text = sprintf([ ...
            'Held-out selection variability ratio is %.3g; distance from ' ...
            'one is the declared compression/inflation diagnostic.'], ...
            value.score);
      case "boundary_behavior"
         text = sprintf([ ...
            'Largest positive admitted boundary-jump rate is %.3g; this ' ...
            'is diagnostic evidence, not a readiness verdict input.'], ...
            value.score);
      case "hourly_disaggregation"
         text = sprintf([ ...
            'Selected proxy span is %s; the lower panel exposes the ' ...
            '15-minute mean-preserving expansion.'], ...
            icemodel.plot.formatDuration( ...
            hours(value.focus_end - value.focus_start)));
      case "residual_missing"
         text = sprintf([ ...
            'Largest residual share is %.3g%%; the marked %s run remains ' ...
            'unfilled by design.'], value.score, ...
            icemodel.plot.formatDuration( ...
            hours(value.focus_end - value.focus_start)));
      case "native_sensor_flat_run"
         text = sprintf([ ...
            'The longest corroborated native sensor-failure signature ' ...
            'lasts %.0f days. The dark trace is the preserved staged ' ...
            'source; the orange trace is the separately reconstructed ' ...
            'product after this channel was excluded from observational ' ...
            'truth.'], value.score);
   end
end

function ledger = emptyFigureLedger()
   %EMPTYFIGURELEDGER Schema-stable empty figure ledger.
   ledger = table('Size', [0 10], 'VariableTypes', {'string', ...
      'string', 'string', 'datetime', 'datetime', 'double', ...
      'double', 'logical', 'logical', 'string'}, ...
      'VariableNames', {'site', 'channel', 'method', 'gap_start', ...
      'gap_end', 'duration_hours', 'segments_shown', ...
      'has_before_context', 'has_after_context', 'file'});
end

function met = loadMet(filename, site, product)
   %LOADMET Load and validate one byte-verified met artifact.
   stem = site + "_" + product;
   if ~isfile(filename)
      error('icemodel:report:buildGapFillReport:missingMet', ...
         'no byte-verified met file for %s: %s', stem, filename);
   end
   [~, basename] = fileparts(filename);
   expected_prefix = "met_" + site + "_" + product + "_";
   if ~startsWith(lower(string(basename)), lower(expected_prefix))
      error('icemodel:report:buildGapFillReport:artifactIdentityMismatch', ...
         'met filename does not identify requested artifact %s: %s', ...
         stem, filename);
   end
   met = icemodel.forcing.reconstruct.loadWidestTimetable(dir(filename));
   if isempty(met)
      error('icemodel:report:buildGapFillReport:invalidMet', ...
         'met file for %s contains no timetable', stem);
   end
   metadata = met.Properties.UserData;
   has_site = isstruct(metadata) && isscalar(metadata) ...
      && isfield(metadata, 'site') ...
      && isscalar(string(metadata.site));
   if ~has_site || ...
         icemodel.forcing.helpers.normalizedFileToken(metadata.site) ~= ...
         icemodel.forcing.helpers.normalizedFileToken(site)
      error('icemodel:report:buildGapFillReport:artifactIdentityMismatch', ...
         'met metadata does not identify requested site %s: %s', ...
         site, filename);
   end
   has_product = isfield(metadata, 'gapfill_product');
   if has_product
      actual_product = string(metadata.gapfill_product);
   else
      actual_product = "";
   end
   if (~has_product && product == "promice_filled") ...
         || (has_product && (~isscalar(actual_product) ...
         || actual_product ~= product))
      error('icemodel:report:buildGapFillReport:artifactIdentityMismatch', ...
         'met metadata does not identify requested product %s: %s', ...
         product, filename);
   end
end

function qmd_file = writeQmd(report_dir, qa_dir, fig_dir, ...
      fill_by_family, readiness, method_diagnostics, ...
      interpretation_catalog, sites, residual_gaps, ...
      readiness_blockers, absent_products, ropts, ...
      context_bounds, max_segments_per_method, ...
      max_detail_figures_per_site)
   %WRITEQMD Emit the fixed-structure Quarto document.
   qmd_file = fullfile(report_dir, 'gapfill-report.qmd');
   fid = fopen(qmd_file, 'w');
   cleaner = onCleanup(@() fclose(fid));
   w = @(varargin) fprintf(fid, [varargin{1} '\n'], varargin{2:end});

   % Every summary count derives from the policy verdict explicitly so
   % the plain-language result and the ledger can never disagree.
   n_ready = nnz(readiness.policy_verdict == "ready");
    n_blocked = nnz(readiness.policy_verdict == "not_forcing_ready");
    n_out = nnz(readiness.policy_verdict == "out_of_policy_window");
    n_low_confidence = 0;
    if ismember('worth_filling', readiness.Properties.VariableNames)
       n_low_confidence = nnz(~readiness.worth_filling);
    end
    n_no_window = nnz(isnat(readiness.window_start));
   w('---');
   w('title: "PROMICE gap-filled forcing: before/after report"');
   w('format:');
   w('  pdf:');
   w('    toc: true');
   w('  html:');
   w('    toc: true');
   w('    embed-resources: false');
   w('---');
   w('');
   w('# Executive Summary');
   w('');
   w(['Gap-filled (`promice_filled`) forcing was produced for %d PROMICE ' ...
       'stations with per-sample provenance. Of %d station-years: %d are ' ...
       'forcing-ready, %d remain not ready, %d fall outside an available ' ...
       'staged-proxy acceptance window on staging-era boundaries alone, ' ...
       'and %d carry the low-native-coverage confidence advisory. ' ...
       'The advisory never changes readiness. No staged-proxy acceptance ' ...
       'window is available for %d station-years; their absolute verdicts ' ...
       'are reported unchanged. Every ' ...
       'verdict and reason is recorded in the readiness ledger. Methods ' ...
       'were admitted per station, channel, and gap-duration stratum ' ...
       'through held-out blocked validation only. %d staged native ' ...
       'stations have no filled product at all; Results names them and ' ...
       'the refusal reason.'], numel(sites), ...
       height(readiness), n_ready, n_blocked, n_out, n_low_confidence, ...
       n_no_window, height(absent_products));
   w('');
   w(['**Plain-language result:** produced products contain %d residual ' ...
      'required-forcing runs and %d blocked station-years. %d staged native ' ...
      'stations have no product. A station appears in the no-product table ' ...
      'only when no validated staged MAR/MERRA proxy window exists or its ' ...
      'native record does not overlap that window; those products cannot ' ...
      'be filled from the available data.'], height(residual_gaps), ...
      height(readiness_blockers), height(absent_products));
   w('');
   w('# Background');
   w('');
   w(['Native PROMICE preserves sensor outages as missing data; the v1.1 ' ...
      'forcing gate showed 48-168 h required-channel gaps block every ' ...
      'formal case. The approved reconstruction policy ' ...
       '(icemodel/+icemodel/+forcing/+reconstruct/POLICY.md) fills the separately named canonical ' ...
      'product: bounded short-gap interpolation, donor station transfer ' ...
      '(PROMICE, K-transect, GC-Net observed samples), calibrated MAR ' ...
      'proxies, and climatology, never overwriting true native samples. ' ...
      'PROMICE albedo produced by the legacy builder is compared against ' ...
      'the raw L3 finite mask so interpolated/winter-stamped values are ' ...
      'not mislabeled as observations.']);
   w('');
   w('# Methods');
   w('');
   % Gate numbers come from the producer-pinned reconstruction options.
   w(sprintf(['Every candidate method was fitted on selection years ' ...
      'and graded on held-out blocked synthetic gaps drawn from each ' ...
      'station''s real gap-duration census; admission required meeting ' ...
        'bias caps, a %.0f%%%% RMSE improvement over persistence for gaps ' ...
        'up to %g h and station day-of-year climatology for longer gaps, ' ...
        'with wind donors additionally required to beat calibrated MAR ' ...
        'wind in the same stratum, ' ...
        'a held-out variability ratio within %.2f--%.2f, zero ' ...
      'physical-bound violations, and %.0f%%%% coverage. A method is not ' ...
      'extrapolated to a real gap longer than the longest held-out gap ' ...
      'that admitted it. Final metrics come from evaluation-year draws. ' ...
        'Each station has at most %d detail figures, selected by ' ...
        'interleaving duration ranks within channel x method-family ' ...
        'groups. Every selected figure contains one method and up to %d ' ...
        'audited gaps spanning that method''s observed duration range; ' ...
        'each panel''s context pad equals its filled period, ' ...
      'floored at %g d and capped at %g d.'], ...
       100 * ropts.rmse_improvement, ropts.cap_hours, ...
       ropts.min_variability_ratio, ropts.max_variability_ratio, ...
       100 * ropts.min_coverage, ...
      max_detail_figures_per_site, max_segments_per_method, ...
      context_bounds(1), context_bounds(2)));
   w('');
   w(['Within each detail panel, only the panel''s own method renders in ' ...
      'the accent color; reconstructions by any other method inside the ' ...
      'context window render muted grey, keyed by the per-sample ' ...
      'provenance registry and the plan audit (POLICY D-31), so a panel ' ...
      'never visually claims a foreign fill.']);
   w('');
   w('### Producer-pinned reconstruction options');
   w('');
   w(['This complete option record is byte-independent producer metadata; ' ...
      'every selected station has an identical record.']);
   w('');
   w('```json');
   w('%s', jsonencode(ropts, PrettyPrint=true));
   w('```');
   w('');
   w(['RMSE and bias are retained as point-error criteria, but they do ' ...
      'not diagnose compressed variability. Correlation and the standard-' ...
      'deviation ratio are therefore reported separately, following the ' ...
      '[component-based model-evaluation rationale](https://doi.org/' ...
      '10.5194/hess-23-4323-2019); no opaque aggregate efficiency score ' ...
      'is used. This also responds to documented variance distortion from ' ...
      '[climate-series gap filling](https://doi.org/10.1175/' ...
      'JCLI-D-19-0244.1).']);
   w('');
   w('# Results');
   w('');
   w('## Cohort readiness');
   w('');
   w(['Policy verdicts per station-year, graded from the producer-pinned ' ...
      'readiness ledger (completeness plus the shared scalar physical ' ...
      'bounds):']);
   w('');
   writeMarkdownTable(w, verdictSummary(readiness));
   w('');
   w('## Why anything remains unfilled');
   w('');
   if isempty(residual_gaps) && isempty(readiness_blockers)
      w(['**No required forcing gap remains unfilled inside any produced ' ...
         'product, and every produced station-year is ready for both ' ...
         'IceModel and the snow model.**']);
   else
      if isempty(residual_gaps)
         w(['No required forcing sample remains missing inside a produced ' ...
            'product. The station-year table below records finite-but-' ...
            'invalid values, snowfall-input limitations, or years outside ' ...
            'the staged proxy window.']);
      else
         w(['Every remaining missing run is listed below. The reason is ' ...
            'the producer''s reconciled final-tier denial; the report ' ...
            'refuses to build if these rows do not exactly match the ' ...
            'shipped missing samples.']);
         w('');
         writeMarkdownTable(w, residual_gaps);
      end
      if ~isempty(readiness_blockers)
         w('');
         w(['These are all produced station-years that are not ready for ' ...
            'IceModel, the snow model, or both:']);
         w('');
         writeMarkdownTable(w, readiness_blockers);
      end
   end
   w('');
   w(['Machine-readable copies: ' ...
      '[residual gaps](%s) and [readiness blockers](%s).'], ...
      relativeFigure(report_dir, qa_dir, ...
      'gapfill_residual_gaps.csv'), ...
      relativeFigure(report_dir, qa_dir, ...
      'gapfill_readiness_blockers.csv'));
   w('');
   w(['Relational diagnostics (POLICY A15/D-28): swd-versus-TOA and ' ...
      'swu-versus-swd exceedances are DIAGNOSTIC evidence reported with ' ...
      'the ledger, never verdict inputs — twilight diffuse light is real ' ...
       'incident energy the TOA model cannot represent, and native sensor ' ...
       'noise must not brand usable years not-ready. The same relational ' ...
       'rules remain hard admission gates for fill candidates except the ' ...
       'explicit D-32 post-final SWD seam bridge. That bridge keeps scalar ' ...
       'bounds hard but reports the solar relation diagnostically because ' ...
       'an instantaneous TOA cap can create the discontinuity it repairs.']);
   w('');
   w('## Fill volume by method family');
   w('');
   w(['Percent of all channel-samples per station, across the ' ...
      'provenance-paired science channels, keyed by the per-sample ' ...
      'provenance registry the engine stamps. Columns that are zero ' ...
      'across the whole cohort are omitted from this view; the linked ' ...
      'CSV keeps the full registry schema.']);
   w('');
   writeMarkdownTable(w, prunedFamilyTable(fill_by_family));
   w('');
   w('## Method admission and held-out skill');
   w('');
   w(['Cohort aggregate of the per-station plan diagnostics: admitted ' ...
      'strata and denied candidates per channel and method family, with ' ...
      'the median held-out RMSE improvement over each stratum''s policy ' ...
      'baseline (100 x (1 - rmse / baseline_rmse) on selection draws). ' ...
      'The complete per-stratum evidence remains in the method ' ...
      'diagnostics CSV.']);
   w('');
   writeMarkdownTable(w, admissionSummary(method_diagnostics));
   w('');
   w('## Stations without filled products');
   w('');
   if isempty(absent_products)
      w(['Every staged native station in the selected data root has a ' ...
         'filled product on disk.']);
   else
      w(['These staged native stations have no filled product. Their exact ' ...
         'native and proxy spans are independently verified below; an ' ...
         'absent product whose spans overlap would stop report generation ' ...
         'as unexplained rather than receive an invented refusal reason.']);
      w('');
      writeMarkdownTable(w, absent_products);
   end
   w('');
   w('## Scientific interpretation and decision-sensitive outcomes');
   w('');
   w(['The category catalog is fixed before station/window selection. ' ...
      'Every selector operates on producer-pinned artifacts; a category ' ...
      'that is not materially present remains an explicit absent row. ' ...
      'Statuses distinguish expected or accepted behavior, diagnostic ' ...
      'limitations that deserve inspection, and unresolved missing data.']);
   w('');
   % Keep the rendered index readable in portrait PDF. Each category's
   % method, score, rule, mechanism, and policy basis appear immediately
   % below, while the linked CSV retains the complete machine-readable row.
   catalog_view = interpretation_catalog(:, {'category', 'present', ...
      'status', 'site', 'channel'});
   writeMarkdownTable(w, catalog_view);
   w('');
   w('Complete interpretation catalog: [%s](%s).', ...
      'gapfill_interpretation_catalog.csv', relativeFigure( ...
      report_dir, qa_dir, 'gapfill_interpretation_catalog.csv'));
   w('');
   for k = 1:height(interpretation_catalog)
      item = interpretation_catalog(k, :);
      heading = strrep(string(item.category), "_", " ");
      w('### %s', heading);
      w('');
      w('**Selection rule.** %s', string(item.selection_rule));
      w('');
      w('**Mechanism.** %s', string(item.mechanism));
      w('');
      if item.present
         w('**Window.** %s to %s UTC.', ...
            string(item.window_start), string(item.window_end));
         w('');
      end
      w('**Interpretation.** %s', string(item.interpretation));
      w('');
      w('**Policy basis.** %s', string(item.policy_basis));
      w('');
      if item.present
         figure_path = relativeFigure(report_dir, fig_dir, ...
            string(item.figure));
         w('[![%s](%s)](%s)', heading, figure_path, figure_path);
         w('');
      end
   end
    w('## Summary tables');
    w('');
    w(['Machine-readable tables: [channel summary](%s), ' ...
       '[fill volume by family](%s), [method diagnostics](%s), ' ...
       '[scientific interpretation catalog](%s), [residual gaps](%s), ' ...
       '[readiness blockers](%s), [absent products](%s), ' ...
       '[readiness ledger](%s), and [figure ledger](%s).'], ...
       relativeFigure(report_dir, qa_dir, 'gapfill_summary.csv'), ...
       relativeFigure(report_dir, qa_dir, ...
       'gapfill_fill_by_family.csv'), ...
       relativeFigure(report_dir, qa_dir, ...
       'gapfill_method_diagnostics.csv'), ...
       relativeFigure(report_dir, qa_dir, ...
       'gapfill_interpretation_catalog.csv'), ...
       relativeFigure(report_dir, qa_dir, ...
       'gapfill_residual_gaps.csv'), ...
       relativeFigure(report_dir, qa_dir, ...
       'gapfill_readiness_blockers.csv'), ...
       relativeFigure(report_dir, qa_dir, ...
       'gapfill_absent_products.csv'), ...
       relativeFigure(report_dir, qa_dir, ...
       'promice_filled_readiness.csv'), ...
       relativeFigure(report_dir, qa_dir, ...
       'gapfill_figure_ledger.csv'));
    w('');
   w('# Summary');
   w('');
    if isempty(residual_gaps) && isempty(readiness_blockers)
       w(['For the selected cohort, all %d produced station-years are ' ...
          'ready for both IceModel and the snow model, with zero required ' ...
          'forcing samples left missing. %d low-native-coverage advisories ' ...
          'remain scientific-confidence notes only.'], ...
          n_ready, n_low_confidence);
    else
       w(['For the selected cohort, the staged ledger records %d ready, %d ' ...
          'not ready, %d outside an available acceptance window, and %d ' ...
          'low-native-coverage advisories. The exhaustive tables above ' ...
          'list %d remaining missing runs and %d blocked station-years ' ...
          'with their exact reasons. %d station-years have no staged-' ...
          'proxy acceptance window and retain their absolute verdict.'], ...
          n_ready, n_blocked, n_out, n_low_confidence, ...
          height(residual_gaps), height(readiness_blockers), n_no_window);
    end
   w('');
   w('# Appendices');
   w('');
   fig_ledger = readTableSafe(fullfile(qa_dir, 'gapfill_figure_ledger.csv'));
   is_overview = string(fig_ledger.method) == "overview";
   appendix_sites = sites;
   w('## A. Station overview figures');
   w('');
   w(['Each station''s full-period overview stacks its provenance-paired ' ...
      'science channels (observed dark, filled accent orange; POLICY ' ...
      'A14/D-19) — the one deliberate full-period exception to the ' ...
      'windowed-figure contract. All %d overview figures live in the ' ...
      'overview/ figure folder and are ledgered in ' ...
      'gapfill_figure_ledger.csv beside the detail set.'], ...
      nnz(is_overview));
   w('');
   w('Stations:');
   w('');
   for site = reshape(appendix_sites, 1, [])
      w('- [%s](#overview-%s)', upper(site), lower(site));
   end
   w('');
   for site = reshape(appendix_sites, 1, [])
      w('### %s {#overview-%s}', upper(site), lower(site));
      w('');
      site_figures = fig_ledger(is_overview ...
         & string(fig_ledger.site) == site, :);
      % Link each embedded preview to its ledgered full-resolution PNG.
      for k = 1:height(site_figures)
         figure_path = relativeFigure(report_dir, fig_dir, ...
            string(site_figures.file(k)));
         w('[![%s %s %s](%s)](%s)', upper(site), ...
            string(site_figures.channel(k)), ...
            string(site_figures.method(k)), ...
            figure_path, figure_path);
         w('');
      end
   end
   w('## B. Method-detail companion report');
   w('');
   w(['All method/gap appendix figures and their per-station channel ' ...
      'summaries and candidate decisions are kept in a separate ' ...
      'companion report so this main report remains responsive: ' ...
      '[PDF](gapfill-detail-report.pdf) and ' ...
      '[HTML](gapfill-detail-report.html).']);
   w('');
   w('## C. Readiness ledger');
   w('');
   writeMarkdownTable(w, readiness);
   w('');
   w('Machine-readable copy: promice_filled_readiness.csv.');
   clear cleaner
end

function qmd_file = writeDetailQmd(report_dir, qa_dir, fig_dir, ...
      summary, method_diagnostics, sites, max_detail_figures_per_site)
   %WRITEDETAILQMD Emit the companion method/gap appendix report.
   %
   % The main report keeps the scientific Results and station overviews
   % responsive. This companion owns every bounded method-detail figure
   % plus compact per-station tables; the complete candidate ledger stays
   % available as CSV instead of expanding to tens of thousands of rows.
   qmd_file = fullfile(report_dir, 'gapfill-detail-report.qmd');
   fid = fopen(qmd_file, 'w');
   cleaner = onCleanup(@() fclose(fid));
   w = @(varargin) fprintf(fid, [varargin{1} '\n'], varargin{2:end});

   % Select only method/gap figures; overview and interpretation figures
   % remain owned by the main scientific report.
   fig_ledger = readTableSafe(fullfile(qa_dir, ...
      'gapfill_figure_ledger.csv'));
   is_overview = string(fig_ledger.method) == "overview";
   is_interpretation = startsWith(string(fig_ledger.method), ...
      "interpretation:");
   detail = fig_ledger(~is_overview & ~is_interpretation, :);

   % PDF comes first as the default review artifact; the builder renders
   % both formats so links from the main report always resolve.
   w('---');
   w('title: "PROMICE gap-filled forcing: method-detail appendix"');
   w('format:');
   w('  pdf:');
   w('    toc: true');
   w('  html:');
   w('    toc: true');
   w('    embed-resources: false');
   w('---');
   w('');
   w('# Method-detail appendix');
   w('');
   w(['This companion contains all %d grouped method/gap figures. Each ' ...
      'station contributes at most %d figures selected by the declared ' ...
      'channel x method-family duration-rank rule. Every figure combines ' ...
      'a bounded representative set of gaps filled by one method; only ' ...
      'that method renders in the accent color, while fills by other ' ...
      'methods in the context window render muted grey (POLICY D-31).'], ...
      height(detail), max_detail_figures_per_site);
   w('');
   w('[Return to the main report](gapfill-report.html).');
   w('');
   w(['Full candidate decision ledger, including every denial and all ' ...
      'selection/evaluation diagnostics: [%s](%s).'], ...
      'gapfill_method_diagnostics.csv', relativeFigure(report_dir, qa_dir, ...
      'gapfill_method_diagnostics.csv'));
   w('');
   w('Stations:');
   w('');
   for site = reshape(sites, 1, [])
      w('- [%s](#detail-%s)', upper(site), lower(site));
   end
   w('');

   % Each station block keeps figures beside the exact channel summary
   % and candidate evidence that explain them.
   for site = reshape(sites, 1, [])
      w('## %s {#detail-%s}', upper(site), lower(site));
      w('');
      site_figures = detail(string(detail.site) == site, :);
      for k = 1:height(site_figures)
         figure_path = relativeFigure(report_dir, fig_dir, ...
            string(site_figures.file(k)));
         w('[![%s %s %s](%s)](%s)', upper(site), ...
            string(site_figures.channel(k)), ...
            string(site_figures.method(k)), ...
            figure_path, figure_path);
         w('');
      end
      w('Channel summary:');
      w('');
      in_summary = string(summary.site) == site;
      writeMarkdownTable(w, summary(in_summary, :));
      w('');
      w('Admitted candidate decisions (compact view):');
      w('');
      in_diagnostics = string(method_diagnostics.site) == site ...
         & string(method_diagnostics.decision) == "admitted";
      compact_columns = {'channel', 'candidate', 'season', 'bucket', ...
         'max_validated_hours', 'selection_n', 'selection_rmse', ...
         'selection_baseline_rmse'};
      writeMarkdownTable(w, method_diagnostics(in_diagnostics, ...
         compact_columns));
      w('');
   end
   w('Machine-readable figure ledger: [%s](%s).', ...
      'gapfill_figure_ledger.csv', relativeFigure(report_dir, qa_dir, ...
      'gapfill_figure_ledger.csv'));
   clear cleaner
end

function pruned = prunedFamilyTable(fill_by_family)
   %PRUNEDFAMILYTABLE Drop family columns that are zero across the cohort.
   % The QMD view stays compact while the CSV keeps the full registry
   % schema; the pruning rule is stated beside the rendered table.
   pruned = fill_by_family;
   if isempty(pruned)
      return
   end
   names = string(pruned.Properties.VariableNames);
   drop = false(size(names));
   for k = 1:numel(names)
      value = pruned.(names(k));
      drop(k) = isnumeric(value) && all(value == 0);
   end
   pruned(:, drop) = [];
end

function t = readTableSafe(filename)
   %READTABLESAFE Read a CSV that this generator wrote earlier in the pass.
   % The delimiter is pinned: our own reason/detail strings carry many
   % semicolons, which can flip readtable's delimiter auto-detection.
   t = readtable(filename, 'Delimiter', ',', 'TextType', 'string');
end

function rel = relativeFigure(report_dir, fig_dir, filename)
   %RELATIVEFIGURE Report-relative figure link (no copying or embedding).
   report_parts = split(string(report_dir), filesep);
   figure_parts = split(string(fig_dir), filesep);
   report_parts(report_parts == "") = [];
   figure_parts(figure_parts == "") = [];
   n_common = 0;
   while n_common < min(numel(report_parts), numel(figure_parts)) ...
         && report_parts(n_common + 1) == figure_parts(n_common + 1)
      n_common = n_common + 1;
   end
   up = repmat("..", numel(report_parts) - n_common, 1);
   down = [figure_parts((n_common + 1):end); string(filename)];
   rel = strjoin([up; down], "/");
end

function writeMarkdownTable(w, t)
   %WRITEMARKDOWNTABLE Emit a small table as GitHub-flavored markdown.
   names = t.Properties.VariableNames;
   w('| %s |', strjoin(names, ' | '));
   w('|%s', repmat('---|', 1, numel(names)));
   for r = 1:height(t)
      cells = strings(1, numel(names));
      for c = 1:numel(names)
         value = t.(names{c})(r);
         if isnumeric(value)
            cells(c) = string(round(double(value), 3));
         else
            cells(c) = string(value);
         end
         % Empty reasons round-trip through CSV as missing strings; a
         % markdown cell renders them as blank, never as an fprintf error.
         if ismissing(cells(c))
            cells(c) = "";
         end
      end
      w('| %s |', strjoin(cells, ' | '));
   end
end

function [pdf_file, html_file] = renderQmd(qmd_file)
   %RENDERQMD Render PDF and HTML report formats beside the preview tree.
   [status, output] = system("quarto render " + ...
      icemodel.verification.helpers.shellQuote(qmd_file) + " --to all");
   if status ~= 0
      error('icemodel:report:buildGapFillReport:renderFailed', ...
         'quarto render failed: %s', output);
   end
   pdf_file = strrep(qmd_file, '.qmd', '.pdf');
   html_file = strrep(qmd_file, '.qmd', '.html');
end

function inputs = verifySiteInputs(site, qa_dir, filled_dir)
   %VERIFYSITEINPUTS Verify the producer-recorded size and SHA-256.
   manifest_file = fullfile(qa_dir, 'plans', ...
      sprintf('%s-report-inputs.json', site));
   if ~isfile(manifest_file)
      error('icemodel:report:buildGapFillReport:missingInputManifest', ...
         'missing report-input manifest for %s: %s', site, manifest_file);
   end
   manifest = jsondecode(fileread(manifest_file));
   if string(manifest.site) ~= site
      error('icemodel:report:buildGapFillReport:inputIdentityMismatch', ...
         'report-input manifest site mismatch for %s', site);
   end

   required = ["native", "filled", "plan", "readiness"];
   roles = string({manifest.artifacts.role});
   inputs = struct();
   if ~isfield(manifest, 'path_base') ...
         || string(manifest.path_base) ~= "selected_data_root"
      error('icemodel:report:buildGapFillReport:inputIdentityMismatch', ...
         'report-input manifest for %s lacks a relocatable path base', site);
   end
   [data_root, ~] = ...
      icemodel.forcing.reconstruct.selectedDataRoot(filled_dir);
   resolved_paths = strings(numel(manifest.artifacts), 1);
   for k = 1:numel(manifest.artifacts)
      record = manifest.artifacts(k);
      relative = string(record.path);
      relative_file = java.io.File(char(relative));
      if relative_file.isAbsolute()
         error('icemodel:report:buildGapFillReport:inputPathOutsideRoot', ...
            '%s report-input path must be relative: %s', site, relative);
      end
      pathname = string(fullfile(data_root, relative));
      if ~icemodel.internal.isPathInside(pathname, data_root)
         error('icemodel:report:buildGapFillReport:inputPathOutsideRoot', ...
            '%s report-input path escapes selected root %s: %s', ...
            site, data_root, relative);
      end
      resolved_paths(k) = pathname;
      role = string(record.role);
       switch role
          case "filled"
             allowed_root = filled_dir;
          case {"native", "proxy_window"}
             allowed_root = string(fileparts(filled_dir));
          case {"plan", "readiness"}
             allowed_root = qa_dir;
          otherwise
             error( ...
                'icemodel:report:buildGapFillReport:inputIdentityMismatch', ...
                'unknown %s report-input role for %s', role, site);
       end
      if ~icemodel.internal.isPathInside(pathname, allowed_root)
         error('icemodel:report:buildGapFillReport:inputPathOutsideRoot', ...
            '%s %s artifact is outside selected root %s: %s', ...
            site, role, allowed_root, pathname);
      end
       [~, artifact_name, artifact_ext] = fileparts(pathname);
       if role == "proxy_window" && ~startsWith( ...
             string(artifact_name) + string(artifact_ext), ...
             "met_" + site + "_")
          error('icemodel:report:buildGapFillReport:inputIdentityMismatch', ...
             'proxy-window artifact does not belong to %s: %s', ...
             site, pathname);
       end
      info = dir(pathname);
      if isempty(info) || info.bytes ~= double(record.bytes) ...
            || lower(icemodel.verification.setup.fileSha256(pathname)) ...
            ~= lower(string(record.sha256))
         error('icemodel:report:buildGapFillReport:inputIdentityMismatch', ...
            'byte identity mismatch for %s %s artifact: %s', ...
            site, string(record.role), pathname);
      end
   end
   for role = required
      match = find(roles == role);
      if ~isscalar(match)
         error('icemodel:report:buildGapFillReport:inputIdentityMismatch', ...
            'report-input manifest for %s must contain one %s artifact', ...
            site, role);
      end
      inputs.(char(role)) = resolved_paths(match);
   end
   if ~isfield(manifest, 'acceptance_window')
      error('icemodel:report:buildGapFillReport:inputIdentityMismatch', ...
         'report-input manifest for %s lacks an acceptance window', site);
   end
   [window_start, window_end] = ...
      icemodel.verification.setup.periodBounds(manifest.acceptance_window);
   window_start.TimeZone = 'UTC';
   window_end.TimeZone = 'UTC';
   inputs.acceptance_window = [window_start, window_end];
   saved = load(inputs.plan, 'plan_record');
   inputs.reconstruction_options = ...
      saved.plan_record.reconstruction_options;
end

function ropts = pinnedReconstructionOptions(inputs)
   %PINNEDRECONSTRUCTIONOPTIONS Require one report policy across the cohort.
   ropts = inputs{1}.reconstruction_options;
   for k = 2:numel(inputs)
      candidate = inputs{k}.reconstruction_options;
      if ~isequaln(ropts, candidate)
         error('icemodel:report:buildGapFillReport:inconsistentPolicy', ...
            'report cohort contains incompatible reconstruction policies');
      end
   end
end

function pathname = defaultDir(value, repo, relative)
   %DEFAULTDIR Resolve one directory default under the repo root.
   pathname = value;
   if pathname == ""
      pathname = string(fullfile(repo, relative));
   end
end
