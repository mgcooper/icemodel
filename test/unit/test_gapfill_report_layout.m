function tests = test_gapfill_report_layout
   %TEST_GAPFILL_REPORT_LAYOUT Layout and Results contract of the gap report.
   %
   % Covers the POLICY D-31 report restructure of
   % icemodel.verification.report.buildGapFillReport:
   %  - overview/detail/interpretation figure folder split and
   %    subfolder-relative ledger paths, including stale cleanup;
   %  - the bounded, deterministic per-station detail selection;
   %  - the overview-only main appendix and separate method-detail report
   %    with clickable station anchors;
   %  - the artifact-derived Results tables (cohort verdicts, fill volume
   %    by provenance family, admission/held-out skill aggregate);
   %  - the explicit statement of staged native stations without any
   %    filled product (windowRecordDisjoint refusals);
   %  - the method-accent layer split (methodFillLayers color logic).
   %
   % Fixtures are hand-built saved artifacts under an isolated temp tree
   % (manifest, filled/native met with provenance channels, plan/audit,
   % readiness ledger); the generator never sees real preview data.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % One fixture tree and one report build feed every layout assertion;
   % the build itself doubles as the fixture-based generator smoke.
   root = string(tempname());
   mkdir(root);
   testCase.TestData.root = root;
   writeLayoutFixture(root);
   seedStaleFigures(root);
   testCase.TestData.report = ...
      icemodel.verification.report.buildGapFillReport( ...
      sites="all", render=false, ...
      filled_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      fig_dir=fullfile(root, 'figures'), ...
      report_dir=fullfile(root, 'report'));
   testCase.TestData.qmd = ...
      string(fileread(testCase.TestData.report.qmd_file));
   testCase.TestData.detail_qmd = ...
      string(fileread(testCase.TestData.report.detail_qmd_file));
end

function teardownOnce(testCase)
   % The fixture tree is fully disposable.
   if isfolder(testCase.TestData.root)
      rmdir(testCase.TestData.root, 's');
   end
end

%% folder split and ledger paths

function test_overview_and_detail_folders_split(testCase)
   % POLICY D-31: overview figures land in overview/, per-method detail
   % figures in detail/, and nothing new lands flat in the figure root.
   root = testCase.TestData.root;
   overview_file = fullfile(root, 'figures', 'overview', ...
      'tsta_overview.png');
   testCase.verifyTrue(isfile(overview_file));
   % A ledgered white canvas is not a rendered overview.
   pixels = imread(overview_file);
   testCase.verifyGreaterThan( ...
      double(max(pixels, [], 'all')) - double(min(pixels, [], 'all')), 0);
   testCase.verifyTrue(isfile(fullfile(root, 'figures', 'detail', ...
      'tsta_tair_bounded_interp_gapfill.png')));
   testCase.verifyTrue(isfile(fullfile(root, 'figures', 'detail', ...
      'tsta_tair_donor_dsta_gapfill.png')));
   testCase.verifyTrue(isfile(fullfile(root, 'figures', ...
      'interpretation', ...
      'tsta_tair_long_gap_limited_context.png')));
   % The pre-split flat location must not receive new exports.
   testCase.verifyFalse(isfile(fullfile(root, 'figures', ...
      'tsta_overview.png')));
   testCase.verifyFalse(isfile(fullfile(root, 'figures', ...
      'tsta_tair_bounded_interp_gapfill.png')));
end

function test_ledger_paths_reconcile_with_folders(testCase)
   % Ledger file values carry the subfolder, stay unique, and every row
   % resolves to a PNG on disk — exact reconciliation over the split.
   root = testCase.TestData.root;
   ledger = testCase.TestData.report.figure_ledger;
   files = string(ledger.file);
   testCase.verifyTrue(all(startsWith(files, "overview/") ...
      | startsWith(files, "detail/") ...
      | startsWith(files, "interpretation/")));
   testCase.verifyEqual(numel(unique(files)), height(ledger));
   for file = reshape(files, 1, [])
      testCase.verifyTrue(isfile(fullfile(root, 'figures', file)));
   end
   % The station overview leads its site block and lives in overview/.
   overview_rows = ledger.method == "overview";
   testCase.verifyEqual(nnz(overview_rows), 1);
   returned = string(ledger.file(overview_rows));
   expected = "overview/tsta_overview.png";
   testCase.verifyEqual(returned, expected);
   site_rows = find(ledger.site == "tsta");
   testCase.verifyEqual(string(ledger.method(site_rows(1))), "overview");
   % Each ledger class stays in its owned folder.
   interpretation_rows = startsWith(ledger.method, "interpretation:");
   detail_rows = ~overview_rows & ~interpretation_rows;
   testCase.verifyTrue(all(startsWith(files(detail_rows), "detail/")));
   testCase.verifyTrue(all(startsWith(files(interpretation_rows), ...
      "interpretation/")));
end

function test_stale_and_foreign_figures_handled(testCase)
   % Selected-site leftovers vanish in the same transaction (including
   % legacy flat-layout PNGs); foreign sites' figures survive untouched.
   root = testCase.TestData.root;
   testCase.verifyFalse(isfile(fullfile(root, 'figures', 'detail', ...
      'tsta_tair_stale_gapfill.png')));
   testCase.verifyFalse(isfile(fullfile(root, 'figures', ...
      'tsta_tair_legacyflat_gapfill.png')));
   testCase.verifyTrue(isfile(fullfile(root, 'figures', ...
      'ghost_legacyflat.png')));
   testCase.verifyTrue(isfile(fullfile(root, 'figures', 'overview', ...
      'ghost_overview.png')));
end

%% appendix sections and anchors

function test_appendix_sections_have_clickable_anchors(testCase)
   % The main appendix embeds only station overviews; all method/gap
   % panels and their station tables live in the companion report.
   text = testCase.TestData.qmd;
   detail_text = testCase.TestData.detail_qmd;
   for section = ["# Executive Summary", "# Background", "# Methods", ...
         "# Results", "# Summary", "# Appendices"]
      testCase.verifyTrue(contains(text, section));
   end
   testCase.verifyTrue(contains(text, "## A. Station overview figures"));
   testCase.verifyTrue(contains(text, ...
      "## B. Method-detail companion report"));
   testCase.verifyTrue(contains(text, "## C. Readiness ledger"));
   testCase.verifyTrue(contains(text, "- [TSTA](#overview-tsta)"));
   testCase.verifyTrue(contains(text, "### TSTA {#overview-tsta}"));
   testCase.verifyFalse(contains(text, "- [TSTA](#detail-tsta)"));
   testCase.verifyTrue(contains(text, "gapfill-detail-report.pdf"));
   testCase.verifyTrue(contains(text, "gapfill-detail-report.html"));
   testCase.verifyTrue(contains(detail_text, ...
      "- [TSTA](#detail-tsta)"));
   testCase.verifyTrue(contains(detail_text, ...
      "## TSTA {#detail-tsta}"));

   % Each figure resolves in the report that owns its display class.
   ledger = testCase.TestData.report.figure_ledger;
   for k = 1:height(ledger)
      file = string(ledger.file(k));
      figure_path = "../figures/" + file;
      is_detail = string(ledger.method(k)) ~= "overview" ...
         && ~startsWith(string(ledger.method(k)), "interpretation:");
      if is_detail
         testCase.verifyTrue(contains(detail_text, ...
            "](" + figure_path + ")](" + figure_path + ")"));
         testCase.verifyFalse(contains(text, figure_path));
      else
         testCase.verifyTrue(contains(text, ...
            "](" + figure_path + ")](" + figure_path + ")"));
         testCase.verifyFalse(contains(detail_text, figure_path));
      end
   end
   % Compact per-station tables move with the figures rather than bloating
   % main; the exhaustive denied-candidate evidence remains a linked CSV.
   testCase.verifyFalse(contains(text, "Channel summary:"));
   testCase.verifyFalse(contains(text, "Admitted candidate decisions"));
   testCase.verifyTrue(contains(detail_text, "Channel summary:"));
   testCase.verifyTrue(contains(detail_text, ...
      "Admitted candidate decisions (compact view):"));
   testCase.verifyTrue(contains(detail_text, ...
      "Full candidate decision ledger"));
   testCase.verifyTrue(contains(detail_text, ...
      "gapfill_method_diagnostics.csv"));
   testCase.verifyFalse(contains(detail_text, ...
      "insufficient staged overlap"));
end

%% Results section

function test_results_tables_present(testCase)
   % The Results body is curated: verdict counts, fill volume by family,
   % admission skill, and rule-stated example figures — not a bare list
   % of uncurated headline images.
   text = testCase.TestData.qmd;
   testCase.verifyTrue(contains(text, "## Cohort readiness"));
   testCase.verifyTrue(contains(text, ...
      "| policy_verdict | station_years | stations |"));
   testCase.verifyTrue(contains(text, "| ready | 1 | 1 |"));
   testCase.verifyTrue(contains(text, "| out_of_policy_window | 1 | 1 |"));
   % The relational-diagnostic story is stated with the verdict table.
   testCase.verifyTrue(contains(text, "Relational diagnostics"));
   testCase.verifyTrue(contains(text, "never verdict inputs"));

   testCase.verifyTrue(contains(text, "## Fill volume by method family"));
   testCase.verifyTrue(contains(text, "| site | observed_pct |"));
   testCase.verifyTrue(contains(text, "bounded_interp_pct"));
   testCase.verifyTrue(contains(text, "station_transfer_pct"));
   % Cohort-wide zero columns are pruned from the rendered view only.
   testCase.verifyFalse(contains(text, "racmo_pct"));

   testCase.verifyTrue(contains(text, ...
      "## Method admission and held-out skill"));
   testCase.verifyTrue(contains(text, "| channel | family | " + ...
      "admitted_strata | denied_candidates | median_improvement_pct |"));
   testCase.verifyTrue(contains(text, "| tair | donor | 1 | 0 | 50 |"));
   testCase.verifyTrue(contains(text, ...
      "| tair | bounded_interp | 1 | 0 | 25 |"));
   % A group with no admitted strata has no median; the shared markdown
   % writer renders the NaN as a blank cell.
   testCase.verifyTrue(contains(text, "| tair | mar | 0 | 1 |  |"));

   % Generic method-detail examples are delegated to the companion; only
   % the interpretation catalog's scientific figures remain in Results.
   testCase.verifyFalse(contains(text, "## Curated example figures"));
   % The machine-readable links include the new per-family CSV.
   testCase.verifyFalse(contains(text, ...
      "Preview of the per-site, per-channel summary"));
   testCase.verifyTrue(contains(text, "../qa/gapfill_fill_by_family.csv"));
end

function test_scientific_interpretation_catalog(testCase)
   % The declared category set is schema-stable: one fixture category is
   % material and rendered, while the others remain explicit absent rows.
   root = testCase.TestData.root;
   catalog = testCase.TestData.report.interpretation_catalog;
   expected = ["swd_fill_geometry"; "long_gap_limited_context"; ...
      "climatology_variability"; "boundary_behavior"; ...
      "hourly_disaggregation"; "residual_missing"; ...
      "native_sensor_flat_run"];
   testCase.verifyEqual(string(catalog.category), expected);
   testCase.verifyEqual(nnz(catalog.present), 1);
   selected = catalog.category == "long_gap_limited_context";
   testCase.verifyTrue(catalog.present(selected));
   testCase.verifyEqual(string(catalog.site(selected)), "tsta");
   testCase.verifyEqual(string(catalog.method(selected)), "donor:dsta");
   testCase.verifyTrue(isfile(fullfile(root, 'figures', ...
      string(catalog.figure(selected)))));
   absent = ~catalog.present;
   testCase.verifyTrue(all(contains( ...
      string(catalog.interpretation(absent)), "No material example")));

   csv = readtable(fullfile(root, 'qa', ...
      'gapfill_interpretation_catalog.csv'), 'Delimiter', ',', ...
      'TextType', 'string');
   testCase.verifyEqual(height(csv), height(catalog));
   text = testCase.TestData.qmd;
   testCase.verifyTrue(contains(text, ...
      "## Scientific interpretation and decision-sensitive outcomes"));
   testCase.verifyTrue(contains(text, "**Selection rule.**"));
   testCase.verifyTrue(contains(text, "**Window.**"));
   testCase.verifyTrue(contains(text, ...
      "| category | present | status | site | channel |"));
   testCase.verifyFalse(contains(text, ...
      "| category | present | status | site | channel | method |"));
   testCase.verifyTrue(contains(text, ...
      "../qa/gapfill_interpretation_catalog.csv"));
end

function test_detail_set_is_bounded_and_deterministic(testCase)
   % A one-figure station budget selects the longest leaders
   % deterministically without mixing the bounded and donor methods.
   root = string(tempname());
   mkdir(root);
   cleaner = onCleanup(@() rmdir(root, 's'));
   writeLayoutFixture(root);
   first = icemodel.verification.report.buildGapFillReport( ...
      sites="all", render=false, max_detail_figures_per_site=1, ...
      filled_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      fig_dir=fullfile(root, 'figures'), ...
      report_dir=fullfile(root, 'report'));
   is_detail = ~startsWith(first.figure_ledger.method, ...
      "interpretation:") & first.figure_ledger.method ~= "overview";
   testCase.verifyEqual(nnz(is_detail), 1);
   testCase.verifyEqual(string(first.figure_ledger.method(is_detail)), ...
      "donor:dsta");
   second = icemodel.verification.report.buildGapFillReport( ...
      sites="all", render=false, max_detail_figures_per_site=1, ...
      filled_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      fig_dir=fullfile(root, 'figures'), ...
      report_dir=fullfile(root, 'report'));
   testCase.verifyEqual(second.figure_ledger.file, ...
      first.figure_ledger.file);
   clear cleaner
end

function test_no_product_sites_stated(testCase)
   % Staged native stations with no filled product are a key result: the
   % section names them and states the independently verified refusal.
   text = testCase.TestData.qmd;
   testCase.verifyTrue(contains(text, ...
      "## Stations without filled products"));
   testCase.verifyTrue(contains(text, ...
      "A station appears in the no-product table only when no validated"));
   testCase.verifyTrue(contains(text, "kant"));
   testCase.verifyTrue(contains(text, "windowRecordDisjoint"));
   absent = testCase.TestData.report.absent_products;
   testCase.verifyEqual(string(absent.site), "kant");
   testCase.verifyEqual(string(absent.reason), "windowRecordDisjoint");
   testCase.verifyLessThan(absent.proxy_end, absent.native_start);
end

function test_gap_explanations_are_exhaustive(testCase)
   % The main Results section and CSVs state every residual run and every
   % non-ready station-year; the fixture has no missing product samples
   % but one explicit out-of-policy station-year.
   root = testCase.TestData.root;
   text = testCase.TestData.qmd;
   report = testCase.TestData.report;
   testCase.verifyTrue(contains(text, ...
      "## Why anything remains unfilled"));
   testCase.verifyEmpty(report.residual_gaps);
   testCase.verifyEqual(height(report.readiness_blockers), 1);
   testCase.verifyEqual(report.readiness_blockers.policy_verdict, ...
      "out_of_policy_window");
   testCase.verifyTrue(contains( ...
      report.readiness_blockers.plain_explanation, ...
      "outside the staged proxy window"));
   for filename = ["gapfill_residual_gaps.csv", ...
         "gapfill_readiness_blockers.csv", ...
         "gapfill_absent_products.csv"]
      testCase.verifyTrue(isfile(fullfile(root, 'qa', filename)));
      testCase.verifyTrue(contains(text, "../qa/" + filename));
   end
end

function test_nonempty_residual_gap_is_exact_and_explained(testCase)
   % Every shipped missing run must appear once with the producer's exact
   % final-tier denial; the report must not infer or summarize it away.
   root = string(tempname());
   mkdir(root);
   cleaner = onCleanup(@() rmdir(root, 's'));
   expected = writeResidualFixture(root, false);
   report = icemodel.verification.report.buildGapFillReport( ...
      sites="all", render=false, ...
      filled_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      fig_dir=fullfile(root, 'figures'), ...
      report_dir=fullfile(root, 'report'));
   row = report.residual_gaps(report.residual_gaps.site == "tsta", :);
   testCase.verifyEqual(height(row), 1);
   testCase.verifyEqual(row.channel, "swd");
   testCase.verifyEqual(row.year, 2020);
   testCase.verifyEqual(row.start_time, expected.start_time);
   testCase.verifyEqual(row.end_time, expected.end_time);
   testCase.verifyEqual(row.duration_hours, expected.duration_hours);
   testCase.verifyEqual(row.reason, string(expected.detail));
   clear cleaner
end

function test_residual_gap_mismatch_refuses_report(testCase)
   % A ledger that omits one shipped missing sample cannot certify a report.
   root = string(tempname());
   mkdir(root);
   cleaner = onCleanup(@() rmdir(root, 's'));
   writeResidualFixture(root, true);
   testCase.verifyError(@() ...
      icemodel.verification.report.buildGapFillReport( ...
      sites="all", render=false, ...
      filled_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      fig_dir=fullfile(root, 'figures'), ...
      report_dir=fullfile(root, 'report')), ...
      'icemodel:report:buildGapFillReport:residualGapMismatch');
   clear cleaner
end

function test_partial_final_day_overlaps_calendar_year(testCase)
   % An exact support window containing only late December 31 samples still
   % overlaps that calendar year; midnight at the day's start cannot exclude it.
   root = string(tempname());
   mkdir(root);
   cleaner = onCleanup(@() rmdir(root, 's'));
   writeLayoutFixture(root);
   manifest_file = fullfile(root, 'qa', 'plans', ...
      'tsta-report-inputs.json');
   manifest = jsondecode(fileread(manifest_file));
   manifest.acceptance_window.start = "2021-12-31 20:00:00";
   manifest.acceptance_window.end = "2021-12-31 23:45:00";
   fid = fopen(manifest_file, 'w');
   file_cleanup = onCleanup(@() fclose(fid));
   fprintf(fid, '%s', jsonencode(manifest));
   clear file_cleanup
   icemodel.verification.report.buildGapFillReport( ...
      sites="all", render=false, ...
      filled_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      fig_dir=fullfile(root, 'figures'), ...
      report_dir=fullfile(root, 'report'));
   readiness = readtable(fullfile(root, 'qa', ...
      'promice_filled_readiness.csv'), 'TextType', 'string');
   row = readiness(readiness.year == 2021, :);
   testCase.verifyEqual(row.policy_verdict, "not_forcing_ready");
   clear cleaner
end

function test_overlapping_absent_product_refuses_report(testCase)
   % An absent product with native/proxy overlap is unexplained; filename
   % absence alone must not mislabel it as windowRecordDisjoint.
   root = string(tempname());
   mkdir(root);
   cleaner = onCleanup(@() rmdir(root, 's'));
   writeLayoutFixture(root);
   disjoint = dir(fullfile(root, 'met', 'mar3.11', ...
      'met_kant_mar3.11_*_15m.mat'));
   delete(fullfile(disjoint.folder, disjoint.name));
   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC'):minutes(15): ...
      datetime(2020, 3, 31, 23, 45, 0, 'TimeZone', 'UTC')).';
   writeAbsentProxy(root, "kant", times);
   testCase.verifyError(@() ...
      icemodel.verification.report.buildGapFillReport( ...
      sites="all", render=false, ...
      filled_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      fig_dir=fullfile(root, 'figures'), ...
      report_dir=fullfile(root, 'report')), ...
      'icemodel:report:buildGapFillReport:unexplainedMissingProduct');
   clear cleaner
end

function test_new_disjoint_absence_is_stated(testCase)
   % Any independently proven disjoint station is reported without a
   % station-name registry, so newly staged coverage needs no code change.
   root = string(tempname());
   mkdir(root);
   cleaner = onCleanup(@() rmdir(root, 's'));
   writeLayoutFixture(root);
   source = dir(fullfile(root, 'met', 'promice', ...
      'met_kant_promice_*_15m.mat'));
   S = load(fullfile(source.folder, source.name), 'met');
   met = S.met;
   met.Properties.UserData.site = "nofl";
   save(fullfile(root, 'met', 'promice', ...
      'met_nofl_promice_20200101_20200331_15m.mat'), 'met');
   writeAbsentProxy(root, "nofl", ...
      (datetime(2019, 1, 1, 'TimeZone', 'UTC'):minutes(15): ...
      datetime(2019, 3, 31, 23, 45, 0, 'TimeZone', 'UTC')).');
   report = icemodel.verification.report.buildGapFillReport( ...
      sites="all", render=false, ...
      filled_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      fig_dir=fullfile(root, 'figures'), ...
      report_dir=fullfile(root, 'report'));
   row = report.absent_products(report.absent_products.site == "nofl", :);
   testCase.verifyEqual(height(row), 1);
   testCase.verifyEqual(row.reason, "windowRecordDisjoint");
   testCase.verifyLessThan(row.proxy_end, row.native_start);
   clear cleaner
end

function test_designed_refusal_without_proxy_window_is_stated(testCase)
   % A vetted station with no validated proxy has no A6 product span. The
   % report must state that exact reason instead of claiming disjoint spans.
   root = string(tempname());
   mkdir(root);
   cleaner = onCleanup(@() rmdir(root, 's'));
   writeLayoutFixture(root);
   proxy = dir(fullfile(root, 'met', 'mar3.11', ...
      'met_kant_mar3.11_*_15m.mat'));
   delete(fullfile(proxy.folder, proxy.name));
   mkdir(fullfile(root, 'figures', 'overview'));
   mkdir(fullfile(root, 'figures', 'detail'));
   stale = [ ...
      string(fullfile(root, 'figures', 'overview', ...
      'kant_overview.png')); ...
      string(fullfile(root, 'figures', 'detail', ...
      'kant_rh.png')); ...
      string(fullfile(root, 'figures', 'kant_old.png'))];
   for file = stale.'
      fid = fopen(file, 'w');
      fprintf(fid, 'stale');
      fclose(fid);
   end
   report = icemodel.verification.report.buildGapFillReport( ...
      sites="all", render=false, ...
      filled_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      fig_dir=fullfile(root, 'figures'), ...
      report_dir=fullfile(root, 'report'));
   row = report.absent_products(report.absent_products.site == "kant", :);
   testCase.verifyEqual(row.reason, "noProxyWindow");
   testCase.verifyTrue(isnat(row.proxy_start) && isnat(row.proxy_end));
   testCase.verifyTrue(contains(row.explanation, ...
      "no fillable product period", 'IgnoreCase', true));
   testCase.verifyFalse(any(arrayfun(@isfile, stale)));
   clear cleaner
end

function test_fill_by_family_from_provenance(testCase)
   % The per-station family percentages derive from the shipped
   % provenance channels: 48 bounded_interp and 96 station_transfer
   % samples over two science channels of the fixture record.
   root = testCase.TestData.root;
   fill_by_family = testCase.TestData.report.fill_by_family;
   testCase.verifyEqual(string(fill_by_family.site), "tsta");
   n_total = 2 * 8736;
   returned = fill_by_family.bounded_interp_pct;
   expected = round(100 * 48 / n_total, 3);
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12);
   returned = fill_by_family.station_transfer_pct;
   expected = round(100 * 96 / n_total, 3);
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12);
   returned = fill_by_family.observed_pct;
   expected = round(100 * (n_total - 144) / n_total, 3);
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12);
   testCase.verifyEqual(fill_by_family.missing_pct, 0);
   % The CSV ships the full registry schema even for pruned families.
   csv = readtable(fullfile(root, 'qa', 'gapfill_fill_by_family.csv'), ...
      'Delimiter', ',', 'TextType', 'string');
   testCase.verifyTrue(all(ismember(["site", "observed_pct", ...
      "bounded_interp_pct", "station_transfer_pct", "mar_pct", ...
      "racmo_pct", "missing_pct"], ...
      string(csv.Properties.VariableNames))));
end

%% method-accent layer split

function test_method_accent_layers_only_own_method(testCase)
   % The panel color logic: own-method fills get the accent color, other
   % methods' fills the muted context color, observed samples the dark
   % color — keyed by provenance plus the audit's own-method mask.
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   values = [1; 2; 3; 4; 5; NaN];
   provenance = double([codes.observed; codes.bounded_interp; ...
      codes.mar; codes.raw_shortwave; codes.missing; ...
      codes.bounded_interp]);
   own = logical([0; 1; 0; 0; 0; 1]);
   layers = icemodel.verification.report.methodFillLayers( ...
      values, provenance, own);
   % Observed layer: true observations plus raw-shortwave measurements.
   testCase.verifyEqual(layers.observed, [1; NaN; NaN; 4; NaN; NaN]);
   % Own layer: only the own-mask fill with a finite value.
   testCase.verifyEqual(layers.own_fill, [NaN; 2; NaN; NaN; NaN; NaN]);
   % Other layer: the foreign (MAR) fill only — never the missing
   % sentinel, never the own-method sample.
   testCase.verifyEqual(layers.other_fill, [NaN; NaN; 3; NaN; NaN; NaN]);
   % Colors come from the single style registry and stay distinct, so a
   % bounded_interp panel can never render a MAR fill in its accent.
   style = icemodel.verification.report.gapfillFigureStyle();
   testCase.verifyEqual(layers.own_color, style.accent);
   testCase.verifyEqual(layers.other_color, style.context);
   testCase.verifyEqual(layers.observed_color, style.observed);
   testCase.verifyNotEqual(layers.own_color, layers.other_color);
   testCase.verifyNotEqual(layers.own_color, layers.observed_color);
end

%% rerender

function test_rerender_is_idempotent(testCase)
   % A second build over the published tree reconciles the same ledger
   % and leaves the split layout and foreign figures intact. This test
   % mutates the shared tree consistently, so it must run last.
   root = testCase.TestData.root;
   returned = icemodel.verification.report.buildGapFillReport( ...
      sites="all", render=false, ...
      filled_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      fig_dir=fullfile(root, 'figures'), ...
      report_dir=fullfile(root, 'report'));
   expected = testCase.TestData.report.figure_ledger.file;
   testCase.verifyEqual(returned.figure_ledger.file, expected);
   testCase.verifyTrue(isfile(fullfile(root, 'figures', 'overview', ...
      'tsta_overview.png')));
   testCase.verifyTrue(isfile(fullfile(root, 'figures', 'overview', ...
      'ghost_overview.png')));
   testCase.verifyTrue(isfile(fullfile(root, 'figures', ...
      'ghost_legacyflat.png')));
end

%% fixture helpers

function writeLayoutFixture(root)
   %WRITELAYOUTFIXTURE Hand-build the saved artifacts the generator reads.
   % Layout: <root>/met/promice (native), <root>/met/promice_filled
   % (product), <root>/qa (plans + ledger). selectedDataRoot resolves the
   % compact <root>/met layout to <root>, so manifest paths are relative
   % to <root>.
   mkdir(fullfile(root, 'met', 'promice'));
   mkdir(fullfile(root, 'met', 'promice_filled'));
   mkdir(fullfile(root, 'met', 'mar3.11'));
   mkdir(fullfile(root, 'qa', 'plans'));
   mkdir(fullfile(root, 'qa', 'ledger'));
   mkdir(fullfile(root, 'report'));

   % One quarter of 2020 at 15 minutes: 8736 samples, small enough for a
   % fast render yet long enough for windowed context on both gaps.
   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC'): minutes(15): ...
      datetime(2020, 3, 31, 23, 45, 0, 'TimeZone', 'UTC')).';
   n = numel(times);
   codes = icemodel.forcing.reconstruct.provenanceCodes();

   % Two audited fills in tair: a 12 h bounded_interp gap and a 24 h
   % donor gap two days later, so each panel's context window contains
   % the OTHER method's fill — the exact POLICY D-31 accent scenario.
   i1 = find(times == datetime(2020, 1, 20, 6, 0, 0, 'TimeZone', 'UTC'));
   gap1 = i1:(i1 + 47);
   i2 = find(times == datetime(2020, 1, 22, 0, 0, 0, 'TimeZone', 'UTC'));
   gap2 = i2:(i2 + 95);

   % Smooth diurnal series keep the figures physically plausible without
   % exercising any science code path.
   tair = 260 + 5 * sin(2 * pi * (0:n - 1).' / 96);
   swd = max(0, 300 * sin(2 * pi * (0:n - 1).' / 96));
   tair_provenance = repmat(double(codes.observed), n, 1);
   tair_provenance(gap1) = double(codes.bounded_interp);
   tair_provenance(gap2) = double(codes.station_transfer);
   swd_provenance = repmat(double(codes.observed), n, 1);

   % Filled product: finite everywhere, provenance-paired channels, and
   % the identity metadata loadMet verifies.
   met = timetable(times, tair, swd, tair_provenance, swd_provenance, ...
      'VariableNames', {'tair', 'swd', 'tair_provenance', ...
      'swd_provenance'});
   met.Properties.UserData = struct('site', "tsta", ...
      'gapfill_product', "promice_filled");
   filled_file = fullfile(root, 'met', 'promice_filled', ...
      'met_tsta_promice_filled_20200101_20200331_15m.mat');
   save(filled_file, 'met');

   % Native counterpart: the audited gap samples are missing at source.
   native_tair = tair;
   native_tair([gap1, gap2]) = NaN;
   met = timetable(times, native_tair, swd, 'VariableNames', ...
      {'tair', 'swd'});
   met.Properties.UserData = struct('site', "tsta");
   native_file = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20200331_15m.mat');
   save(native_file, 'met');

   % A native-only station with NO filled product: the Results section
   % must surface it as a windowRecordDisjoint refusal.
   refusal_location = struct('site', "kant", 'lat', 60, 'lon', -45, ...
      'elev', 500);
   met = timetable(times(1:96), tair(1:96), 'VariableNames', {'tair'});
   met.Properties.UserData = refusal_location;
   save(fullfile(root, 'met', 'promice', ...
      'met_kant_promice_20200101_20200331_15m.mat'), 'met');
   writeAbsentProxy(root, "kant", ...
      (datetime(2019, 1, 1, 'TimeZone', 'UTC'):minutes(15): ...
      datetime(2019, 3, 31, 23, 45, 0, 'TimeZone', 'UTC')).');

   % Plan + audit artifact: the audit maps the two methods to their
   % sample spans; the plan carries producer-pinned options and the
   % admission evidence the diagnostics/Results tables aggregate.
   audit_record = cell2table({ ...
      'tair', times(gap1(1)), times(gap1(end)), 12, 'bounded_interp', ...
      '48 samples; bounded linear interpolation', 'bounded_interp'; ...
      'tair', times(gap2(1)), times(gap2(end)), 24, 'donor:dsta', ...
      '96 samples; donor transfer', 'donor:dsta'}, ...
      'VariableNames', {'channel', 'start_time', 'end_time', ...
      'duration_hours', 'method', 'detail', 'context_id'});
   plan_record = layoutPlanRecord();
   plan_file = fullfile(root, 'qa', 'plans', 'tsta-plan.mat');
   save(plan_file, 'plan_record', 'audit_record');

   % Readiness ledger: one ready year inside the acceptance window and
   % one not-ready year outside it, so the policy view produces two
   % distinct verdict rows for the cohort table.
   readiness = table(["tsta"; "tsta"], [2020; 2021], [0.95; 0.90], ...
      [0; 0.1], ["ready"; "not_forcing_ready"], ...
      [""; "tair 10.0% missing"], [0; 0.1], ...
      ["ready"; "not_forcing_ready"], [""; "tair 10.0% missing"], ...
      [1; 1], 'VariableNames', {'site', 'year', ...
      'native_core_coverage', 'worst_residual_invalid', ...
      'verdict_icemodel', 'reason_icemodel', ...
      'worst_residual_snowmodel', 'verdict_snowmodel', ...
      'reason_snowmodel', 'n_admitted_methods'});
   readiness_file = fullfile(root, 'qa', 'ledger', 'tsta-readiness.csv');
   writetable(readiness, readiness_file);

   % Producer manifest: every consumed artifact byte-pinned, paths
   % relative to the selected data root, plus the acceptance window.
   artifacts = [ ...
      manifestArtifact("native", root, native_file); ...
      manifestArtifact("filled", root, filled_file); ...
      manifestArtifact("plan", root, plan_file); ...
      manifestArtifact("readiness", root, readiness_file)];
   manifest = struct('site', "tsta", ...
      'path_base', "selected_data_root", 'artifacts', artifacts, ...
      'acceptance_window', struct('start', "2020-01-01 00:00:00", ...
      'end', "2020-03-31 23:45:00"));
   fid = fopen(fullfile(root, 'qa', 'plans', ...
      'tsta-report-inputs.json'), 'w');
   cleaner = onCleanup(@() fclose(fid));
   fprintf(fid, '%s', jsonencode(manifest));
   clear cleaner
end

function expected = writeResidualFixture(root, mismatch)
   %WRITERESIDUALFIXTURE Add one real gap and its reconciled denial row.
   writeLayoutFixture(root);
   filled_file = fullfile(root, 'met', 'promice_filled', ...
      'met_tsta_promice_filled_20200101_20200331_15m.mat');
   S = load(filled_file, 'met');
   met = S.met;
   gap = 500:503;
   met.swd(gap) = NaN;
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   met.swd_provenance(gap) = double(codes.missing);
   save(filled_file, 'met');

   % Build the producer row with the same canonical segment helper used by
   % production. The mismatch case deliberately omits the last sample.
   missing = false(height(met), 1);
   ledger_gap = gap;
   if mismatch
      ledger_gap = gap(1:end - 1);
   end
   missing(ledger_gap) = true;
   detail = "no admissible final-tier candidate in fixture";
   segments = icemodel.forcing.reconstruct.auditSegments( ...
      met.Properties.RowTimes, missing, "swd", "unfilled", detail);
   plan_file = fullfile(root, 'qa', 'plans', 'tsta-plan.mat');
   P = load(plan_file, 'plan_record', 'audit_record');
   residual = cell2table(vertcat(segments{:}), ...
      'VariableNames', P.audit_record.Properties.VariableNames);
   audit_record = [P.audit_record; residual];
   plan_record = P.plan_record;
   save(plan_file, 'plan_record', 'audit_record');

   % Refresh only the fixture's byte pins after changing its product and
   % plan; production uses the same manifest verification at report entry.
   manifest_file = fullfile(root, 'qa', 'plans', ...
      'tsta-report-inputs.json');
   manifest = jsondecode(fileread(manifest_file));
   for k = 1:numel(manifest.artifacts)
      filename = fullfile(root, string(manifest.artifacts(k).path));
      info = dir(filename);
      manifest.artifacts(k).bytes = info.bytes;
      manifest.artifacts(k).sha256 = ...
         icemodel.verification.setup.fileSha256(filename);
   end
   fid = fopen(manifest_file, 'w');
   file_cleanup = onCleanup(@() fclose(fid));
   fprintf(fid, '%s', jsonencode(manifest));
   clear file_cleanup
   expected = residual(1, :);
end

function record = manifestArtifact(role, root, filename)
   %MANIFESTARTIFACT Pin one fixture artifact by role, bytes, and SHA-256.
   info = dir(filename);
   record = struct('role', role, ...
      'path', extractAfter(string(filename), root + filesep), ...
      'bytes', info.bytes, ...
      'sha256', icemodel.verification.setup.fileSha256(filename));
end

function plan_record = layoutPlanRecord()
   %LAYOUTPLANRECORD Minimal schema-correct plan for the report reader.
   % Shapes mirror stationMethodPlan's persisted plan: per-channel
   % methods with selection/evaluation strata and a denials table. The
   % numbers are chosen so the admission aggregate is exactly checkable:
   % donor 1.2/2.4 -> 50% improvement, bounded_interp 1.5/2.0 -> 25%.
   defaults = icemodel.forcing.reconstruct.setopts();
   ropts = struct('rmse_improvement', 0.17, 'cap_hours', 6, ...
      'min_variability_ratio', 0.63, 'max_variability_ratio', 1.37, ...
      'min_coverage', 0.82, 'blend_hours', 24, 'donor_sites', "dsta", ...
      'use_ktransect', false, 'use_gcnet', false, ...
      'proxy_catalog', defaults.proxy_catalog);
   donor = planMethod("donor:dsta", uint8(2), 1.2, 2.4);
   interp = planMethod("bounded_interp", uint8(1), 1.5, 2.0);
   methods = [donor; interp];
   denials = table({'mar'}, 24, {'all'}, ...
      {'insufficient staged overlap'}, 'VariableNames', ...
      {'candidate', 'bucket', 'season', 'reasons'});
   tair_channel = struct('channel', "tair", 'methods', {methods}, ...
      'denials', {denials});
   % swd has no admissions or denials; empty-but-typed members keep the
   % struct-array fields consistent.
   swd_channel = struct('channel', "swd", ...
      'methods', {methods([], 1)}, 'denials', {denials([], :)});
   plan_record = struct('station', "tsta", ...
      'reconstruction_options', ropts, ...
      'channels', [tair_channel; swd_channel]);
end

function writeAbsentProxy(root, site, times)
   %WRITEABSENTPROXY Write one identity-valid staged MAR proxy window.
   % The caller controls its dates so tests cover disjoint and overlapping
   % absent-product evidence with the exact production validator.
   tair = repmat(260, numel(times), 1);
   met = timetable(times, tair, 'VariableNames', {'tair'});
   met.Properties.UserData = struct('site', site, ...
      'lat_wgs84', 60, 'lon_wgs84', -45, ...
      'mar_qc_status', "not_applicable", ...
      'source_files', "MARv3.11-fixture.nc");
   filename = sprintf('met_%s_mar3.11_%s_%s_15m.mat', site, ...
      string(times(1), 'yyyyMMdd'), string(times(end), 'yyyyMMdd'));
   save(fullfile(root, 'met', 'mar3.11', filename), 'met');
end

function method = planMethod(name, code, rmse, baseline_rmse)
   %PLANMETHOD One admitted plan method with a single 24 h "all" stratum.
   selection = metricStratumRow(rmse);
   selection.baseline = "persistence";
   selection.baseline_rmse = baseline_rmse;
   evaluation = struct('bucket', 24, ...
      'by_stratum', metricStratumRow(rmse));
   method = struct('name', name, 'code', code, 'seasons', "all", ...
      'buckets', 24, 'max_validated_hours', 48, ...
      'uncertainty', "empirical", 'selection', {selection}, ...
      'evaluation', {evaluation});
end

function row = metricStratumRow(rmse)
   %METRICSTRATUMROW One bucket/season metric row in the persisted schema.
   row = table(24, "all", 20, 1, 0.1, rmse, 0.95, 1.0, 0.4, 0, 0, ...
      0.7, 0.96, 1, 'VariableNames', {'bucket', 'season', 'n', ...
      'coverage', 'bias', 'rmse', 'correlation', 'variability_ratio', ...
      'within_gap_observed_spread', 'bound_violations', ...
      'boundary_jump_rate', 'sigma1_coverage', 'sigma2_coverage', ...
      'provenance_accounting'});
end

function seedStaleFigures(root)
   %SEEDSTALEFIGURES Plant stale and foreign PNGs the install must sort.
   % Selected-site leftovers (flat legacy layout and a stale detail PNG)
   % must vanish; foreign-station figures must survive everywhere.
   mkdir(fullfile(root, 'figures', 'overview'));
   mkdir(fullfile(root, 'figures', 'detail'));
   placeholders = [ ...
      fullfile(root, 'figures', 'tsta_tair_legacyflat_gapfill.png'); ...
      fullfile(root, 'figures', 'ghost_legacyflat.png'); ...
      fullfile(root, 'figures', 'detail', 'tsta_tair_stale_gapfill.png'); ...
      fullfile(root, 'figures', 'overview', 'ghost_overview.png')];
   for k = 1:numel(placeholders)
      fid = fopen(placeholders(k), 'w');
      fprintf(fid, 'stale');
      fclose(fid);
   end
end
