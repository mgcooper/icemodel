function tests = test_report_helpers
   %TEST_REPORT_HELPERS Cover the shared report helpers.
   %
   % Three groups. Both report builders render saved metadata through the
   % markdown helpers. A value containing markup, HTML, or backticks must not
   % be able to change the document structure. The metric helpers score what
   % the reports show. The provenance helpers decide whether a saved cohort
   % still supports a report. They are the content digests, the
   % report-consumed channel list, the schema compatibility test, and the
   % physics fingerprint.
   %
   % test_physics_fingerprint_is_stable_and_covers_physics resolves the model
   % defaults through icemodel.setopts, so that one case needs a provisioned
   % workspace. The rest read only namelists and in-memory values.
   tests = functiontests(localfunctions);
end

function test_code_span_fence_grows_past_backticks(testCase)
   % A value containing backticks must not be able to close the span early.
   returned = icemodel.verification.report.markdownCode("a``b");

   testCase.verifyTrue(startsWith(returned, "```"))
   testCase.verifyTrue(endsWith(returned, "```"))
   testCase.verifyTrue(contains(returned, "a``b"))
end

function test_code_span_pads_when_the_value_touches_a_backtick(testCase)
   % A leading or trailing backtick needs a space, or the fence merges with it.
   returned = icemodel.verification.report.markdownCode("`v1`");

   testCase.verifyTrue(contains(returned, " `v1` "))
end

function test_code_span_neutralizes_html_delimiters(testCase)
   % No rendered report may contain a raw script sequence, even inside a span.
   returned = icemodel.verification.report.markdownCode("<script>x</script>");

   testCase.verifyFalse(contains(returned, "<script>"))
   testCase.verifyTrue(contains(returned, "&lt;script&gt;"))
end

function test_control_characters_are_collapsed(testCase)
   % A newline inside saved metadata would otherwise break the line it sits on.
   returned = icemodel.verification.report.markdownCode( ...
      sprintf("a\r\nb"));

   testCase.verifyFalse(contains(returned, newline))
   testCase.verifyTrue(contains(returned, "a b"))
end

function test_escaped_text_cannot_introduce_markup(testCase)
   % Every ASCII punctuation character is escaped, backslash first.
   returned = icemodel.verification.report.escapeMarkdownText("a**b**|c");

   testCase.verifyEqual(returned, "a\*\*b\*\*\|c")
   testCase.verifyEqual( ...
      icemodel.verification.report.escapeMarkdownText("a\b"), "a\\b")
end

function test_markdown_table_emits_a_header_rule(testCase)
   % A table without its rule renders as a paragraph.
   returned = icemodel.verification.report.markdownTable( ...
      table("a", 1, 'VariableNames', {'name', 'value'}));

   testCase.verifyTrue(any(contains(returned, "---")))
   testCase.verifyTrue(any(contains(returned, "name")))
end

function test_single_pair_still_reports_an_exact_error(testCase)
   % One pair defines no spread, so NSE is undefined, but the error itself is
   % known exactly and must not be thrown away.
   returned = icemodel.verification.helpers.residualMetrics([3 NaN], [1 2]);

   testCase.verifyEqual(returned.n_pairs, 1)
   testCase.verifyEqual(returned.bias, 2, AbsTol=1e-12)
   testCase.verifyEqual(returned.rmse, 2, AbsTol=1e-12)
   testCase.verifyEqual(returned.mae, 2, AbsTol=1e-12)
   testCase.verifyTrue(isnan(returned.nse))
end

function test_no_finite_pairs_reports_nothing(testCase)
   % With no overlap there is no error to report.
   returned = icemodel.verification.helpers.residualMetrics([NaN NaN], [1 2]);

   testCase.verifyEqual(returned.n_pairs, 0)
   testCase.verifyTrue(isnan(returned.rmse))
   testCase.verifyTrue(isnan(returned.bias))
end

function test_residual_metrics_values(testCase)
   % Bias, MAE, RMSE, and max error on a known residual.
   returned = icemodel.verification.helpers.residualMetrics([2 4 6], [1 2 3]);

   testCase.verifyEqual(returned.bias, 2, AbsTol=1e-12)
   testCase.verifyEqual(returned.mae, 2, AbsTol=1e-12)
   testCase.verifyEqual(returned.max_abs_error, 3, AbsTol=1e-12)
   testCase.verifyEqual(returned.n_pairs, 3)
end

function test_residual_metrics_nse_is_undefined_without_variance(testCase)
   % A constant observed series has nothing to explain, so NSE is NaN rather
   % than infinite.
   returned = icemodel.verification.helpers.residualMetrics([1 2 3], [2 2 2]);

   testCase.verifyTrue(isnan(returned.nse))
   testCase.verifyFalse(isnan(returned.rmse))
end

function test_residual_metrics_perfect_model_scores_one(testCase)
   % A model equal to the observations explains all of the variance.
   returned = icemodel.verification.helpers.residualMetrics( ...
      [1 2 3], [1 2 3]);

   testCase.verifyEqual(returned.nse, 1, AbsTol=1e-12)
   testCase.verifyEqual(returned.rmse, 0, AbsTol=1e-12)
end

function test_escaped_html_cannot_form_a_tag(testCase)
   % Saved metadata reaches Markdown table cells through this helper, so an
   % angle bracket must not survive as the start of a tag.
   returned = icemodel.verification.report.escapeMarkdownText( ...
      "<script>alert(1)</script>");

   testCase.verifyFalse(contains(returned, "<script"))
   testCase.verifyTrue(contains(returned, "script"))
end

function test_sample_quantile_interpolates_between_samples(testCase)
   % The percentile gates that admit a site-year run through this, so the
   % interpolation and the empty guard are worth pinning.
   returned = icemodel.verification.helpers.sampleQuantile([1 2 3 4], 0.5);
   testCase.verifyEqual(returned, 2.5, AbsTol=1e-12)

   % Endpoints land on the samples themselves.
   testCase.verifyEqual( ...
      icemodel.verification.helpers.sampleQuantile([1 2 3 4], 0), 1)
   testCase.verifyEqual( ...
      icemodel.verification.helpers.sampleQuantile([1 2 3 4], 1), 4)
end

function test_sample_quantile_sorts_and_handles_empty(testCase)
   % Callers must not have to sort first, and an empty sample has no quantile.
   testCase.verifyEqual( ...
      icemodel.verification.helpers.sampleQuantile([4 1 3 2], 0.5), 2.5, ...
      AbsTol=1e-12)
   testCase.verifyTrue(isnan( ...
      icemodel.verification.helpers.sampleQuantile([], 0.5)))
end

function test_format_value_renders_each_type(testCase)
   % Every Markdown table cell goes through this.
   testCase.verifyEqual( ...
      icemodel.verification.report.formatValue(NaN), "NA")
   testCase.verifyEqual( ...
      icemodel.verification.report.formatValue(true), "true")
   testCase.verifyTrue(contains( ...
      icemodel.verification.report.formatValue(2.5), "2"))

   % A multi-column table variable arrives as a vector.
   testCase.verifyTrue(contains( ...
      icemodel.verification.report.formatValue([1 2]), "1"))
end

function test_safe_label_collapses_control_characters(testCase)
   % A newline in an axis label breaks the figure.
   returned = icemodel.verification.report.safeLabel(sprintf("a\tb\nc"));

   testCase.verifyFalse(contains(returned, newline))
   testCase.verifyEqual(returned, "a b c")
end

function test_evaluation_season_bounds_are_utc(testCase)
   % Readiness and the runner must resolve the same window.
   policy = icemodel.verification.namelists.promiceAblationPolicy();
   [first, last] = ...
      icemodel.verification.helpers.evaluationSeason(2019, policy);

   testCase.verifyEqual(first.TimeZone, 'UTC')
   testCase.verifyEqual(year(first), 2019)
   testCase.verifyLessThan(first, last)
end

function test_ablation_ledger_increments_signs(testCase)
   % Sign conventions: loss of solid is positive, deposition is negative.
   ledger = struct( ...
      'mass_budget_phase_solid_mwe', -0.3, ...
      'mass_budget_vapor_solid_mwe', -0.1, ...
      'mass_budget_top_export_solid_mwe', 0.2, ...
      'mass_budget_top_export_liquid_mwe', 0.05);
   returned = ...
      icemodel.verification.helpers.ablationLedgerIncrements(ledger);

   testCase.verifyEqual(returned.solid_balance, 0.4, AbsTol=1e-12)
   testCase.verifyEqual(returned.surface_loss, 0.25, AbsTol=1e-12)
   testCase.verifyEqual(returned.solid_vapor_loss, 0.1, AbsTol=1e-12)
end

function test_sha256_digests_match_the_published_vectors(testCase)
   % One hashing implementation serves the file, byte, and text digests, so a
   % published vector is the check that keeps all three honest.

   empty_digest = ...
      "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855";
   abc_digest = ...
      "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad";

   testCase.verifyEqual( ...
      icemodel.verification.setup.bytesSha256(uint8([])), empty_digest)
   testCase.verifyEqual( ...
      icemodel.verification.setup.textSha256(""), empty_digest)
   testCase.verifyEqual( ...
      icemodel.verification.setup.textSha256("abc"), abc_digest)

   % A row and a column of the same bytes are the same content.
   bytes = uint8('abc');
   testCase.verifyEqual( ...
      icemodel.verification.setup.bytesSha256(bytes), abc_digest)
   testCase.verifyEqual( ...
      icemodel.verification.setup.bytesSha256(bytes'), abc_digest)

   % A file of those bytes must hash to the same value, which is what makes
   % the in-memory digest comparable to a stored artifact digest.
   fixture = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   scratch = fullfile(fixture.Folder, 'abc.txt');
   fid = fopen(scratch, 'w');
   testCase.assertGreaterThan(fid, 0)
   fwrite(fid, bytes);
   fclose(fid);
   testCase.verifyEqual( ...
      icemodel.verification.setup.fileSha256(scratch), abc_digest)
end

function test_ablation_report_channels_partition_by_report_table(testCase)
   % The gate compares the whole list, and each report table checks only the
   % channels it reads, so the groups must together be the whole list.

   components = ...
      icemodel.verification.namelists.ablationReportChannels('components');
   grid = icemodel.verification.namelists.ablationReportChannels('grid');
   all_fields = icemodel.verification.namelists.ablationReportChannels();

   testCase.verifyEqual(all_fields, [components, grid])
   testCase.verifyEqual(all_fields, ...
      icemodel.verification.namelists.ablationReportChannels('all'))
   testCase.verifyEqual(numel(unique(all_fields)), numel(all_fields))
   testCase.verifyError( ...
      @() icemodel.verification.namelists.ablationReportChannels('other'), ...
      'icemodel:verification:namelists:ablationReportChannels:kind')

   % Every channel the report reads must be one the runner records, or the
   % schema gate would reject every cohort the runner can produce.
   policy = icemodel.verification.namelists.promiceAblationPolicy();
   testCase.verifyTrue( ...
      all(ismember(all_fields, string(policy.required_model_fields))))
end

function test_model_schema_compatibility_is_one_directional(testCase)
   % The saved cohort must carry every channel the report reads. Extra
   % current channels are additive and only reported as unavailable.

   consumed = icemodel.verification.namelists.ablationReportChannels();
   policy = icemodel.verification.namelists.promiceAblationPolicy();
   current = string(policy.required_model_fields);

   % A saved schema one non-consumed channel short is a cohort that ran
   % before that channel existed. It stays usable, and the missing channel is
   % named as unavailable.
   droppable = setdiff(current, consumed, 'stable');
   testCase.assertNotEmpty(droppable)
   saved = setdiff(current, droppable(1), 'stable');
   returned = ...
      icemodel.verification.helpers.validateAblationModelSchema( ...
      saved, current);
   testCase.verifyEqual(returned, droppable(1))

   % An identical pair leaves nothing unavailable.
   testCase.verifyEmpty( ...
      icemodel.verification.helpers.validateAblationModelSchema( ...
      current, current))

   % A saved cohort missing a channel the report reads cannot be reported on.
   testCase.verifyError(@() ...
      icemodel.verification.helpers.validateAblationModelSchema( ...
      setdiff(current, consumed(1), 'stable'), current), ...
      'icemodel:verification:report:incompatibleModelSchema')

   % A channel the current namelists do not define is a removal or a
   % redefinition, so the saved values do not mean what the report says.
   testCase.verifyError(@() ...
      icemodel.verification.helpers.validateAblationModelSchema( ...
      current, setdiff(current, consumed(1), 'stable')), ...
      'icemodel:verification:report:incompatibleModelSchema')

   % When BOTH lists lack the channel, the advice must not be to rerun: a
   % rerun records the same reduced schema and the owner loops. Both branches
   % raise one identifier, so the message text is what distinguishes them.
   reduced = setdiff(current, consumed(1), 'stable');
   try
      icemodel.verification.helpers.validateAblationModelSchema( ...
         reduced, reduced);
      testCase.verifyFail('a schema missing a consumed channel must raise')
   catch err
      testCase.verifyEqual(err.identifier, ...
         'icemodel:verification:report:incompatibleModelSchema')
      testCase.verifySubstring(err.message, 'Restore the channel')
      testCase.verifyFalse(contains(err.message, 'Rerun the cohort'))
   end

   % The strict policy comparison excludes this field. A saved schema of the
   % wrong type must therefore raise invalidAblationPolicy, not a bare
   % conversion error.
   malformed = {struct('a', 1), 42, {"ok", 7}, [current, ""]};
   for k = 1:numel(malformed)
      testCase.verifyError(@() ...
         icemodel.verification.helpers.validateAblationModelSchema( ...
         malformed{k}, current), ...
         'icemodel:verification:report:invalidAblationPolicy')
   end
end

function test_physics_stamp_shape_test_rejects_every_damaged_shape(testCase)
   % The report asks this helper whether a saved stamp is usable, and never
   % raises. Every rejected shape is checked here rather than through a report
   % build, so the coverage costs no rendering.

   good = struct('opts_sha256', string(repmat('a', 1, 64)), ...
      'icemodel_version', "1.2.3", 'excluded_fields', "pathinput");
   testCase.verifyTrue( ...
      icemodel.verification.helpers.isPhysicsFingerprint(good))

   % A char digest is text too, so it must pass.
   char_digest = good;
   char_digest.opts_sha256 = repmat('a', 1, 64);
   testCase.verifyTrue( ...
      icemodel.verification.helpers.isPhysicsFingerprint(char_digest))

   % Not a struct, not scalar, a missing field, then each bad value shape on
   % each field. Reading a field of the first three would raise, so this also
   % pins the short-circuit order.
   wrong_container = {42, "text", {good}, [good, good], struct([]), ...
      rmfield(good, 'opts_sha256'), rmfield(good, 'icemodel_version')};
   bad_values = {[], "", '', ["a", "b"], 42, {"a"}, string(missing)};
   n_fixed = numel(wrong_container);
   damaged = cell(1, n_fixed + 2 * numel(bad_values));
   damaged(1:n_fixed) = wrong_container;
   for k = 1:numel(bad_values)
      bad_digest = good;
      bad_digest.opts_sha256 = bad_values{k};
      bad_version = good;
      bad_version.icemodel_version = bad_values{k};
      damaged{n_fixed + 2 * k - 1} = bad_digest;
      damaged{n_fixed + 2 * k} = bad_version;
   end
   for k = 1:numel(damaged)
      testCase.verifyFalse( ...
         icemodel.verification.helpers.isPhysicsFingerprint(damaged{k}), ...
         sprintf('shape %d must be rejected', k))
   end
end

function test_physics_fingerprint_is_stable_and_covers_physics(testCase)
   % The digest must repeat for unchanged code, or every report would warn.
   returned = icemodel.verification.helpers.physicsFingerprint();
   testCase.verifyEqual( ...
      icemodel.verification.helpers.physicsFingerprint(), returned)
   testCase.verifyEqual(strlength(returned.opts_sha256), 64)
   testCase.verifyEqual(returned.icemodel_version, ...
      string(icemodel.internal.version()))

   % The exclusions must drop the machine-dependent workspace paths, or the
   % same code would fingerprint differently on two computers.
   testCase.verifyTrue(all(ismember( ...
      ["pathinput", "pathoutput"], returned.excluded_fields)))

   % A physics option must reach the digest, and a workspace path must not.
   % Without both checks the gate could warn on every machine, or never warn
   % on a real physics change.
   defaults = icemodel.setopts("icemodel", "kanm", 2016, "kanm");
   moved = defaults;
   moved.f_ice_min = defaults.f_ice_min + 0.01;
   relocated = defaults;
   relocated.pathoutput = '/somewhere/else';

   baseline = icemodel.verification.helpers.physicsFingerprint(defaults);
   testCase.verifyNotEqual( ...
      icemodel.verification.helpers.physicsFingerprint(moved).opts_sha256, ...
      baseline.opts_sha256)
   testCase.verifyEqual( ...
      icemodel.verification.helpers.physicsFingerprint( ...
      relocated).opts_sha256, baseline.opts_sha256)

   % An appended output channel must leave the digest alone. configureRun
   % expands output_profile into vars1 and vars2, and warning about a channel
   % that changes no computed value is the false alarm this gate removes.
   channels = defaults;
   channels.vars1 = [defaults.vars1, {'a_new_diagnostic_channel'}];
   channels.vars2 = [defaults.vars2, {'another_new_channel'}];
   testCase.verifyEqual( ...
      icemodel.verification.helpers.physicsFingerprint( ...
      channels).opts_sha256, baseline.opts_sha256)

   % The instrument heights come from the reference station, so the digest
   % must not depend on which station resolved the defaults.
   restationed = defaults;
   restationed.z_tair = defaults.z_tair + 1;
   restationed.z_wind = NaN;
   testCase.verifyEqual( ...
      icemodel.verification.helpers.physicsFingerprint( ...
      restationed).opts_sha256, baseline.opts_sha256)
end
