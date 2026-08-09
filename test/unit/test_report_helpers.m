function tests = test_report_helpers
   %TEST_REPORT_HELPERS Cover the shared report markdown and metric helpers.
   %
   % Both report builders render saved metadata through these, so a value
   % containing markup, HTML, or backticks must not be able to change the
   % document structure.
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
