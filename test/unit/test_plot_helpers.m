function tests = test_plot_helpers
   %TEST_PLOT_HELPERS Cover the shared report-figure helpers.
   %
   % Report builders share these so figures look the same and stay out of the
   % legend. markTimeSpan gained a fill style for panels that highlight many
   % spans at once, where boundary lines would be unreadable.
   tests = functiontests(localfunctions);
end

function ax = scratchAxes(testCase)
   %SCRATCHAXES Make an invisible axes that closes when the test ends.
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));
   ax = axes(fig);
end

function test_line_style_draws_two_boundary_lines(testCase)
   % The default marks both ends and stays out of the legend.
   ax = scratchAxes(testCase);
   t0 = datetime(2019, 7, 1);
   returned = icemodel.plot.markTimeSpan(ax, t0, t0 + days(3));

   testCase.verifyNumElements(returned, 2)
   testCase.verifyEqual(string({returned.HandleVisibility}), ["off", "off"])
end

function test_fill_style_draws_one_shaded_region(testCase)
   % The fill style shades the interval instead, and is also hidden from the
   % legend so overlay labels stay clean.
   ax = scratchAxes(testCase);
   t0 = datetime(2019, 7, 1);
   returned = icemodel.plot.markTimeSpan(ax, t0, t0 + days(3), ...
      style="fill", color=[0.8 0.9 0.8], face_alpha=0.25);

   testCase.verifyNumElements(returned, 1)
   testCase.verifyEqual(returned.HandleVisibility, 'off')
   testCase.verifyEqual(returned.FaceAlpha, 0.25, AbsTol=1e-12)
end

function test_unknown_style_is_rejected(testCase)
   % The style set is closed, so a typo fails instead of silently drawing the
   % default.
   ax = scratchAxes(testCase);
   t0 = datetime(2019, 7, 1);
   testCase.verifyError(@() icemodel.plot.markTimeSpan(ax, t0, ...
      t0 + days(1), style="shade"), 'MATLAB:validators:mustBeMember')
end

function test_new_figure_uses_the_requested_size(testCase)
   % Report figures must export at a stable pixel size.
   fig = icemodel.plot.newFigure(width=800, height=300);
   testCase.addTeardown(@() close(fig));

   testCase.verifyEqual(fig.Position(3:4), [800 300])
   testCase.verifyEqual(fig.Color, [1 1 1])
end
