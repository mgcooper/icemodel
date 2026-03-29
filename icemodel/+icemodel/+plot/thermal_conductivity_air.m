function thermal_conductivity_air(Tc, reference)
   %THERMAL_CONDUCTIVITY_AIR Plot dry-air thermal conductivity diagnostics.
   %
   %  icemodel.plot.thermal_conductivity_air()
   %  icemodel.plot.thermal_conductivity_air(Tc)
   %  icemodel.plot.thermal_conductivity_air(Tc, reference)
   %
   % Description:
   %  Plots dry-air thermal conductivity, its temperature derivative, and
   %  its ratio to the repo reference constant over the Celsius vector Tc.
   %  One or more kernel formulations can be compared on the same axes. The
   %  snow-relevant range from -60 to 0 C is shaded.
   %
   % Inputs:
   %  Tc        - Temperature [°C] (default: -60:10)
   %  reference - One or more conductivity formulations supported by
   %              icemodel.kernels.thermal_conductivity_air
   %
   % See also: icemodel.kernels.thermal_conductivity_air

   arguments
      Tc {mustBeNumeric} = (-60:10).'
      reference (1, :) string {mustBeMember(reference, ...
         ["lemmon_jacobsen_2004", "sutherland"])} = [ ...
         "lemmon_jacobsen_2004", "sutherland"]
   end

   T = Tc + 273.15;  % K

   k_ref = icemodel.physicalConstant('k_air'); % W m^-1 K^-1
   xshade = [-60, 0];

   f = figure( ...
      'Name', 'Dry-air thermal conductivity diagnostics', ...
      'Color', 'w', ...
      'Position', [10 10 891 787]);

   tl = tiledlayout(f, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
   title(tl, 'Dry-air thermal conductivity diagnostics')

   ax1 = nexttile(tl, 1);
   hold(ax1, 'on')
   xregion(ax1, xshade(1), xshade(2), ...
      'FaceColor', [0.92 0.92 0.92], ...
      'EdgeColor', 'none')
   xline(ax1, 0, '-', 'LineWidth', 1.5)
   yline(ax1, k_ref, '-', 'k_{ref}', ...
      'LabelHorizontalAlignment', 'left', ...
      'LabelVerticalAlignment', 'top', ...
      'FontSize', 16, ...
      'LineWidth', 1.5)
   xlim(ax1, [-60, 10])
   box(ax1, 'on')
   grid(ax1, 'on')
   ylabel(ax1, 'k [W m^{-1} K^{-1}]')
   title(ax1, 'Thermal conductivity')

   ax1.XAxis.TickLength = ax1.XAxis.TickLength / 2;
   ax1.YAxis.TickLength = ax1.YAxis.TickLength / 2;
   ax1.XAxis.TickDirection = "in";
   ax1.YAxis.TickDirection = "in";

   ax2 = nexttile(tl, 2);
   hold(ax2, 'on')
   xregion(ax2, xshade(1), xshade(2), ...
      'FaceColor', [0.92 0.92 0.92], ...
      'EdgeColor', 'none')
   xline(ax2, 0, '-', 'LineWidth', 1.5)
   xlim(ax2, [-60, 10])
   box(ax2, 'on')
   grid(ax2, 'on')
   ylabel(ax2, 'dk/dT [W m^{-1} K^{-2}]')
   title(ax2, 'Temperature derivative')

   ax2.XAxis.TickLength = ax2.XAxis.TickLength / 2;
   ax2.YAxis.TickLength = ax2.YAxis.TickLength / 2;
   ax2.XAxis.TickDirection = "in";
   ax2.YAxis.TickDirection = "in";

   ax3 = nexttile(tl, 3);
   hold(ax3, 'on')
   xregion(ax3, xshade(1), xshade(2), ...
      'FaceColor', [0.92 0.92 0.92], ...
      'EdgeColor', 'none')
   xline(ax3, 0, '-', 'LineWidth', 1.5)
   yline(ax3, 1.0, '-', 'unity assumption', ...
      'LabelHorizontalAlignment', 'left', ...
      'LabelVerticalAlignment', 'top', ...
      'FontSize', 16, ...
      'LineWidth', 1.5)
   xlim(ax3, [-60, 10])
   box(ax3, 'on')
   grid(ax3, 'on')
   xlabel(ax3, 'Temperature [^oC]')
   ylabel(ax3, 'k / k_{ref} [-]')
   title(ax3, 'Ratio used in snow conductivity')

   ax3.XAxis.TickLength = ax3.XAxis.TickLength / 2;
   ax3.YAxis.TickLength = ax3.YAxis.TickLength / 2;
   ax3.XAxis.TickDirection = "in";
   ax3.YAxis.TickDirection = "in";

   solution_handles = gobjects(numel(reference), 1);
   for i = 1:numel(reference)
      ref = reference(i);
      [k, dkdT] = icemodel.kernels.thermal_conductivity_air(T, ref);
      kratio = k ./ k_ref;

      solution_handles(i) = plot(ax1, Tc, k, 'LineWidth', 1.5);
      plot(ax2, Tc, dkdT, 'LineWidth', 1.5)
      plot(ax3, Tc, kratio, 'LineWidth', 1.5)
   end

   legend(ax1, solution_handles, strrep(reference, "_", " "), ...
      'Location', 'southeast', 'Interpreter', 'none')
   linkaxes([ax1, ax2, ax3], 'x');
end
