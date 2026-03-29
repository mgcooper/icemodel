function thermal_conductivity_water(Tc, rho, reference)
   %THERMAL_CONDUCTIVITY_WATER Plot liquid-water conductivity diagnostics.
   %
   %  icemodel.plot.thermal_conductivity_water()
   %  icemodel.plot.thermal_conductivity_water(Tc)
   %  icemodel.plot.thermal_conductivity_water(Tc, rho)
   %  icemodel.plot.thermal_conductivity_water(Tc, rho, reference)
   %
   % Description:
   %  Plots liquid-water conductivity and its constant-density temperature
   %  derivative over the Celsius vector Tc using
   %  icemodel.kernels.thermal_conductivity_water. Multiple density curves
   %  can be compared on the same axes.
   %
   % Inputs:
   %  Tc        - Temperature [°C] (default: -20:0.5:20)
   %  rho       - One or more liquid-water densities [kg m-3]
   %  reference - Conductivity formulation supported by
   %              icemodel.kernels.thermal_conductivity_water
   %
   % See also: icemodel.kernels.thermal_conductivity_water

   arguments
      Tc {mustBeNumeric} = (-20:0.5:20).'
      rho {mustBeNumeric} = 999.8395
      reference (1, :) string {mustBeMember(reference, ...
         ["iapws_2011_linear_0c", "iapws_2011"])} = [ ...
         "iapws_2011_linear_0c", "iapws_2011"]
   end

   T = Tc + 273.15;

   f = figure( ...
      'Name', 'Liquid-water thermal conductivity diagnostics', ...
      'Color', 'w', ...
      'Position', [200 200 880 680]);

   tl = tiledlayout(f, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
   title(tl, 'Liquid-water thermal conductivity diagnostics')

   ax1 = nexttile(tl, 1);
   hold(ax1, 'on')
   box(ax1, 'on')
   xline(ax1, 0, 'LineWidth', 1.0)
   grid(ax1, 'on')
   ylabel(ax1, 'k [W m^{-1} K^{-1}]')
   title(ax1, 'Thermal conductivity')

   ax1.XAxis.TickLength = ax1.XAxis.TickLength / 2;
   ax1.YAxis.TickLength = ax1.YAxis.TickLength / 2;
   ax1.XAxis.TickDirection = "in";
   ax1.YAxis.TickDirection = "in";

   ax2 = nexttile(tl, 2);
   hold(ax2, 'on')
   box(ax2, 'on')
   xline(ax2, 0, 'LineWidth', 1.0)
   grid(ax2, 'on')
   xlabel(ax2, 'Temperature [^oC]')
   ylabel(ax2, 'dk/dT [W m^{-1} K^{-2}]')
   title(ax2, 'Constant-density temperature derivative')

   ax2.XAxis.TickLength = ax2.XAxis.TickLength / 2;
   ax2.YAxis.TickLength = ax2.YAxis.TickLength / 2;
   ax2.XAxis.TickDirection = "in";
   ax2.YAxis.TickDirection = "in";

   solution_handles = gobjects(numel(reference) * numel(rho), 1);
   labels = strings(numel(reference) * numel(rho), 1);
   idx = 0;
   for ref = reference(:)'
      if ref == "iapws_2011_linear_0c"
         rho_iter = rho(1);
      else
         rho_iter = rho(:).';
      end

      for rhoi = rho_iter
         idx = idx + 1;
         [k, dkdT] = icemodel.kernels.thermal_conductivity_water(T, rhoi, ref);
         solution_handles(idx) = plot(ax1, Tc(:), k(:), 'LineWidth', 1.5);
         plot(ax2, Tc(:), dkdT(:), 'LineWidth', 1.5)
         if ref == "iapws_2011_linear_0c"
            labels(idx) = strrep(ref, "_", " ");
         else
            labels(idx) = string(sprintf('%s, rho = %.1f kg m^{-3}', ...
               char(strrep(ref, "_", " ")), rhoi));
         end
      end
   end

   solution_handles = solution_handles(1:idx);
   labels = labels(1:idx);
   legend(ax1, solution_handles, labels, 'Location', 'best', 'Interpreter', 'none')
   linkaxes([ax1, ax2], 'x');
end
