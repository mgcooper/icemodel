function thermal_conductivity_firn(Tc, rho_firn, reference, depth)
   %THERMAL_CONDUCTIVITY_FIRN Plot firn conductivity references.
   %
   %  icemodel.plot.thermal_conductivity_firn()
   %  icemodel.plot.thermal_conductivity_firn(Tc, rho_firn, reference, depth)
   %
   % Description:
   %  Primarily plots firn conductivity as a function of temperature at a
   %  fixed firn density, so temperature-dependent firn references can be
   %  compared directly. The Oster and Albert depth relation is available as
   %  a special case when "oster_2020_depth" is included in REFERENCE.
   %
   % Inputs:
   %  Tc        - Temperature [°C]
   %  rho_firn  - Firn density [kg m-3]
   %  reference - One or more firn conductivity references
   %  depth     - Depth vector [m] for the special-case depth relation
   %
   % See also: icemodel.kernels.thermal_conductivity_firn

   arguments
      Tc {mustBeNumeric} = (-40:1:0).'
      rho_firn (1, 1) {mustBeNumeric} = 600
      reference (1, :) string {mustBeMember(reference, [ ...
         "calonne_2019_eq5", "calonne_2019_eq1", "schwerdtfeger_1963", ...
         "yen_1981", "van_dusen_1929", "oster_2020_density", ...
         "oster_2020_depth"])} = [ ...
         "calonne_2019_eq5", "calonne_2019_eq1", "schwerdtfeger_1963", ...
         "yen_1981", "van_dusen_1929", "oster_2020_density"]
      depth {mustBeNumeric} = (0:0.25:30).'
   end

   depth_reference = reference == "oster_2020_depth";
   temperature_reference = reference(~depth_reference);
   T = Tc + 273.15;

   if any(depth_reference)
      f = figure( ...
         'Name', 'Firn thermal conductivity diagnostics', ...
         'Color', 'w', ...
         'Position', [200 200 980 420]);

      tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
      title(tl, 'Firn thermal conductivity diagnostics')

      ax1 = nexttile(tl, 1);
      hold(ax1, 'on')
      box(ax1, 'on')
      grid(ax1, 'on')
      xlabel(ax1, 'Temperature [^oC]')
      ylabel(ax1, 'k [W m^{-1} K^{-1}]')
      title(ax1, ...
         sprintf('Temperature comparison at \\rho = %.0f kg m^{-3}', rho_firn))

      ax2 = nexttile(tl, 2);
      hold(ax2, 'on')
      box(ax2, 'on')
      grid(ax2, 'on')
      xlabel(ax2, 'Depth [m]')
      ylabel(ax2, 'k [W m^{-1} K^{-1}]')
      title(ax2, 'Special-case depth relation')
   else
      f = figure( ...
         'Name', 'Firn thermal conductivity diagnostics', ...
         'Color', 'w', ...
         'Position', [200 200 760 460]);

      tl = tiledlayout(f, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
      title(tl, 'Firn thermal conductivity diagnostics')

      ax1 = nexttile(tl, 1);
      hold(ax1, 'on')
      box(ax1, 'on')
      grid(ax1, 'on')
      xlabel(ax1, 'Temperature [^oC]')
      ylabel(ax1, 'k [W m^{-1} K^{-1}]')
      title(ax1, ...
         sprintf('Temperature comparison at \\rho = %.0f kg m^{-3}', rho_firn))
      ax2 = [];
   end

   solution_handles = gobjects(numel(temperature_reference), 1);
   for i = 1:numel(temperature_reference)
      k = icemodel.kernels.thermal_conductivity_firn( ...
         T, rho_firn, temperature_reference(i));
      if isscalar(k)
         k = repmat(k, size(Tc));
      end
      solution_handles(i) = plot(ax1, Tc(:), k(:), 'LineWidth', 1.5);
   end

   if ~isempty(solution_handles)
      legend(ax1, solution_handles, strrep(temperature_reference, "_", " "), ...
         'Location', 'northeast', 'NumColumns', 2)
   end

   ax1.XAxis.TickDirection = "in";
   ax1.YAxis.TickDirection = "in";
   ax1.XAxis.TickLength = ax1.XAxis.TickLength / 2;
   ax1.YAxis.TickLength = ax1.YAxis.TickLength / 2;

   if ~isempty(ax2)
      k_depth = icemodel.kernels.thermal_conductivity_firn( ...
         T(1), rho_firn, "oster_2020_depth", depth);
      plot(ax2, depth(:), k_depth(:), 'LineWidth', 1.5)
      ax2.XAxis.TickDirection = "in";
      ax2.YAxis.TickDirection = "in";
      ax2.XAxis.TickLength = ax2.XAxis.TickLength / 2;
      ax2.YAxis.TickLength = ax2.YAxis.TickLength / 2;
   end
end
