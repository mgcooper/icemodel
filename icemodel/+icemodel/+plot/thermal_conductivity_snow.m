function thermal_conductivity_snow(rho_dry, T, reference)
   %THERMAL_CONDUCTIVITY_SNOW Plot archived snow conductivity options.
   %
   %  icemodel.plot.thermal_conductivity_snow()
   %  icemodel.plot.thermal_conductivity_snow(rho_dry, T, reference)
   %
   % Description:
   %  Plots one or more archived conductivity references from the
   %  icemodel.kernels.thermal_conductivity_snow helper against dry snow
   %  density. This is intended as a diagnostic/reference plot rather than
   %  the production conductivity path.
   %
   % Inputs:
   %  rho_dry - Dry snow density [kg m-3]
   %  T       - Temperature [K] used by the temperature-dependent options
   %  reference - One or more conductivity references to plot
   %
   % See also: icemodel.kernels.thermal_conductivity_snow,
   %  icemodel.column.firn_thermal_conductivity

   arguments
      rho_dry {mustBeNumeric} = (50:10:900).'
      T (1, 1) {mustBeNumeric} = 263.15
      reference (1, :) string {mustBeMember(reference, [ ...
         "sturm_1997", "anderson_1976", "aggarwal_2009", "yen_1981", ...
         "calonne_2017", "van_dusen_1929", "schwerdtfeger_1963", ...
         "calonne_2019_eq1", "calonne_2019_eq5"])} = [ ...
         "sturm_1997", "anderson_1976", "aggarwal_2009", "yen_1981", ...
         "calonne_2017", "van_dusen_1929", "schwerdtfeger_1963", ...
         "calonne_2019_eq1", "calonne_2019_eq5"]
   end

   ro_ice = icemodel.physicalConstant('ro_ice');
   f_ice = rho_dry ./ ro_ice;

   figure('Name', 'Snow thermal conductivity', ...
      'Color', 'w', 'Position', [200 200 860 520])
   hold on
   box on

   solution_handles = gobjects(numel(reference), 1);
   for i = 1:numel(reference)
      ksnow = icemodel.kernels.thermal_conductivity_snow( ...
         T, f_ice, ro_ice, rho_dry, reference(i));
      solution_handles(i) = plot(rho_dry(:), ksnow(:), 'LineWidth', 1.5);
   end

   grid on
   xlabel('Dry snow density [kg m^{-3}]')
   ylabel('k [W m^{-1} K^{-1}]')
   title(sprintf('Snow conductivity references at T = %.2f K', T))
   legend(solution_handles, strrep(reference, "_", " "), ...
      'Location', 'northwest', 'Interpreter', 'none')

   ax = gca;
   ax.XAxis.TickLength = ax.XAxis.TickLength / 2;
   ax.YAxis.TickLength = ax.YAxis.TickLength / 2;
   ax.XAxis.TickDirection = "in";
   ax.YAxis.TickDirection = "in";
end
