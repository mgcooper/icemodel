function thermal_conductivity_ice(T, reference)
   %THERMAL_CONDUCTIVITY_ICE Plot ice conductivity reference curves.
   %
   %  icemodel.plot.thermal_conductivity_ice()
   %  icemodel.plot.thermal_conductivity_ice(T)
   %  icemodel.plot.thermal_conductivity_ice(T, reference)
   %
   % Description:
   %  Plots one or more ice thermal-conductivity reference curves from
   %  icemodel.kernels.thermal_conductivity_ice over the temperature vector
   %  T [K].
   %
   % Inputs:
   %  T         - Temperature [K] (default: 240:280)
   %  reference - One or more conductivity references supported by
   %              icemodel.kernels.thermal_conductivity_ice
   %
   % See also: icemodel.kernels.thermal_conductivity_ice

   arguments (Input)
      T {mustBeNumeric} = (240:280).'
      reference (1, :) string {mustBeMember(reference, ...
         ["andersson_2005", "yen_1981", "rabin_2000", "petrenko_1999", ...
         "slack_1980", "engineering_toolbox"])} = ["andersson_2005", ...
         "yen_1981", "rabin_2000", "petrenko_1999", "slack_1980", ...
         "engineering_toolbox"]
   end

   Tf = 273.15;

   figure("Position", [200 200 700 500])
   hold on
   box on

   for ref = reference(:)'
      plot(T-Tf, icemodel.kernels.thermal_conductivity_ice(T, ref));
   end

   axis tight
   xline(0, 'LineWidth', 1.5)
   legend(strrep(reference, "_", " "), "NumColumns", 2, "Location", "southwest")
   ylabel('k [W m^{-1} K^{-1}]')
   xlabel('Temperature [^oC]')
   title("Thermal Conductivity of Ice Ih")

   ax = gca;
   ax.XAxis.TickLength = ax.XAxis.TickLength / 2;
   ax.YAxis.TickLength = ax.YAxis.TickLength / 2;
   ax.XAxis.TickDirection = "in";
   ax.YAxis.TickDirection = "in";
end
