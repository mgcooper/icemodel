function albedo = fillPromiceAlbedo(albedo, Time, kwargs)
   %FILLPROMICEALBEDO Fill PROMICE albedo gaps with the winter-fill policy.
   %
   %  albedo = icemodel.forcing.fillPromiceAlbedo(albedo, Time)
   %  albedo = ... fillPromiceAlbedo(_, fillwinter=false, winter_albedo=0.8)
   %
   % PROMICE albedo (Albedo_theta<70d) is undefined through the polar
   % winter when the sun is below the 70-degree zenith cutoff, so the raw
   % series has months-long gaps every year. This applies the legacy
   % winter-fill policy, year by year:
   %
   %  1. Values outside [0, 1] are set missing.
   %  2. Leading/trailing gaps in each calendar year fill with the first/
   %     last valid value of that year; interior gaps fill linearly.
   %  3. With fillwinter=true (default), this function sets winter-month
   %     samples to winter_albedo, a dry-snow value. The constant and the
   %     month window both come from icemodel.parameterLookup, so the
   %     reconstruction can detect this exact stamp. Without this step the
   %     first and last valid values around the polar night, often low
   %     late-summer values, would back-fill the winter.
   %
   % Inputs
   %  albedo - albedo series with winter gaps [-]
   %  Time   - datetimes of the samples
   %
   % Outputs
   %  albedo - gap-free albedo series [-]
   %
   % Legacy: this function reimplements runoff/functions/fillPromiceAlbedo.m.
   % The original stays unchanged as the legacy reference. This version takes a
   % vector and Time instead of the legacy matrix and header, fixes the legacy
   % `nyears` undefined-variable bug, and has no plot option.
   %
   % See also: icemodel.forcing.buildPromiceMet,
   %  icemodel.forcing.helpers.metchecks

   arguments
      albedo (:, 1) double
      Time (:, 1) datetime
      kwargs.fillwinter (1, 1) logical = true
      kwargs.winter_albedo (1, 1) double = ...
         icemodel.parameterLookup('promice_winter_albedo')
      kwargs.winter_months (1, :) double = ...
         icemodel.parameterLookup('promice_winter_albedo_months')
   end

   assert(numel(albedo) == numel(Time), ...
      'albedo and Time must have equal lengths')

   albedo(~icemodel.forcing.promiceAlbedoSourceValid(albedo)) = NaN;

   for yyyy = unique(year(Time))'
      inyear = year(Time) == yyyy;
      filled = fillmissing(albedo(inyear), 'linear', 'EndValues', 'nearest');
      if kwargs.fillwinter
         m = month(Time(inyear));
         filled(ismember(m, kwargs.winter_months)) = kwargs.winter_albedo;
      end
      albedo(inyear) = filled;
   end
end
