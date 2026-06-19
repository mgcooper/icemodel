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
   %  3. With fillwinter=true (default), November-February samples are
   %     set to winter_albedo (default 0.8, a dry-snow value), because
   %     the first/last valid values bounding the polar night are often
   %     low late-summer values that would otherwise back-fill the
   %     winter.
   %
   % Inputs
   %  albedo - albedo series with winter gaps [-]
   %  Time   - datetimes of the samples
   %
   % Outputs
   %  albedo - gap-free albedo series [-]
   %
   % Legacy: reimplements runoff/functions/fillPromiceAlbedo.m (retained,
   % unchanged, as the legacy reference). Operates on a vector + Time instead
   % of the legacy matrix + header; fixes the legacy `nyears` undefined-
   % variable bug; drops the legacy plot option.
   %
   % See also: icemodel.forcing.buildPromiceMet,
   %  icemodel.forcing.helpers.metchecks

   arguments
      albedo (:, 1) double
      Time (:, 1) datetime
      kwargs.fillwinter (1, 1) logical = true
      kwargs.winter_albedo (1, 1) double = 0.8
   end

   assert(numel(albedo) == numel(Time), ...
      'albedo and Time must have equal lengths')

   albedo(albedo < 0 | albedo > 1) = NaN;

   for yyyy = unique(year(Time))'
      inyear = year(Time) == yyyy;
      filled = fillmissing(albedo(inyear), 'linear', 'EndValues', 'nearest');
      if kwargs.fillwinter
         m = month(Time(inyear));
         filled(m <= 2 | m >= 11) = kwargs.winter_albedo;
      end
      albedo(inyear) = filled;
   end
end
