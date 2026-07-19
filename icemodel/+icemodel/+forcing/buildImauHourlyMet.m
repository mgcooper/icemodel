function [met, metadata, Data] = buildImauHourlyMet(station, kwargs)
   %BUILDIMAUHOURLYMET Build native met from IMAU hourly AWS data.
   %
   %  [met, metadata] = icemodel.forcing.buildImauHourlyMet("S21")
   %  [met, metadata, Data] = ... buildImauHourlyMet("S21", source_dir=...)
   %
   % Builds canonical IMAU Data with buildImauHourlyData, then converts it
   % through the shared data2met contract. Missing precipitation channels stay
   % explicit NaN placeholders when fillwithmissing=true.
   %
   % See also: icemodel.forcing.buildImauHourlyData,
   %  icemodel.forcing.data2met

   arguments
      station (1, 1) string
      kwargs.source_dir (1, 1) string = ""
      kwargs.startdate = ""
      kwargs.enddate = ""
      kwargs.fillgaps (1, 1) logical = false
      kwargs.fillwithmissing (1, 1) logical = true
   end

   [Data, ~] = icemodel.forcing.buildImauHourlyData(station, ...
      source_dir=kwargs.source_dir, startdate=kwargs.startdate, ...
      enddate=kwargs.enddate, fillgaps=kwargs.fillgaps);

   % Use the shared collection-aware conversion path; this source returns one
   % native hourly timetable, so the default blank dt_out is an exact no-op.
   [met, metadata] = icemodel.forcing.helpers.data2metCollection(Data, ...
      fillwithmissing=kwargs.fillwithmissing);
end
