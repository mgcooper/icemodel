function [met, metadata, Data] = buildKtransectMet(station, kwargs)
   %BUILDKTRANSECTMET Build native met from K-transect 30-minute AWS data.
   %
   %  [met, metadata] = icemodel.forcing.buildKtransectMet("AWS9")
   %  [met, metadata, Data] = ... buildKtransectMet("AWS9", source_dir=...)
   %
   % Builds canonical K-transect Data with buildKtransectData, then converts it
   % through the shared data2met contract. Missing precipitation channels stay
   % explicit NaN placeholders when fillwithmissing=true.
   %
   % See also: icemodel.forcing.buildKtransectData,
   %  icemodel.forcing.data2met

   arguments
      station (1, 1) string
      kwargs.source_dir (1, 1) string = ""
      kwargs.startdate = ""
      kwargs.enddate = ""
      kwargs.fillgaps (1, 1) logical = false
      kwargs.fillwithmissing (1, 1) logical = true
   end

   [Data, ~] = icemodel.forcing.buildKtransectData(station, ...
      source_dir=kwargs.source_dir, startdate=kwargs.startdate, ...
      enddate=kwargs.enddate, fillgaps=kwargs.fillgaps);

   % Use the shared collection-aware conversion path; this source returns one
   % native 30-minute timetable, so the default blank dt_out is an exact no-op.
   [met, metadata] = icemodel.forcing.helpers.data2metCollection(Data, ...
      fillwithmissing=kwargs.fillwithmissing);
end
