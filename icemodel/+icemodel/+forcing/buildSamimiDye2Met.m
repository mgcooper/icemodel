function [met, metadata, Data] = buildSamimiDye2Met(kwargs)
   %BUILDSAMIMIDYE2MET Build RetMIP Dye-2 2016 native met from Samimi AWS.
   %
   %  [met, metadata] = icemodel.forcing.buildSamimiDye2Met()
   %  [met, metadata, Data] = ... buildSamimiDye2Met(source_dir=...)
   %
   % Builds canonical Samimi Dye-2 Data with buildSamimiDye2Data, then converts
   % it through the shared data2met contract. Missing precipitation channels stay
   % explicit NaN placeholders when fillwithmissing=true.
   %
   % See also: icemodel.forcing.buildSamimiDye2Data,
   %  icemodel.forcing.data2met

   arguments
      kwargs.source_dir (1, 1) string = ""
      kwargs.startdate = ""
      kwargs.enddate = ""
      kwargs.fillgaps (1, 1) logical = false
      kwargs.fillwithmissing (1, 1) logical = true
   end

   [Data, ~] = icemodel.forcing.buildSamimiDye2Data( ...
      source_dir=kwargs.source_dir, startdate=kwargs.startdate, ...
      enddate=kwargs.enddate, fillgaps=kwargs.fillgaps);

   % Use the shared collection-aware conversion path; this source returns one
   % native half-hourly timetable, so the default blank dt_out is an exact no-op.
   [met, metadata] = icemodel.forcing.helpers.data2metCollection(Data, ...
      fillwithmissing=kwargs.fillwithmissing);
end
