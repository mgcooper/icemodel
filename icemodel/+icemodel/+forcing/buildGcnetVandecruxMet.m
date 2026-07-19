function [met, metadata, Data] = buildGcnetVandecruxMet(station, kwargs)
   %BUILDGCNETVANDECRUXMET Build native met from GC-Net/Vandecrux surface data.
   %
   %  [met, metadata] = icemodel.forcing.buildGcnetVandecruxMet(station)
   %  [met, metadata, Data] = ... buildGcnetVandecruxMet(station, source_dir=...)
   %
   % Builds a canonical Data timetable with
   % icemodel.forcing.buildGcnetVandecruxData, then converts it through the
   % shared data2met contract. Missing required channels become explicit NaN
   % placeholders when fillwithmissing=true, not fabricated zeros.
   %
   % See also: icemodel.forcing.buildGcnetVandecruxData,
   %  icemodel.forcing.data2met, icemodel.forcing.helpers.writemet

   arguments
      station (1, 1) string
      kwargs.source_dir (1, 1) string = ""
      kwargs.startdate = ""
      kwargs.enddate = ""
      kwargs.fillgaps (1, 1) logical = false
      kwargs.fillwithmissing (1, 1) logical = true
   end

   [Data, ~] = icemodel.forcing.buildGcnetVandecruxData(station, ...
      source_dir=kwargs.source_dir, startdate=kwargs.startdate, ...
      enddate=kwargs.enddate, fillgaps=kwargs.fillgaps);

   % data2met derives ppt from rainf+snowf. Because rainf is intentionally
   % all-NaN for this source, ppt stays all-NaN and validates as an explicit
   % missing-source placeholder when fillwithmissing is enabled. Use the same
   % collection-aware conversion path as every Data-backed met builder.
   [met, metadata] = icemodel.forcing.helpers.data2metCollection(Data, ...
      fillwithmissing=kwargs.fillwithmissing);
end
