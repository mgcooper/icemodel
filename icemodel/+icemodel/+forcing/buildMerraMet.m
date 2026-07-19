function [met, metadata, Data] = buildMerraMet(location, years, kwargs)
   %BUILDMERRAMET Build an icemodel met timetable from MERRA-2 data.
   %
   %  [met, metadata] = icemodel.forcing.buildMerraMet(location, years)
   %  [met, metadata, Data] = ... buildMerraMet(location, years)
   %  [met, metadata] = ... buildMerraMet(_, source_dir=..., modis_dir=..., ...
   %     dt_out="15m")
   %
   % Extracts the MERRA-2 forcing at a point (or polygon average) and
   % converts it to the icemodel met contract: the Data-channel extraction
   % of icemodel.forcing.buildMerraData followed by icemodel.forcing.data2met
   % (MERRA already carries a total-precipitation channel, so ppt passes
   % through directly). Save the result with
   % icemodel.forcing.helpers.writemet.
   %
   % Inputs
   %  location - [lat lon] point, polyshape (EPSG:3413 m), or an Nx2 [lat lon]
   %             list of points. A point list returns a 1xN cell of met
   %             timetables (one per point); a single location returns one.
   %  years    - calendar years to extract
   %
   % Name-value
   %  source_dir, modis_dir, method, remap : see
   %      icemodel.forcing.buildMerraData (method = point sampling;
   %      remap = "conservative" (default, exact overlap-area weighting via
   %      exactremap); "equal" = plain in-polygon cell-centre mean
   %  fillgaps : opt in to legacy metchecks gap filling (default false)
   %  dt_out : output timestep resampled with the interval-support shared helper;
   %           default "15m", pass "" for native hourly
   %  fillwithmissing : add absent required channels as NaN placeholders
   %                    before met validation (default true)
   %
   % Outputs
   %  met      - met-contract timetable (hourly, or dt_out)
   %  metadata - finalized met metadata; exactly met.Properties.UserData
   %  Data     - source Data timetable before conversion/resampling
   %
   % MERRA provides a structurally valid met-contract timetable, while native
   % source gaps (including SNICEALB/albedo gaps) remain explicit by default.
   % Callers that require direct-run forcing must separately establish strict
   % forcing readiness or select a provenance-bearing repair policy. The
   % companion Data builder remains available for auxiliary channels.
   %
   % See also: icemodel.forcing.buildMerraData, icemodel.forcing.data2met,
   %  icemodel.forcing.buildMarMet, icemodel.forcing.helpers.writemet,
   %  icemodel.forcing.helpers.resampleMetTimestep

   arguments
      location
      years (1, :) double {mustBeInteger}
      kwargs.source_dir (1, 1) string = ""
      kwargs.modis_dir (1, 1) string = ""
      kwargs.fillgaps (1, 1) logical = false
      kwargs.method (1, 1) string {mustBeMember(kwargs.method, ...
         ["nearest", "natural"])} = "nearest"
      kwargs.remap (1, 1) string {mustBeMember(kwargs.remap, ...
         ["equal", "conservative"])} = "conservative"
      kwargs.dt_out (1, 1) string = "15m"
      kwargs.fillwithmissing (1, 1) logical = true
   end

   % buildMerraData accepts a single location (returns one Data timetable) or a
   % list of N points (Nx2 [lat lon], returns a 1xN cell of Data timetables).
   % Build the source Data first, then apply the shared Data-to-met collection
   % conversion. The returned shape mirrors Data: one timetable for one
   % location and a 1xN cell for a point list.
   [Data, ~] = icemodel.forcing.buildMerraData(location, years, ...
      source_dir=kwargs.source_dir, modis_dir=kwargs.modis_dir, ...
      fillgaps=kwargs.fillgaps, method=kwargs.method, remap=kwargs.remap);

   [met, metadata] = icemodel.forcing.helpers.data2metCollection(Data, ...
      dt_out=kwargs.dt_out, fillwithmissing=kwargs.fillwithmissing);
end
