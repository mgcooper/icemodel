function [met, metadata] = buildMarMet(location, years, kwargs)
   %BUILDMARMET Build an icemodel met timetable from MAR v3.11 data.
   %
   %  [met, metadata] = icemodel.forcing.buildMarMet(location, years)
   %  [met, metadata] = ... buildMarMet(_, source_dir=..., modis_dir=..., ...
   %     dt_out="15m")
   %
   % Extracts the MAR forcing at a point (or polygon average) and
   % converts it to the icemodel met contract: the Data-channel
   % extraction of icemodel.forcing.buildMarData followed by
   % icemodel.forcing.data2met (ppt = snow + rain). Save the result with
   % icemodel.forcing.helpers.writemet.
   %
   % Inputs
   %  location - [lat lon] point or polyshape in EPSG:3413 meters
   %  years    - calendar years to extract
   %
   % Name-value
   %  source_dir, modis_dir, fillgaps, method, remap : see
   %      icemodel.forcing.buildMarData (method = point sampling;
   %      remap = "conservative" (default, exact overlap-area weighting via
   %      exactremap); "equal" = plain in-polygon cell-centre mean
   %  dt_out : optional output timestep ("15m") resampled per calendar
   %           year with icemodel.interpmet; default keeps hourly
   %
   % Outputs
   %  met      - met-contract timetable (hourly, or dt_out)
   %  metadata - provenance from buildMarData
   %
   % Legacy: reimplements runoff/functions/makeMarMetfile.m (the original
   % retained, unchanged, as the legacy reference workflow). The legacy
   % per-point loop, per-variable rounding, and in-function saving are not
   % reproduced (saving is icemodel.forcing.helpers.writemet's job).
   %
   % See also: icemodel.forcing.buildMarData, icemodel.forcing.data2met,
   %  icemodel.forcing.helpers.writemet, icemodel.interpmet

   arguments
      location
      years (1, :) double {mustBeInteger}
      kwargs.source_dir (1, 1) string = ""
      kwargs.modis_dir (1, 1) string = ""
      kwargs.fillgaps (1, 1) logical = true
      kwargs.method (1, 1) string {mustBeMember(kwargs.method, ...
         ["nearest", "natural"])} = "nearest"
      kwargs.remap (1, 1) string {mustBeMember(kwargs.remap, ...
         ["equal", "conservative"])} = "conservative"
      kwargs.dt_out (1, 1) string = ""
   end

   [Data, metadata] = icemodel.forcing.buildMarData(location, years, ...
      source_dir=kwargs.source_dir, modis_dir=kwargs.modis_dir, ...
      fillgaps=kwargs.fillgaps, method=kwargs.method, remap=kwargs.remap);

   met = icemodel.forcing.data2met(Data);

   if kwargs.dt_out ~= ""
      % icemodel.interpmet resamples one calendar year at a time.
      parts = cell(numel(years), 1);
      for n = 1:numel(years)
         parts{n} = icemodel.interpmet( ...
            met(year(met.Time) == years(n), :), char(kwargs.dt_out));
      end
      met = vertcat(parts{:});
   end
end
