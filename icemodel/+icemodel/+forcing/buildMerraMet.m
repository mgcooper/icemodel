function [met, metadata] = buildMerraMet(location, years, kwargs)
   %BUILDMERRAMET Build an icemodel met timetable from MERRA-2 data.
   %
   %  [met, metadata] = icemodel.forcing.buildMerraMet(location, years)
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
   %  location - [lat lon] point or polyshape in EPSG:3413 meters
   %  years    - calendar years to extract
   %
   % Name-value
   %  source_dir, modis_dir, fillgaps, method, remap : see
   %      icemodel.forcing.buildMerraData (method = point sampling;
   %      remap = "conservative" (default, exact overlap-area weighting via
   %      exactremap); "equal" = plain in-polygon cell-centre mean
   %  dt_out : optional output timestep ("15m") resampled per calendar year
   %           with icemodel.interpmet; default keeps hourly
   %
   % Outputs
   %  met      - met-contract timetable (hourly, or dt_out)
   %  metadata - provenance from buildMerraData
   %
   % Legacy: there was no MERRA met builder in runoff/ (MERRA was used for
   % Data files only); this mirrors icemodel.forcing.buildMarMet so MERRA can
   % drive a run directly.
   %
   % See also: icemodel.forcing.buildMerraData, icemodel.forcing.data2met,
   %  icemodel.forcing.buildMarMet, icemodel.forcing.helpers.writemet,
   %  icemodel.interpmet

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

   [Data, metadata] = icemodel.forcing.buildMerraData(location, years, ...
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
