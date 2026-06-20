function [forcing, metadata] = buildSumupForcing(point, years, kwargs)
   %BUILDSUMUPFORCING Build co-located MAR/RACMO forcing at a SUMup point.
   %
   %  [forcing, metadata] = ...
   %     icemodel.verification.setup.buildSumupForcing([lat lon], years)
   %  [forcing, metadata] = ...
   %     icemodel.verification.setup.buildSumupForcing([lat lon], years, ...
   %        mar_dir=..., racmo_dir=...)
   %
   %  Builds the co-located forcing/Data bundle staged alongside a SUMup firn
   %  observation point. SUMup points are not weather stations, so - unlike
   %  the PROMICE anchors - there is no station met: the forcing is the MAR
   %  point met (whose albedo is swapped downstream) and the reference Data is
   %  the co-located RACMO SMB/eval Data. This mirrors the MAR/RACMO legs of
   %  importPromiceSites' co-located bundle, reusing the same forcing builders
   %  so the conversion path is identical.
   %
   %  Inputs
   %    point : [lat lon]  WGS84 SUMup point coordinates.
   %    years : numeric vector  calendar years to extract.
   %
   %  Name-value
   %    mar_dir    : string ("")  MAR source dir (passed to buildMarMet).
   %    racmo_dir  : string ("")  RACMO source dir (passed to buildRacmoData).
   %    modis_dir  : string ("")  optional MODIS albedo dir.
   %    method     : string ("nearest")  point-sampling method.
   %    dt_out     : string ("")  optional output timestep for the met builder.
   %    window_start / window_end : datetime or "" optional clamp window.
   %
   %  Returns
   %    forcing  : struct  with fields:
   %                 mar_met    timetable  MAR point met
   %                 racmo_data timetable  RACMO point Data (SMB/eval)
   %    metadata : struct  provenance + sample method.
   %
   %  Role
   %    Reusable per-point SUMup forcing builder, symmetric with
   %    buildSumupObservations. Used by importSumup during staging.
   %
   % See also: icemodel.verification.setup.buildSumupObservations,
   %  icemodel.verification.setup.importSumup,
   %  icemodel.forcing.buildMarMet, icemodel.forcing.buildRacmoData

   arguments
      point (1, 2) double
      years (1, :) double
      kwargs.mar_dir (1, 1) string = ""
      kwargs.racmo_dir (1, 1) string = ""
      kwargs.modis_dir (1, 1) string = ""
      kwargs.method (1, 1) string = "nearest"
      kwargs.dt_out (1, 1) string = ""
      kwargs.window_start = ""
      kwargs.window_end = ""
   end

   % MAR point met (albedo swapped downstream like the PROMICE-anchor MAR leg).
   mar_met = icemodel.forcing.buildMarMet(point, years, ...
      source_dir=kwargs.mar_dir, modis_dir=kwargs.modis_dir, ...
      method=kwargs.method, dt_out=kwargs.dt_out);
   mar_met = windowSubset(mar_met, kwargs.window_start, kwargs.window_end);

   % RACMO point Data - the SMB/eval reference (never a met source).
   [racmo_data, ~] = icemodel.forcing.buildRacmoData(point, years, ...
      source_dir=kwargs.racmo_dir, modis_dir=kwargs.modis_dir, ...
      method=kwargs.method, dt="1hr");
   racmo_data = windowSubset(racmo_data, kwargs.window_start, kwargs.window_end);

   forcing = struct('mar_met', mar_met, 'racmo_data', racmo_data);

   metadata = icemodel.verification.setup.metadataStruct({ ...
      'forcing_source', 'MAR point met (albedo swapped downstream)'
      'reference_source', 'RACMO2.3p3 point Data (SMB/eval only)'
      'sample_method', char(kwargs.method)
      'point_lat_wgs84', point(1)
      'point_lon_wgs84', point(2)});
end

%% Local helpers
function tt = windowSubset(tt, t1, t2)
   %WINDOWSUBSET Clamp a timetable to [t1, t2] on a UTC-aware axis (no-op blank).
   if strcmp(string(t1), "") || strcmp(string(t2), "")
      return
   end
   t1 = icemodel.verification.setup.ensureUtc(t1);
   t2 = icemodel.verification.setup.ensureUtc(t2);
   t = tt.Time;
   if isempty(t.TimeZone)
      t.TimeZone = 'UTC';
   end
   keep = t >= t1 & t <= t2;
   tt = tt(keep, :);
end
