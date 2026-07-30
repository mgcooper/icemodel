function metadata = geusModisCoverageMetadata(requested_years, coverage_years)
   %GEUSMODISCOVERAGEMETADATA Build canonical GEUS MODIS coverage provenance.
   %
   %  metadata = ...
   %     icemodel.forcing.helpers.geusModisCoverageMetadata( ...
   %     requested_years, coverage_years)
   %
   % REQUESTED_YEARS are the calendar years on the target artifact axis.
   % COVERAGE_YEARS are the requested years backed by one unambiguous GEUS
   % Greenland Reflectivity 5 km C6 source file and physical target values.
   % The returned flat struct is the single product/status/year contract shared
   % by fresh RCM builders, metadata repair, saved payloads, and artifact QA.

   arguments
      requested_years double {mustBeFinite, mustBeInteger}
      coverage_years double {mustBeFinite, mustBeInteger}
   end

   % Normalize orientation, ordering, and duplicates before classifying. Sorted
   % years make top-level/UserData equality independent of caller vector shape.
   requested_years = unique(requested_years(:))';
   coverage_years = unique(coverage_years(:))';
   outside = setdiff(coverage_years, requested_years);
   if ~isempty(outside)
      error('icemodel:forcing:geusModisCoverageMetadata:outsideRequest', ...
         'MODIS coverage year(s) are outside the requested axis: %s', ...
         strjoin(string(outside), ', '))
   end

   % Empty coverage is explicit absence; otherwise distinguish full coverage
   % from a source-backed subset without inventing values for missing years.
   if isempty(coverage_years)
      status = "no_source_coverage";
   elseif isequal(coverage_years, requested_years)
      status = "source_coverage";
   else
      status = "partial_source_coverage";
   end

   metadata = struct( ...
      'modis_product', 'GEUS Greenland Reflectivity 5km C6', ...
      'modis_status', char(status), ...
      'modis_coverage_years', coverage_years(:));
end
