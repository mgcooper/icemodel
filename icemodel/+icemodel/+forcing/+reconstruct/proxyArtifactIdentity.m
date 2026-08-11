function tf = proxyArtifactIdentity(metadata, site, location, product)
   %PROXYARTIFACTIDENTITY Verify one staged proxy's target and producer.
   %
   %  tf = icemodel.forcing.reconstruct.proxyArtifactIdentity( ...
   %     metadata, site, location, product)
   %
   % The filename is only a catalog hint. The staged timetable metadata must
   % identify the requested target point. It must also identify the exact
   % producer that the catalog storage token names.

   arguments
      metadata
      site (1, 1) string {icemodel.forcing.reconstruct.mustBeStationToken}
      location (1, 1) struct
      product (1, 1) string
   end

   tf = false;
   if ~isstruct(metadata) || ~isscalar(metadata) ...
         || ~all(isfield(metadata, ["lat_wgs84", "lon_wgs84"]))
      return
   end
   has_target = isequal(double(metadata.lat_wgs84), ...
      double(location.lat_wgs84)) ...
      && isequal(double(metadata.lon_wgs84), double(location.lon_wgs84));

   % Some canonical proxy artifacts carry no site token. When the token is
   % present it must agree. When it is absent, the exact target coordinates
   % give the identity.
   has_site = true;
   if isfield(metadata, 'site')
      has_site = (ischar(metadata.site) ...
         || (isstring(metadata.site) && isscalar(metadata.site))) ...
         && strcmpi(string(metadata.site), site);
   end

   has_source = false;
   switch product
      case "mar3.11"
         required = ["mar_qc_status", "source_files"];
         if all(isfield(metadata, required))
            status = string(metadata.mar_qc_status);
            files = string(metadata.source_files);
            has_source = isscalar(status) ...
               && ismember(status, ["applied", "not_applicable"]) ...
               && ~isempty(files) && all(contains(files, "MARv3.11-"));
         end
      case "merra2"
         has_source = ...
            icemodel.forcing.helpers.hasCanonicalMerraTimeSupport(metadata);
   end
   tf = has_target && has_site && has_source;
end
