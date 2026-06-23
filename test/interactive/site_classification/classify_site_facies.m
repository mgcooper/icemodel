function T = classify_site_facies(varargin)
   %CLASSIFY_SITE_FACIES Authoritative site classification from real datasets.
   %
   %  T = classify_site_facies()
   %  T = classify_site_facies("write_review", true)
   %
   %  Spatially samples three authoritative reference datasets at every PROMICE
   %  AWS site (and the ESM-SnowMIP sites for permafrost) to derive an
   %  AUTHORITATIVE surface_zone + permafrost_zone classification, replacing the
   %  elevation-band FIRST-PASS heuristic in
   %  icemodel.verification.helpers.promicesiteinfo.
   %
   %  THIS SCRIPT REQUIRES /Volumes/S03 (mounted) plus the staged SUMup_2025
   %  density file under data/verification/sumup. It is the analysis tool; its
   %  RESULTS are hard-coded into promicesiteinfo so the committed catalog has NO
   %  S03 runtime dependency. Re-run this when the datasets or sites change, then
   %  paste the emitted table into promicesiteinfo.
   %
   %  DATASETS
   %  --------
   %  1. MODIS end-of-summer BARE-ICE EXTENT (PRIMARY: ablation vs accumulation):
   %       /Volumes/S03/DATA/greenland/racmo2p3/annual/
   %         MODIS_Bare_Ice_Extent_2000-2018.nc
   %     Variable BareIceExtent is a BINARY (0/1) annual end-of-summer bare-ice
   %     mask on the native RACMO2.3p2 FGRN055 ~1 km grid (1496x2700), 19 annual
   %     slices 2000-2018, with 2-D LON/LAT. A "1" means the cell was exposed
   %     bare ice at the end of the melt season (the snowline had retreated past
   %     it). The mean late-summer snowline ~ the equilibrium line, so this is the
   %     proper end-of-summer ELA proxy. We compute the bare-ice FREQUENCY
   %     f_bare = N(year bare) / N(valid years) per cell and take the MAX over a
   %     3x3 native-cell window (margin sites can have their nearest single cell
   %     land off-ice and read 0; the window-max recovers the ablation signal one
   %     cell inland without inflating interior sites, which are ~0 window-wide).
   %     f_bare >= 0.50 (bare ice in a MAJORITY of years) =>
   %     ablation; below => the site is in the accumulation area at the surface.
   %     This REPRODUCES the KAN anchors: KAN_L f_bare=1.00, KAN_M=0.95 =>
   %     ablation; KAN_U=0.00 => accumulation-at-surface (snow-covered all years;
   %     KAN_U is the lower PERCOLATION zone by firn-core truth, pinned below).
   %
   %  2. Obu et al. (2019) permafrost ZONES (permafrost_zone, OFF-ICE sites):
   %       /Volumes/S03/DATA/interface/permafrost/obu/UiO_PEX_PERZONES_5/wgs/
   %         UiO_PEX_PERZONES_5.0_20181128_2000_2016_NH.shp
   %     ESA GlobPermafrost / UiO permafrost-extent zones, WGS84 polygons,
   %     attribute EXTENT in {Cont, Discon, Spora, Isol} ->
   %     continuous/discontinuous/sporadic/isolated (matching
   %     activelayer.readobuzones' parsing). We point-in-polygon test each OFF-ICE
   %     site (tundra/bedrock/not-Greenland/local-glacier) and each ESM-SnowMIP
   %     site directly on lon/lat (the WGS variant needs no projection). On-ice
   %     ice-sheet sites sit on ice, not permafrost ground -> "none".
   %
   %     The Obu shapefile is read through activelayer.readobuzones (the
   %     production reader), invoked with variant="wgs" to load the WGS84 lon/lat
   %     geometry and asshapefile/asgeostruct=true to return the raw geostruct
   %     for the point-in-polygon query. activelayer + its matfunclib helper
   %     dependencies are placed on the path by
   %     icemodel.test.helpers.bootstrapTestEnvironment (the readers earlier
   %     "failed" only because matfunclib was off the path). EXTENT->zone mapping
   %     (Cont/Discon/Spora/Isol) matches readobuzones' parsing exactly.
   %
   %  3. SUMup_2025 GrIS density profiles (FIRN facies within accumulation):
   %       data/verification/sumup/SUMup_2025_density_greenland.nc
   %     /DATA/latitude,/DATA/longitude carry ~2M firn/snow density measurement
   %     points. Within the accumulation area (not bare ice), a SUMup density
   %     profile co-located within ~15 km is direct EVIDENCE of firn presence ->
   %     percolation/firn facies. Otherwise we leave the honest coarse value
   %     "accumulation" (the surface data cannot resolve percolation vs dry snow).
   %
   %  surface_zone LOGIC
   %  ------------------
   %    off-ice (tundra)                       -> tundra, Obu permafrost_zone.
   %    off-ice (bedrock/not Greenland/snowmip) -> land,   Obu permafrost_zone.
   %    local glacier                          -> ablation (marginal).
   %    ice sheet -> by MODIS f_bare:
   %        f_bare >= 0.50                      -> ablation     (frequently bare).
   %        f_bare <  0.50 & SUMup within 15 km -> percolation  (firn observed).
   %        f_bare <  0.50 otherwise            -> accumulation (honest coarse).
   %    (A former elev >= 2500 m & f_bare==0 -> dry_snow branch was removed: the
   %     hard elevation cutoff did not generalize. dry_snow is now never emitted;
   %     accumulation-area facies is percolation via SUMup density else accumulation.)
   %      MODIS no-data at the cell -> fall back to elevation band, flagged.
   %
   %  permafrost_zone is ORTHOGONAL to surface_zone: ice-sheet/glacier sites get
   %  "none"; off-ice sites get the Obu extent class.
   %
   %  KAN ANCHORS (KAN_L/M=ablation, KAN_U=percolation) are USER-AUTHORITATIVE and
   %  pin over the method as the sanity check. KAN_U reading f_bare=0 / accumulation
   %  at the surface while being percolation by firn-core truth is EXPECTED (surface
   %  remote sensing cannot see subsurface percolation).
   %
   % See also: icemodel.verification.helpers.promicesiteinfo,
   %  icemodel.verification.namelists.surfacezone,
   %  icemodel.verification.namelists.permafrostzone,
   %  activelayer.readobuzones

   arguments (Repeating)
      varargin
   end
   p = inputParser;
   p.addParameter("write_review", true, @(x) islogical(x) || isnumeric(x));
   p.addParameter("sumup_radius_km", 15, @isnumeric);
   p.addParameter("ablation_freq", 0.50, @isnumeric);
   p.parse(varargin{:});
   opts = p.Results;

   S03 = "/Volumes/S03";
   assert(isfolder(S03), "S03 is not mounted at %s", S03);

   here = fileparts(mfilename("fullpath"));
   repo = string(icemodel.internal.fullpath("data"));
   csv = fullfile(repo, "verification", "promice", "AWS_sites_metadata.csv");

   % ---- Assemble the site list (PROMICE + ESM-SnowMIP for permafrost) --------
   sites = readSiteList(csv);

   % ---- Load the dataset handles --------------------------------------------
   modis = openModis(S03);
   obu = openObu(S03);
   sumup = openSumup(repo);

   % ---- Classify every site --------------------------------------------------
   n = numel(sites);
   T = table('Size', [n 13], ...
      'VariableTypes', ["string","double","double","double","string", ...
         "double","double","logical","double","string","string","string","string"], ...
      'VariableNames', ["site","lat","lon","elev","location_type", ...
         "modis_fbare","modis_nyr","sumup_firn","sumup_km", ...
         "obu_extent","surface_zone","eval_target","permafrost_zone"]);
   note = strings(n, 1);

   for k = 1:n
      s = sites(k);
      T.site(k) = s.site;
      T.lat(k) = s.lat;
      T.lon(k) = s.lon;
      T.elev(k) = s.elev;
      T.location_type(k) = s.location_type;

      if ~isfinite(s.lat) || ~isfinite(s.lon)
         T.surface_zone(k) = "unknown";
         T.eval_target(k) = "";
         T.permafrost_zone(k) = "unknown";
         T.modis_fbare(k) = NaN; T.modis_nyr(k) = 0;
         T.sumup_firn(k) = false; T.sumup_km(k) = NaN;
         T.obu_extent(k) = "n/a";
         note(k) = "no coordinates in metadata";
         continue
      end

      % --- MODIS bare-ice frequency (end-of-summer) ---
      [fbare, nyr] = sampleModis(modis, s.lat, s.lon);
      T.modis_fbare(k) = fbare;
      T.modis_nyr(k) = nyr;

      % --- SUMup firn evidence (nearest density profile) ---
      skm = sampleSumup(sumup, s.lat, s.lon);
      T.sumup_km(k) = skm;
      T.sumup_firn(k) = skm <= opts.sumup_radius_km;

      % --- Obu permafrost extent (off-ice only; ice-sheet -> none) ---
      ext = sampleObu(obu, s.lat, s.lon);
      T.obu_extent(k) = ext;

      % --- Assign surface_zone + eval_target + permafrost_zone ---
      [zone, target, pfz, nt] = assign(s, fbare, nyr, skm, ext, opts);

      % KAN anchors are USER-AUTHORITATIVE (firn-core ground truth). MODIS
      % surface bare-ice extent cannot see subsurface percolation, so KAN_U
      % (snow-covered surface, f_bare~0) is pinned to "percolation".
      anchor = kanAnchor(s.site);
      if anchor.is_anchor
         if zone ~= anchor.zone
            nt = sprintf("%s | PINNED to user anchor %s (method said %s)", ...
               nt, anchor.zone, zone);
         end
         zone = anchor.zone;
         target = anchor.target;
      end

      T.surface_zone(k) = zone;
      T.eval_target(k) = strjoin(target, ";");
      T.permafrost_zone(k) = pfz;
      note(k) = nt;
   end
   T.note = note;

   % ---- Sanity-check the KAN anchors -----------------------------------------
   checkKanAnchors(T);

   % ---- Emit review table (generated tables go to artifacts/, not figures/) --
   disp(T)
   if opts.write_review
      artdir = fullfile(here, "artifacts");
      icemodel.helpers.ensureDirExists(artdir);
      writeReview(T, fullfile(artdir, "site_facies_classification.md"), opts);
      writetable(T, fullfile(artdir, "site_facies_classification.csv"));
      fprintf("Review written to %s\n", ...
         fullfile(artdir, "site_facies_classification.md"));
   end
end

%% ------------------------------------------------------------------ site list
function sites = readSiteList(csv)
   Tp = readtable(csv, "TextType", "string");
   sites = struct("site", {}, "lat", {}, "lon", {}, "elev", {}, ...
      "location_type", {}, "kind", {});
   for n = 1:height(Tp)
      site = string(Tp.site_id(n));
      lat = todouble(Tp.latitude_installation(n));
      lon = todouble(Tp.longitude_installation(n));
      elev = todouble(Tp.altitude_installation(n));

      % Some stations (e.g. the ZAC GlacioBasis transect) carry no
      % latitude_installation/longitude_installation in the CSV. Source the
      % site lon/lat from the per-station L3 NetCDF (global latitude/longitude
      % attributes that readPromiceAws also reads) before falling back to the
      % CSV last-valid coordinates, so these sites still classify.
      if ~isfinite(lat) || ~isfinite(lon)
         [lat_nc, lon_nc, elev_nc] = ncCoords(site);
         if isfinite(lat_nc) && isfinite(lon_nc)
            lat = lat_nc; lon = lon_nc;
            if ~isfinite(elev), elev = elev_nc; end
         else
            lat = todouble(Tp.latitude_last_valid(n));
            lon = todouble(Tp.longitude_last_valid(n));
            if ~isfinite(elev)
               elev = todouble(Tp.altitude_last_valid(n));
            end
         end
      end

      sites(end+1) = struct( ...
         "site", site, ...
         "lat", lat, ...
         "lon", lon, ...
         "elev", elev, ...
         "location_type", lower(strtrim(string(Tp.location_type(n)))), ...
         "kind", "promice"); %#ok<AGROW>
   end
   % ESM-SnowMIP sites: permafrost descriptor only (off-ice land surfaces).
   smip = icemodel.verification.helpers.snowmipinfo();
   coords = struct( ... % lat lon from Menard et al. 2019 site table
      "cdp",[45.30 5.77], "oas",[53.629 -106.198], "obs",[53.987 -105.118], ...
      "ojp",[53.916 -104.692], "rme",[43.064 -116.755], "sap",[43.08 141.34], ...
      "snb",[37.907 -107.726], "sod",[67.362 26.633], "swa",[37.907 -107.711], ...
      "wfj",[46.83 9.81]);
   for n = 1:numel(smip)
      sc = char(smip(n).sitename);
      ll = coords.(sc);
      sites(end+1) = struct("site", upper(string(sc)), "lat", ll(1), ...
         "lon", ll(2), "elev", NaN, "location_type", "snowmip_land", ...
         "kind", "snowmip"); %#ok<AGROW>
   end
end

function [lat, lon, elev] = ncCoords(site)
   %NCCOORDS Read a station's lon/lat/elev from its per-station L3 NetCDF.
   %
   % The PROMICE/GC-Net L3 hourly files carry global latitude/longitude (and,
   % when present, altitude) attributes - the same coordinates readPromiceAws
   % reads. This recovers coordinates for stations whose CSV installation
   % coordinates are blank. Returns NaNs when the file or attributes are absent.
   lat = NaN; lon = NaN; elev = NaN;
   ncfile = string(fullfile(icemodel.internal.fullpath("data"), ...
      "verification", "promice", "hour", site + "_hour.nc"));
   if ~isfile(ncfile)
      return
   end
   info = ncinfo(char(ncfile));
   attrs = string({info.Attributes.Name});
   if ismember("latitude", attrs)
      lat = scalarAttr(ncreadatt(char(ncfile), "/", "latitude"));
   end
   if ismember("longitude", attrs)
      lon = scalarAttr(ncreadatt(char(ncfile), "/", "longitude"));
   end
   for ealt = ["altitude", "elevation"]
      if ismember(ealt, attrs)
         elev = scalarAttr(ncreadatt(char(ncfile), "/", char(ealt)));
         break
      end
   end
end

function v = scalarAttr(x)
   %SCALARATTR Reduce a NetCDF attribute to one finite scalar (or NaN).
   %
   % L3 files store the global latitude/longitude as CHAR strings, so parse
   % char attributes with str2double rather than casting their ASCII codes.
   if ischar(x) || isstring(x)
      v = str2double(string(x));
      return
   end
   x = double(x(:));
   x = x(isfinite(x));
   if isempty(x)
      v = NaN;
   else
      v = x(1);
   end
end

function v = todouble(x)
   if isstring(x) || ischar(x)
      v = str2double(x);
   else
      v = double(x);
   end
   if isempty(v), v = NaN; end
end

%% --------------------------------------------------------------------- MODIS
function m = openModis(S03)
   f = fullfile(S03, "DATA", "greenland", "racmo2p3", "annual", ...
      "MODIS_Bare_Ice_Extent_2000-2018.nc");
   assert(isfile(f), "MODIS bare-ice extent file not found: %s", f);
   mv = -1.000000015047466e+30;
   m.lat = ncread(f, "LAT");   % 1496x2700
   m.lon = ncread(f, "LON");
   m.bie = ncread(f, "BareIceExtent");  % 1496x2700x19, 0/1
   m.lat(m.lat == mv) = NaN;
   m.lon(m.lon == mv) = NaN;
   m.bie(m.bie == mv) = NaN;
   m.sz = size(m.lat);
end

function [fbare, nyr] = sampleModis(m, lat, lon)
   %SAMPLEMODIS Multi-year bare-ice frequency over a 3x3 native-cell window.
   %
   % A site exactly at the ice margin can have its single NEAREST cell land on an
   % off-ice (always-0) cell, hiding a clear ablation signal one cell inland
   % (verified at KAN_T/RED_L: nearest cell freq=0 but window reaches 1.0). To be
   % robust to this margin sampling artifact we take the MAXIMUM annual bare-ice
   % frequency over a 3x3 window of native ~1 km cells. Interior accumulation
   % sites have ~0 bare ice across the whole window, so the window-max does not
   % inflate them; it only recovers the ablation signal for margin sites.
   d = (m.lat - lat).^2 + (m.lon - lon).^2;
   [~, ix] = min(d(:));
   [ir, ic] = ind2sub(m.sz, ix);
   rr = max(1, ir-1):min(m.sz(1), ir+1);
   cc = max(1, ic-1):min(m.sz(2), ic+1);
   w = m.bie(rr, cc, :);
   nyr_cell = sum(~isnan(w), 3);
   freq_cell = sum(w == 1, 3) ./ max(nyr_cell, 1);
   freq_cell(nyr_cell == 0) = NaN;
   if all(isnan(freq_cell(:)))
      fbare = NaN; nyr = 0;
   else
      fbare = max(freq_cell(:), [], "omitnan");
      nyr = max(nyr_cell(:));
   end
end

%% ----------------------------------------------------------------- Obu zones
function b = openObu(S03)
   %OPENOBU Read the Obu (UiO PEX) permafrost zones via activelayer.readobuzones.
   %
   % Reads the WGS84 lon/lat variant of the UiO PEX permafrost-zones shapefile
   % through the activelayer toolbox reader (the production Obu reader), instead
   % of duplicating shaperead. readobuzones(variant="wgs", asshapefile=true,
   % asgeostruct=true) returns the raw geostruct (X/Y/BoundingBox/EXTENT) that
   % the point-in-polygon query in sampleObu needs. activelayer + matfunclib are
   % placed on the path by icemodel.test.helpers.bootstrapTestEnvironment.
   pn = fullfile(S03, "DATA", "interface", "permafrost", "obu", ...
      "UiO_PEX_PERZONES_5");
   fn = "UiO_PEX_PERZONES_5.0_20181128_2000_2016_NH.shp";
   assert(isfile(fullfile(pn, "wgs", fn)), ...
      "Obu zones shapefile not found: %s", fullfile(pn, "wgs", fn));
   assert(~isempty(which("activelayer.readobuzones")), ...
      ["activelayer.readobuzones not on path; run " ...
       "icemodel.test.helpers.bootstrapTestEnvironment to add " ...
       "projects/activelayer/toolbox + projects/matfunclib."]);
   b.S = activelayer.readobuzones(pathname=char(pn), filename=char(fn), ...
      variant="wgs", asshapefile=true, asgeostruct=true);
   b.ext = strings(numel(b.S), 1);
   for j = 1:numel(b.S)
      b.ext(j) = string(b.S(j).EXTENT);
   end
end

function ext = sampleObu(b, lat, lon)
   % WGS84 polygons: point-in-polygon directly on lon/lat.
   for j = 1:numel(b.S)
      bb = b.S(j).BoundingBox;
      if lon >= bb(1,1) && lon <= bb(2,1) && lat >= bb(1,2) && lat <= bb(2,2)
         if inpolygon(lon, lat, b.S(j).X, b.S(j).Y)
            ext = obuCode(b.ext(j)); return
         end
      end
   end
   ext = "none-mapped";  % not inside any Obu permafrost polygon
end

function name = obuCode(c)
   %OBUCODE Map Obu EXTENT to the permafrostzone vocabulary (per readobuzones).
   switch upper(strtrim(c))
      case "CONT",   name = "continuous";
      case "DISCON", name = "discontinuous";
      case "SPORA",  name = "sporadic";
      case "ISOL",   name = "isolated";
      otherwise,     name = "unknown";
   end
end

%% --------------------------------------------------------------------- SUMup
function s = openSumup(repo)
   f = fullfile(repo, "verification", "sumup", "SUMup_2025_density_greenland.nc");
   assert(isfile(f), "SUMup density file not found: %s", f);
   s.lat = ncread(f, "/DATA/latitude");
   s.lon = ncread(f, "/DATA/longitude");
end

function km = sampleSumup(s, lat, lon)
   % great-circle-ish nearest distance (small-angle, fine at GrIS scale)
   dk = 111 .* hypot(s.lat - lat, (s.lon - lon) .* cosd(lat));
   km = min(dk);
end

%% ------------------------------------------------------------------- assign
function [zone, target, pfz, note] = assign(s, fbare, nyr, skm, ext, opts)
   seasonal_ice  = ["seasonal_snow"; "bare_ice"];
   seasonal_firn = ["seasonal_snow"; "firn"];
   target = strings(0,1);
   pfz = "none";

   % An off-ice site outside any Obu permafrost polygon ("none-mapped") sits on
   % permafrost-free ground -> the schema value is "none".
   pfz_office = ext;
   if pfz_office == "none-mapped"
      pfz_office = "none";
   end

   switch s.location_type
      case "tundra"
         zone = "tundra";
         pfz = pfz_office;
         note = sprintf("off-ice tundra; Obu EXTENT=%s", ext);
         return
      case "snowmip_land"
         zone = "land";
         pfz = pfz_office;
         note = sprintf("ESM-SnowMIP off-ice land; Obu EXTENT=%s", ext);
         return
      case {"bedrock", "not greenland"}
         zone = "land";
         pfz = pfz_office;
         note = sprintf("off-ice land (%s); Obu EXTENT=%s", s.location_type, ext);
         return
      case "local glacier"
         zone = "ablation";
         target = seasonal_ice;
         note = "marginal local glacier -> ablation (eval seasonal_snow+bare_ice)";
         return
   end

   % ---- ice sheet: MODIS-driven, refined by SUMup + elevation ----
   if nyr == 0 || isnan(fbare)
      if ~isfinite(s.elev)
         zone = "unknown"; target = strings(0,1);
         note = "ice sheet, MODIS no-data, no elevation -> unknown";
      elseif s.elev < 1500
         zone = "ablation"; target = seasonal_ice;
         note = "ice sheet, MODIS no-data -> elev-band ablation (FALLBACK)";
      elseif s.elev < 2000
         zone = "percolation"; target = seasonal_firn;
         note = "ice sheet, MODIS no-data -> elev-band percolation (FALLBACK)";
      else
         zone = "accumulation"; target = seasonal_firn;
         note = "ice sheet, MODIS no-data -> elev-band accumulation (FALLBACK)";
      end
      return
   end

   if fbare >= opts.ablation_freq
      zone = "ablation"; target = seasonal_ice;
      note = sprintf("ice sheet, MODIS bare-ice freq=%.2f -> ablation", fbare);
      return
   end

   % accumulation area at the surface: refine the firn facies
   if skm <= opts.sumup_radius_km
      zone = "percolation"; target = seasonal_firn;
      note = sprintf("ice sheet, freq=%.2f, SUMup density %.1f km <= %.0f km -> percolation (firn observed)", ...
         fbare, skm, opts.sumup_radius_km);
   else
      zone = "accumulation"; target = seasonal_firn;
      note = sprintf("ice sheet, freq=%.2f, no SUMup firn within %.0f km -> accumulation (facies unresolved)", ...
         fbare, opts.sumup_radius_km);
   end
end

%% --------------------------------------------------------------- KAN anchor
function a = kanAnchor(site)
   %KANANCHOR User-authoritative KAN anchor zone/target, or not-an-anchor.
   a = struct("is_anchor", false, "zone", "", "target", strings(0,1));
   switch site
      case "KAN_L"
         a = struct("is_anchor", true, "zone", "ablation", ...
            "target", ["seasonal_snow"; "bare_ice"]);
      case "KAN_M"
         a = struct("is_anchor", true, "zone", "ablation", ...
            "target", ["seasonal_snow"; "bare_ice"]);
      case "KAN_U"
         a = struct("is_anchor", true, "zone", "percolation", ...
            "target", ["seasonal_snow"; "firn"]);
   end
end

%% --------------------------------------------------------------- KAN sanity
function checkKanAnchors(T)
   fprintf("\n--- KAN anchor sanity check (method vs user-authoritative) ---\n");
   anchors = struct("KAN_L","ablation","KAN_M","ablation","KAN_U","percolation");
   fn = fieldnames(anchors);
   for k = 1:numel(fn)
      row = T(T.site == string(fn{k}), :);
      if isempty(row), continue; end
      expected = anchors.(fn{k});
      got = row.surface_zone;
      flag = "OK"; if got ~= expected, flag = "MISMATCH"; end
      fprintf("  %-6s expected=%-12s final=%-12s bare-ice-freq=%.2f  [%s]\n", ...
         fn{k}, expected, got, row.modis_fbare, flag);
   end
   fprintf("\n");
end

%% ------------------------------------------------------------- review writer
function writeReview(T, path, opts)
   fid = fopen(path, "w");
   c = onCleanup(@() fclose(fid));
   fprintf(fid, "# Site facies classification (authoritative)\n\n");
   fprintf(fid, "Derived by `test/interactive/classify_site_facies.m` from:\n");
   fprintf(fid, "- MODIS end-of-summer **Bare Ice Extent** 2000-2018 (ablation vs accumulation; ablation when bare-ice freq >= %.2f),\n", opts.ablation_freq);
   fprintf(fid, "- **SUMup_2025** GrIS density profiles (firn evidence within %.0f km),\n", opts.sumup_radius_km);
   fprintf(fid, "- **Obu et al. (2019)** permafrost zones (off-ice permafrost_zone).\n\n");
   fprintf(fid, "KAN_L/M/U are USER-AUTHORITATIVE anchors and pin over the method.\n\n");
   fprintf(fid, "| site | lat | lon | elev | loc_type | bare_ice_freq | n_yr | sumup_firn | sumup_km | obu_extent | surface_zone | eval_target | permafrost | note |\n");
   fprintf(fid, "|------|-----|-----|------|----------|---------------|------|-----------|----------|------------|--------------|-------------|------------|------|\n");
   for k = 1:height(T)
      fprintf(fid, "| %s | %.3f | %.3f | %.0f | %s | %.2f | %d | %d | %.1f | %s | %s | %s | %s | %s |\n", ...
         T.site(k), T.lat(k), T.lon(k), T.elev(k), T.location_type(k), ...
         T.modis_fbare(k), double(T.modis_nyr(k)), double(T.sumup_firn(k)), ...
         T.sumup_km(k), T.obu_extent(k), T.surface_zone(k), T.eval_target(k), ...
         T.permafrost_zone(k), T.note(k));
   end
end
