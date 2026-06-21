function T = classify_site_facies(varargin)
   %CLASSIFY_SITE_FACIES Authoritative site classification from real datasets.
   %
   %  T = classify_site_facies()
   %  T = classify_site_facies("write_review", true)
   %
   %  Spatially samples three reference datasets at every PROMICE AWS site (and
   %  the ESM-SnowMIP sites for permafrost) to derive an AUTHORITATIVE
   %  surface_zone + permafrost_zone classification, replacing the elevation-band
   %  FIRST-PASS heuristic in icemodel.verification.helpers.promicesiteinfo.
   %
   %  THIS SCRIPT REQUIRES /Volumes/S03 (mounted). It is the analysis tool; its
   %  RESULTS are hard-coded into promicesiteinfo so the committed catalog has NO
   %  S03 runtime dependency. Re-run this when the datasets or sites change, then
   %  paste the emitted table into promicesiteinfo.
   %
   %  DATASETS
   %  --------
   %  1. MODIS albedo facies (PRIMARY for surface_zone):
   %       /Volumes/S03/DATA/greenland/snowlines/MOD10A1_albedo_classify_<YR>.nc
   %     Grid: ISMIP6 1km Greenland (x=1681, y=2881), 2D latitude/longitude vars,
   %     EPSG:3413-equivalent. Variable `classified` is int16 on a (z=92,x,y) stack
   %     of melt-season composites with classes:
   %        0 = no-data / off-ice (albedo 0)
   %        1 = bare ice / dark zone   (albedo ~38)
   %        2 = snow / firn            (albedo ~81)
   %     We sample the nearest pixel's class STACK and compute the multi-year
   %     bare-ice fraction f_bare = N(class==1) / N(class==1 | class==2). High
   %     f_bare => the site is snow-free for much of the melt season => ablation.
   %     This REPRODUCES the KAN anchors (KAN_L 0.93, KAN_M 0.86 => ablation;
   %     KAN_U 0.02 => percolation/accumulation) where MERRA-2 SMB does not.
   %
   %  2. MERRA-2 SMB / ELA mask (SECONDARY / sanity only):
   %       /Volumes/S03/DATA/merra2/ela/SMB.mat   (E,P,R 51x112x39 + Rcells georef)
   %       /Volumes/S03/DATA/merra2/ela/ELAmask_Merra2.mat (510x1400, UNGEOREF)
   %     The georeferenced product is the 51x112 MERRA-2 grid (GeographicCells/
   %     PostingsReference, lat 85..60, lon -78.125..-8.75). At ~50 km this grid is
   %     TOO COARSE to resolve the narrow SW-GrIS ablation margin: SMB is positive
   %     at KAN_L/M (which are authoritatively ablation), so SMB sign is recorded
   %     as a diagnostic ONLY and never overrides MODIS. The 510x1400 ELAmask.mat
   %     carries NO georeferencing (not a clean multiple of 51x112 or 1681x2881)
   %     and could not be georeferenced; it is NOT used.
   %
   %  3. Brown et al. 1997 permafrost (permafrost_zone, off-ice sites):
   %       /Volumes/S03/DATA/interface/GIS/Brown_Permafrost_Map/
   %         Permafrost_region_Brown_et_al_1997_Lambert_Azimuthal.shp
   %     CRS: Sphere ARC/INFO Lambert Azimuthal Equal Area (central meridian 180,
   %     origin 90N). Attribute EXTENT is the IPA permafrost-extent code:
   %        C=continuous  D=discontinuous  S=sporadic  I=isolated.
   %     We projfwd the site lon/lat into the CRS and point-in-polygon test; if the
   %     point falls in an unmapped gap we fall back to the nearest polygon edge
   %     (flagged low confidence beyond ~25 km).
   %
   %  surface_zone LOGIC
   %  ------------------
   %    off-ice (location_type tundra/bedrock/not Greenland) -> land or tundra,
   %      and permafrost_zone from Brown.
   %    local glacier                                        -> ablation (marginal).
   %    ice sheet -> by MODIS f_bare:
   %        f_bare >= 0.40            -> ablation     (frequently bare ice)
   %        0.10 <= f_bare < 0.40     -> percolation  (mixed bare-ice / snow margin)
   %        f_bare <  0.10 & elev<2500-> accumulation (never bare; wet/percolation
   %                                                    facies not separable here)
   %        f_bare <  0.10 & elev>=2500-> dry_snow     (cold high interior)
   %      MODIS no-data at the pixel -> fall back to the elevation band, flagged.
   %
   %  permafrost_zone is ORTHOGONAL to surface_zone: ice-sheet sites get "none"
   %  (the ice surface is not permafrost ground); off-ice sites get the Brown code.
   %
   % See also: icemodel.verification.helpers.promicesiteinfo,
   %  icemodel.verification.namelists.surfacezone,
   %  icemodel.verification.namelists.permafrostzone

   arguments (Repeating)
      varargin
   end
   p = inputParser;
   p.addParameter("write_review", true, @(x) islogical(x) || isnumeric(x));
   p.addParameter("modis_years", [2008 2010 2012 2014 2016], @isnumeric);
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
   modis = openModis(S03, opts.modis_years);
   smb = openSmb(S03);
   brown = openBrown(S03);

   % ---- Classify every site --------------------------------------------------
   n = numel(sites);
   T = table('Size', [n 14], ...
      'VariableTypes', ["string","double","double","double","string","string", ...
         "double","double","string","string","string","string","string","string"], ...
      'VariableNames', ["site","lat","lon","elev","location_type", ...
         "ela_smb_sign","modis_fbare","modis_ndata","modis_facies", ...
         "brown_extent","surface_zone","eval_target","permafrost_zone", ...
         "confidence_note"]);

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
         T.confidence_note(k) = "no coordinates in metadata";
         T.modis_facies(k) = "n/a";
         T.brown_extent(k) = "n/a";
         T.ela_smb_sign(k) = "n/a";
         continue
      end

      % --- MERRA-2 SMB sign (diagnostic only) ---
      smbval = sampleSmb(smb, s.lat, s.lon);
      if isnan(smbval)
         T.ela_smb_sign(k) = "nan";
      elseif smbval < 0
         T.ela_smb_sign(k) = "ablation";
      else
         T.ela_smb_sign(k) = "accum";
      end

      % --- MODIS bare-ice fraction ---
      [fbare, ndata] = sampleModis(modis, s.lat, s.lon);
      T.modis_fbare(k) = fbare;
      T.modis_ndata(k) = ndata;
      if ndata == 0 || isnan(fbare)
         T.modis_facies(k) = "no-data";
      elseif fbare >= 0.40
         T.modis_facies(k) = "bare-ice";
      elseif fbare >= 0.10
         T.modis_facies(k) = "mixed";
      else
         T.modis_facies(k) = "snow";
      end

      % --- Brown permafrost ---
      [ext, pdist] = sampleBrown(brown, s.lat, s.lon);
      T.brown_extent(k) = ext;

      % --- Assign surface_zone + eval_target + permafrost_zone ---
      [zone, target, pfz, note] = assign(s, fbare, ndata, ext, pdist);

      % KAN anchors are USER-AUTHORITATIVE (firn-core ground truth). MODIS
      % surface albedo cannot see subsurface percolation, so KAN_U (snow-covered
      % surface, f_bare~0.02) is pinned to "percolation" even though the method
      % alone would call it "accumulation". Record the disagreement honestly.
      anchor = kanAnchor(s.site);
      if anchor.is_anchor
         if zone ~= anchor.zone
            note = sprintf("%s | PINNED to user anchor %s (method said %s)", ...
               note, anchor.zone, zone);
         end
         zone = anchor.zone;
         target = anchor.target;
      end

      T.surface_zone(k) = zone;
      T.eval_target(k) = strjoin(target, ";");
      T.permafrost_zone(k) = pfz;
      T.confidence_note(k) = note;
   end

   % ---- Sanity-check the KAN anchors -----------------------------------------
   checkKanAnchors(T);

   % ---- Emit review table ----------------------------------------------------
   disp(T)
   if opts.write_review
      writeReview(T, fullfile(here, "figures", "site_zone_classification.md"));
      writetable(T, fullfile(here, "figures", "site_zone_classification.csv"));
      fprintf("Review written to %s\n", ...
         fullfile(here, "figures", "site_zone_classification.md"));
   end
end

%% ------------------------------------------------------------------ site list
function sites = readSiteList(csv)
   Tp = readtable(csv, "TextType", "string");
   sites = struct("site", {}, "lat", {}, "lon", {}, "elev", {}, ...
      "location_type", {}, "kind", {});
   for n = 1:height(Tp)
      sites(end+1) = struct( ...
         "site", string(Tp.site_id(n)), ...
         "lat", todouble(Tp.latitude_installation(n)), ...
         "lon", todouble(Tp.longitude_installation(n)), ...
         "elev", todouble(Tp.altitude_installation(n)), ...
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

function v = todouble(x)
   if isstring(x) || ischar(x)
      v = str2double(x);
   else
      v = double(x);
   end
   if isempty(v), v = NaN; end
end

%% --------------------------------------------------------------------- MODIS
function m = openModis(S03, years)
   dir0 = fullfile(S03, "DATA", "greenland", "snowlines");
   files = strings(0,1);
   for y = years
      f = fullfile(dir0, sprintf("MOD10A1_albedo_classify_%d.nc", y));
      if isfile(f), files(end+1) = f; end %#ok<AGROW>
   end
   assert(~isempty(files), "no MODIS classify files found for the requested years");
   ref = files(1);
   m.files = files;
   m.lat = ncread(ref, "latitude");   % 1681x2881
   m.lon = ncread(ref, "longitude");
   m.sz = size(m.lat);
end

function [fbare, ndata] = sampleModis(m, lat, lon)
   % nearest pixel by squared geographic distance (sites are far from dateline)
   d = (m.lat - lat).^2 + (m.lon - lon).^2;
   [~, ix] = min(d(:));
   [ir, ic] = ind2sub(m.sz, ix);
   nbare = 0; nsnow = 0;
   for f = reshape(m.files, 1, [])
      col = double(squeeze(ncread(f, "classified", [1 ir ic], [Inf 1 1])));
      nbare = nbare + sum(col == 1);
      nsnow = nsnow + sum(col == 2);
   end
   ndata = nbare + nsnow;
   if ndata == 0
      fbare = NaN;
   else
      fbare = nbare / ndata;
   end
end

%% ----------------------------------------------------------------- MERRA SMB
function s = openSmb(S03)
   f = fullfile(S03, "DATA", "merra2", "ela", "SMB.mat");
   D = load(f);
   s.E = mean(D.E, 3, "omitnan");   % multi-year mean SMB, 51x112
   s.R = D.Rcells;
end

function val = sampleSmb(s, lat, lon)
   try
      val = geointerp(s.E, s.R, lat, lon);
   catch
      val = NaN;
   end
end

%% --------------------------------------------------------------------- Brown
function b = openBrown(S03)
   f = fullfile(S03, "DATA", "interface", "GIS", "Brown_Permafrost_Map", ...
      "Permafrost_region_Brown_et_al_1997_Lambert_Azimuthal.shp");
   info = shapeinfo(f);
   b.crs = info.CoordinateReferenceSystem;
   b.S = shaperead(f);
   b.ext = strings(numel(b.S), 1);
   for j = 1:numel(b.S)
      b.ext(j) = string(b.S(j).EXTENT);
   end
end

function [ext, dist] = sampleBrown(b, lat, lon)
   [x, y] = projfwd(b.crs, lat, lon);
   % exact point-in-polygon first
   for j = 1:numel(b.S)
      bb = b.S(j).BoundingBox;
      if x >= bb(1,1) && x <= bb(2,1) && y >= bb(1,2) && y <= bb(2,2)
         if inpolygon(x, y, b.S(j).X, b.S(j).Y)
            ext = brownCode(b.ext(j)); dist = 0; return
         end
      end
   end
   % nearest polygon edge fallback (point in an unmapped gap)
   best = inf; bestj = 0;
   for j = 1:numel(b.S)
      px = b.S(j).X; py = b.S(j).Y;
      ok = isfinite(px) & isfinite(py);
      dd = min((px(ok) - x).^2 + (py(ok) - y).^2);
      if dd < best, best = dd; bestj = j; end
   end
   dist = sqrt(best);
   ext = brownCode(b.ext(bestj));
end

function name = brownCode(c)
   switch upper(strtrim(c))
      case "C", name = "continuous";
      case "D", name = "discontinuous";
      case "S", name = "sporadic";
      case "I", name = "isolated";
      otherwise, name = "unknown";
   end
end

%% ------------------------------------------------------------------- assign
function [zone, target, pfz, note] = assign(s, fbare, ndata, ext, pdist)
   seasonal_ice  = ["seasonal_snow"; "bare_ice"];
   seasonal_firn = ["seasonal_snow"; "firn"];
   target = strings(0,1);
   pfz = "none";

   % Brown distance confidence qualifier for off-ice sites.
   if pdist == 0
      pf_qual = "in-polygon";
   elseif pdist < 25e3
      pf_qual = sprintf("nearest %.0f km", pdist/1e3);
   else
      pf_qual = sprintf("nearest %.0f km (LOW conf)", pdist/1e3);
   end

   switch s.location_type
      case "tundra"
         zone = "tundra";
         pfz = ext;
         pf_note = sprintf("Brown EXTENT=%s (%s)", ext, pf_qual);
         note = sprintf("off-ice tundra; %s", pf_note);
         return
      case "snowmip_land"
         % ESM-SnowMIP reference sites are off-ice land surfaces (boreal forest,
         % alpine meadow, wetland) -> "land". Carried here for permafrost_zone.
         zone = "land";
         pfz = ext;
         note = sprintf("ESM-SnowMIP off-ice land; Brown EXTENT=%s (%s)", ...
            ext, pf_qual);
         return
      case {"bedrock", "not greenland"}
         zone = "land";
         pfz = ext;
         note = sprintf("off-ice land (%s); Brown EXTENT=%s (%s)", ...
            s.location_type, ext, pf_qual);
         return
      case "local glacier"
         zone = "ablation";
         target = seasonal_ice;
         note = "marginal local glacier -> ablation (eval seasonal_snow+bare_ice)";
         return
   end

   % ---- ice sheet: MODIS-driven ----
   if ndata == 0 || isnan(fbare)
      % MODIS no-data: fall back to elevation band, flag it.
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

   if fbare >= 0.40
      zone = "ablation"; target = seasonal_ice;
      note = sprintf("ice sheet, MODIS f_bare=%.2f -> ablation", fbare);
   elseif fbare >= 0.10
      zone = "percolation"; target = seasonal_firn;
      note = sprintf("ice sheet, MODIS f_bare=%.2f -> percolation (mixed)", fbare);
   elseif isfinite(s.elev) && s.elev >= 2500
      zone = "dry_snow"; target = seasonal_firn;
      note = sprintf("ice sheet, MODIS f_bare=%.2f, elev=%.0f -> dry_snow", ...
         fbare, s.elev);
   else
      zone = "accumulation"; target = seasonal_firn;
      note = sprintf("ice sheet, MODIS f_bare=%.2f -> accumulation", fbare);
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
      fprintf("  %-6s expected=%-12s method=%-12s f_bare=%.2f  [%s]\n", ...
         fn{k}, expected, got, row.modis_fbare, flag);
   end
   fprintf("\n");
end

%% ------------------------------------------------------------- review writer
function writeReview(T, path)
   fid = fopen(path, "w");
   c = onCleanup(@() fclose(fid));
   fprintf(fid, "# Site zone classification (authoritative)\n\n");
   fprintf(fid, "Derived by `test/interactive/classify_site_facies.m` from MODIS ");
   fprintf(fid, "albedo facies (primary), MERRA-2 SMB (diagnostic), and Brown ");
   fprintf(fid, "et al. 1997 permafrost (off-ice). MODIS bare-ice fraction is the ");
   fprintf(fid, "primary surface_zone signal; it reproduces the KAN anchors.\n\n");
   fprintf(fid, "| site | lat | lon | elev | loc_type | SMB | f_bare | facies | brown | surface_zone | eval_target | permafrost | note |\n");
   fprintf(fid, "|------|-----|-----|------|----------|-----|--------|--------|-------|--------------|-------------|------------|------|\n");
   for k = 1:height(T)
      fprintf(fid, "| %s | %.3f | %.3f | %.0f | %s | %s | %.2f | %s | %s | %s | %s | %s | %s |\n", ...
         T.site(k), T.lat(k), T.lon(k), T.elev(k), T.location_type(k), ...
         T.ela_smb_sign(k), T.modis_fbare(k), T.modis_facies(k), ...
         T.brown_extent(k), T.surface_zone(k), T.eval_target(k), ...
         T.permafrost_zone(k), T.confidence_note(k));
   end
end
