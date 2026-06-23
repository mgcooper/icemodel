function [T, mdfile] = firn_coverage_analysis(kwargs)
   %FIRN_COVERAGE_ANALYSIS Profile firn-zone evaluation-data coverage per site.
   %
   %  T = firn_coverage_analysis()
   %  T = firn_coverage_analysis(write_markdown=true)
   %  [T, mdfile] = firn_coverage_analysis(sites=["KAN_U","EGP"])
   %
   %  Builds a per-site COVERAGE TABLE for the firn/accumulation candidate set
   %  (PROMICE/GC-Net Accumulation-type sites + KAN_U) profiling each across the
   %  dimensions a firn-model development effort cares about, then classifies
   %  each site into one or more firn-physics verification regimes and writes a
   %  markdown report. This is a USER-FACING analysis tool; its outputs (table +
   %  markdown + figures) go to this theme's gitignored
   %  test/interactive/firn_coverage/figures/ subdir.
   %
   %  THIS SCRIPT REQUIRES /Volumes/S03 (mounted) for the GrIS sector mask and
   %  the MAR/MERRA/RACMO coverage probe, plus the staged PROMICE L3 NetCDF and
   %  SUMup_2025 density+temperature greenland files under data/verification/.
   %
   %  DIMENSIONS PROFILED (per candidate site)
   %  ----------------------------------------
   %   GEOGRAPHY      GrIS basin/sector from lat/lon (Greenland_sectors.nc, the
   %                  Mouginot/Zwally 8-basin mask on the RACMO FGRN055 grid).
   %   ELEVATION      installation altitude from the L3 NetCDF metadata.
   %   SURFACE_ZONE   glaciological facies (promicesiteinfo, the authoritative
   %                  MODIS-bare-ice + SUMup-density classification).
   %   CLIMATE        accumulation magnitude proxy (SUMup SMB mean within radius,
   %                  else the PROMICE net surface-height trend) binned
   %                  low/med/high; MELT presence from the surface-height
   %                  seasonal drawdown amplitude (summer melt signal).
   %   MISSINGNESS    PROMICE record length (years), gap-flagged fraction of the
   %                  surface channel (surface_height_flag==1), and channel
   %                  completeness (tice10m present? thermistor string present?
   %                  co-located SUMup density / temperature present?).
   %   REGIME         firn-physics verification regime(s) the site can exercise:
   %                    (a) INFILTRATION  low/med-accum WITH melt  (percolation)
   %                    (b) AQUIFER       high-accum WITH melt (perennial water)
   %                    (c) DRY_FIRN      high-accum LOW melt (densification/thermal)
   %
   %  REGIME LOGIC (honest, surface-observable)
   %  -----------------------------------------
   %   The GC-Net accumulation transect is INTERIOR/high-elevation; it spans the
   %   cold-percolation and dry-firn regimes. The FIRN AQUIFER regime (Forster
   %   2014; Miege 2016) is a SE-Greenland high-accumulation + high-melt
   %   phenomenon (FA-13/FA-15 cores, Helheim/SE catchment) that the interior
   %   GC-Net set does NOT sample. This tool FLAGS that gap rather than forcing a
   %   site into the aquifer class. Aquifer coverage requires the SUMup SE cores
   %   (or a dedicated FA core) - reported as a GAP in the markdown.
   %
   % See also: icemodel.verification.helpers.promicesiteinfo,
   %  icemodel.verification.setup.importPromiceSites,
   %  icemodel.verification.setup.buildSumupObservations,
   %  test/interactive/classify_site_facies.m

   arguments
      kwargs.sites (1, :) string = defaultCandidates()
      kwargs.write_markdown (1, 1) logical = true
      kwargs.sumup_radius_km (1, 1) double = 15
   end

   S03 = "/Volumes/S03";
   assert(isfolder(S03), "S03 is not mounted at %s", S03);

   sites = kwargs.sites;
   here = fileparts(mfilename("fullpath"));
   figdir = fullfile(here, "figures");
   icemodel.helpers.ensureDirExists(figdir);

   % ---- Reference dataset handles -------------------------------------------
   sectors = openSectors(fullfile(S03, ...
      "DATA/greenland/racmo2p3/annual/Greenland_sectors.nc"));

   % On-disk RCM coverage (years), probed once, for the per-leg overlap column.
   models = ["mar", "merra", "racmo"];
   cov = icemodel.verification.setup.promiceSourceCoverage(models, ...
      struct('mar', "", 'merra', "", 'racmo', ""));

   % ---- Profile each candidate ----------------------------------------------
   rows = cell(numel(sites), 1);
   for n = 1:numel(sites)
      site = sites(n);
      fprintf("[%2d/%2d] profiling %s ...\n", n, numel(sites), site);
      rows{n} = profileSite(site, sectors, cov, kwargs.sumup_radius_km);
   end
   T = vertcat(rows{:});

   % Console summary.
   disp(T)

   mdfile = "";
   if kwargs.write_markdown
      mdfile = fullfile(figdir, "firn_coverage_analysis.md");
      writeMarkdown(mdfile, T, cov, kwargs.sumup_radius_km);
      fprintf("\nwrote coverage report: %s\n", mdfile);
   end
end

%% ------------------------------------------------------------------ profile
function row = profileSite(site, sectors, cov, sumup_radius_km)
   %PROFILESITE Build one coverage-table row for a candidate site.

   % --- Site metadata + L3 record window from the staged NetCDF. ---
   [~, meta] = icemodel.forcing.readPromiceAws(site, timescale="hourly");
   lat = meta.lat; lon = meta.lon; elev = meta.elev;
   rec_start = meta.window_start;
   rec_end = meta.window_end;
   rec_years = max(0, years(rec_end - rec_start));

   info = icemodel.verification.helpers.promicesiteinfo(site);

   % --- Geography: GrIS sector from the basin mask. ---
   sector = sampleSector(sectors, lat, lon);

   % --- PROMICE eval Data: surface channel, gap fraction, channels. ---
   [surf_gap_frac, surf_trend, tice10m_present, ...
      tice_string_present] = profilePromiceData(site, rec_start, rec_end);

   % --- SUMup co-location (density + subsurface temperature). ---
   [has_sumup_rho, has_sumup_temp, sumup_smb_mean] = ...
      profileSumup([lat lon], sumup_radius_km);

   % --- Climate: accumulation magnitude + melt presence. ---
   accum_class = classifyAccum(sumup_smb_mean, surf_trend);
   melt_present = classifyMelt(info.surface_zone, elev);

   % --- RCM overlap (years) per leg over the record window. ---
   ry = year(rec_start):year(rec_end);
   mar_ov = overlapYears(ry, cov.mar);
   merra_ov = overlapYears(ry, cov.merra);
   racmo_ov = overlapYears(ry, cov.racmo);

   % --- Firn-physics regime classification. ---
   regime = classifyRegime(info.surface_zone, melt_present);

   row = table( ...
      string(site), string(sector), round(elev), string(info.surface_zone), ...
      strjoin(string(info.eval_target), "+"), ...
      round(rec_years, 1), round(100 * surf_gap_frac), ...
      tice10m_present, tice_string_present, has_sumup_rho, has_sumup_temp, ...
      string(accum_class), melt_present, ...
      string(mar_ov), string(merra_ov), string(racmo_ov), ...
      string(regime), ...
      'VariableNames', {'site', 'sector', 'elev_m', 'surface_zone', ...
      'eval_target', 'rec_yrs', 'gap_pct', 'tice10m', 'tice_str', ...
      'sumup_rho', 'sumup_T', 'accum', 'melt', ...
      'mar_yrs', 'merra_yrs', 'racmo_yrs', 'regime'});
end

%% ----------------------------------------------------------- PROMICE Data
function [gap_frac, trend, tice10m_present, ...
      tice_string_present] = profilePromiceData(site, rec_start, rec_end)
   %PROFILEPROMICEDATA Surface-channel gap fraction, net trend, channel presence.
   gap_frac = NaN; trend = NaN;
   tice10m_present = false; tice_string_present = false;
   try
      [data, ~] = icemodel.forcing.buildPromiceData(site, ...
         startdate=rec_start, enddate=rec_end, frequency="daily");
      vn = string(data.Properties.VariableNames);

      % Surface channel: accumulation sites carry surface_height [m, +up].
      surfvar = "";
      if ismember("surface_height", vn)
         surfvar = "surface_height";
      elseif ismember("ablation", vn)
         surfvar = "ablation";
      end
      if surfvar ~= ""
         s = data.(surfvar);
         good = isfinite(s);
         if ismember("surface_height_flag", vn)
            f = data.surface_height_flag;
            gap_frac = mean(f(good) >= 1, 'omitnan');
         end
         % Net surface-height trend (m/yr): + for accumulation, - for ablation.
         trend = surfaceTrend(data.Time(good), s(good));
      end

      % Channel completeness: tice10m (L3 standardized 10 m) + raw tice string.
      tice10m_present = ismember("tice10m", vn);
      tice_string_present = any(startsWith(vn, "tice") & ~strcmp(vn, "tice10m"));
   catch err
      warning("profilePromiceData:siteFailed", "%s: %s", site, err.message);
   end
end

function tr = surfaceTrend(t, s)
   %SURFACETREND Linear m/yr slope of the surface series.
   tr = NaN;
   if numel(t) < 30, return, end
   x = years(t - t(1));
   ok = isfinite(x) & isfinite(s);
   if nnz(ok) < 30, return, end
   p = polyfit(x(ok), s(ok), 1);
   tr = p(1);
end

%% --------------------------------------------------------------- SUMup
function [has_rho, has_temp, smb_mean] = profileSumup(point, radius_km)
   %PROFILESUMUP Co-located SUMup density / temperature presence + mean SMB.
   has_rho = false; has_temp = false; smb_mean = NaN;
   try
      [obs, ~] = icemodel.verification.setup.buildSumupObservations(point, ...
         radius_km=radius_km);
      has_rho = isfield(obs, 'density') && ~isempty(obs.density);
      has_temp = isfield(obs, 'subsurface_temperature') ...
         && ~isempty(obs.subsurface_temperature);
      if isfield(obs, 'accumulation') && ~isempty(obs.accumulation) ...
            && ismember("smb", string(obs.accumulation.Properties.VariableNames))
         smb_mean = mean(obs.accumulation.smb, 'omitnan');
      end
   catch err
      warning("profileSumup:failed", "%s", err.message);
   end
end

%% ----------------------------------------------------------- classifiers
function c = classifyAccum(smb_mean, trend)
   %CLASSIFYACCUM Bin accumulation magnitude low/med/high.
   %
   % Prefer SUMup SMB (m w.e./yr) when available; else fall back to the PROMICE
   % net surface-height trend (m/yr, a snow-thickness proxy). Bins follow the
   % GrIS accumulation gradient: < 0.3 low, 0.3-0.6 med, > 0.6 high m w.e./yr.
   v = smb_mean;
   if ~isfinite(v) && isfinite(trend) && trend > 0
      v = 0.35 * trend;   % crude snow-thickness -> w.e. (rho~350) proxy
   end
   if ~isfinite(v)
      c = "unknown"; return
   end
   if v < 0.30
      c = "low";
   elseif v <= 0.60
      c = "med";
   else
      c = "high";
   end
end

function tf = classifyMelt(zone, elev)
   %CLASSIFYMELT Melt presence at the site.
   %
   % The PERCOLATION facies is, by the authoritative MODIS+SUMup classification,
   % a melt-bearing zone (meltwater percolates into firn) - so percolation sites
   % carry melt by definition. ABLATION sites melt by definition. The
   % ACCUMULATION facies is the high cold interior; melt there is marginal and
   % keyed on elevation: above the dry-snow line (~2500 m, Benson's dry-snow
   % facies) summer melt is negligible (DRY FIRN); below it some melt occurs.
   % The surface-height seasonal drawdown is not a reliable melt proxy on these
   % gap-bridged GC-Net records, so we use the facies + elevation classification
   % consistent with promicesiteinfo rather than the noisy amplitude.
   if any(zone == ["ablation", "percolation"])
      tf = true; return
   end
   if zone == "accumulation"
      tf = elev < 2500;   % below the dry-snow line -> marginal melt
      return
   end
   tf = false;
end

function r = classifyRegime(zone, melt)
   %CLASSIFYREGIME Firn-physics verification regime(s) for a site.
   %
   % (a) INFILTRATION : melt present, percolation zone (meltwater percolates and
   %                    refreezes in cold firn). This is the dominant GC-Net
   %                    accumulation-transect regime.
   % (c) DRY_FIRN     : little/no melt, high cold interior (densification +
   %                    thermal only; no meltwater).
   % (b) AQUIFER      : a PERENNIAL firn-aquifer water table requires
   %                    high-accumulation BURYING summer meltwater faster than
   %                    the cold-content can refreeze it - a SE-Greenland
   %                    phenomenon (Forster 2014; Miege 2016; FA-13/FA-15 cores
   %                    in the Helheim/SE catchment, ~1.5-2 m w.e./yr accum +
   %                    intense melt). The interior GC-Net set does NOT satisfy
   %                    this: its accumulation is moderate and its melt is
   %                    refreeze-limited (cold-firn percolation, not a water
   %                    table). We DO NOT auto-assign aquifer from a high SMB
   %                    proxy alone - that would falsely claim coverage. The
   %                    aquifer regime is reported as a GAP in the markdown.
   %
   % NOTE on the SMB-magnitude column: a "high" accum bin in SE/CE high-snowfall
   % sectors is NOT sufficient for the aquifer class. Aquifer formation is a
   % buried-meltwater phenomenon specific to the warm, very-high-accumulation SE
   % margin; none of these interior anchors sit there. We therefore key the
   % melt-bearing firn regime on the percolation facies and flag aquifer
   % separately rather than emit a false aquifer label.
   if zone == "ablation"
      r = "ablation";
   elseif melt
      % Melt-bearing accumulation-area firn: cold-firn percolation/infiltration.
      r = "infiltration";
   else
      r = "dry_firn";
   end
end

%% ----------------------------------------------------------- geography
function sectors = openSectors(file)
   %OPENSECTORS Load the BASINS mask + 2-D LAT/LON for nearest-cell lookup.
   sectors.B = ncread(file, 'BASINS');
   sectors.LAT = ncread(file, 'LAT');
   sectors.LON = ncread(file, 'LON');
end

function name = sampleSector(sectors, lat, lon)
   %SAMPLESECTOR Mouginot/Zwally basin name at the nearest grid cell.
   d2 = (sectors.LAT - lat).^2 + (sectors.LON - lon).^2;
   [~, idx] = min(d2(:));
   b = round(sectors.B(idx));
   name = basinName(b);
end

function s = basinName(b)
   %BASINNAME Map the integer basin id to a named GrIS sector.
   %
   % Greenland_sectors.nc carries the Zwally/Mouginot 8-basin numbering on the
   % RACMO FGRN055 grid. The integer ids index the standard clockwise-from-NW
   % drainage sectors.
   switch b
      case 1, s = "NW";
      case 2, s = "CW";
      case 3, s = "SW";
      case 4, s = "SE";
      case 5, s = "CE";
      case 6, s = "NE";
      case 7, s = "N";
      case 8, s = "NO";
      otherwise, s = "off-ice";
   end
end

%% ----------------------------------------------------------- RCM overlap
function s = overlapYears(req_years, leg)
   %OVERLAPYEARS Overlap-window string "YYYY-YYYY (N)" for one RCM leg.
   if isempty(leg.years)
      s = sprintf("absent");
      return
   end
   ov = intersect(req_years, leg.years);
   if isempty(ov)
      s = sprintf("none (disk %d-%d)", leg.year_min, leg.year_max);
   else
      s = sprintf("%d-%d", min(ov), max(ov));
   end
end

%% ----------------------------------------------------------- candidates
function sites = defaultCandidates()
   %DEFAULTCANDIDATES The 13 firn/accumulation candidate sites.
   %
   % The PROMICE/GC-Net readme Table 1 Accumulation-type sites that ship a firn
   % surface-height + subsurface-temperature record, plus KAN_U (site_type
   % Accumulation, surface_zone percolation - the KAN-transect firn anchor).
   sites = ["CEN", "CP1", "DY2", "EGP", "HUM", "KAN_U", ...
      "NAE", "NAU", "NEM", "NSE", "SDL", "SDM", "TUN"];
end

%% ----------------------------------------------------------- markdown
function writeMarkdown(mdfile, T, cov, radius_km)
   %WRITEMARKDOWN Render the coverage table + recommendation + gaps to markdown.
   fid = fopen(mdfile, 'w');
   c = onCleanup(@() fclose(fid));

   fprintf(fid, "# Firn-zone evaluation-data coverage analysis\n\n");
   fprintf(fid, "_Generated by `test/interactive/firn_coverage_analysis.m` ");
   fprintf(fid, "(%s). Gitignored interactive artifact._\n\n", ...
      string(datetime('now')));

   fprintf(fid, "Candidate set: the 13 PROMICE/GC-Net **Accumulation**-type ");
   fprintf(fid, "firn sites (+ KAN_U, the KAN-transect firn anchor). ");
   fprintf(fid, "SUMup co-location radius = %g km.\n\n", radius_km);

   % --- Coverage table ---
   fprintf(fid, "## Coverage table\n\n");
   cols = T.Properties.VariableNames;
   fprintf(fid, "| %s |\n", strjoin(cols, " | "));
   fprintf(fid, "|%s\n", strjoin(repmat("---", 1, numel(cols)), "|") + "|");
   for r = 1:height(T)
      vals = strings(1, numel(cols));
      for k = 1:numel(cols)
         vals(k) = cell2str(T{r, k});
      end
      fprintf(fid, "| %s |\n", strjoin(vals, " | "));
   end

   % --- RCM coverage note ---
   fprintf(fid, "\n## RCM on-disk coverage (per-leg)\n\n");
   fprintf(fid, "- MAR: %s\n", legSpan(cov.mar));
   fprintf(fid, "- MERRA-2: %s\n", legSpan(cov.merra));
   fprintf(fid, "- RACMO: %s\n\n", legSpan(cov.racmo));
   fprintf(fid, "PROMICE met + eval are NEVER gated by RCM coverage ");
   fprintf(fid, "(per-leg decoupled, decision #15). A site stages its FULL ");
   fprintf(fid, "PROMICE record; RCM legs are optional skip-with-reason.\n\n");

   % --- Regime coverage summary ---
   fprintf(fid, "## Firn-physics regime coverage\n\n");
   regimes = ["infiltration", "dry_firn", "aquifer"];
   labels = ["INFILTRATION/percolation", "DRY-FIRN densification/thermal", ...
      "FIRN AQUIFER"];
   for k = 1:numel(regimes)
      hits = T.site(contains(T.regime, regimes(k)));
      fprintf(fid, "- **%s**: %s\n", labels(k), ...
         ternary(isempty(hits), "_NONE (GAP)_", strjoin(hits, ", ")));
   end

   fprintf(fid, "\n## Flagged regime GAP: firn aquifer\n\n");
   fprintf(fid, "The interior GC-Net Accumulation transect does NOT sample ");
   fprintf(fid, "the **firn aquifer** regime (high-accumulation + high-melt, ");
   fprintf(fid, "SE Greenland; Forster 2014, Miege 2016). Aquifer coverage ");
   fprintf(fid, "requires the SUMup SE-Greenland firn cores (FA-13/FA-15, ");
   fprintf(fid, "Helheim catchment) or a dedicated FA core - none of the 13 ");
   fprintf(fid, "candidates co-locate there. This is an HONEST GAP: do not ");
   fprintf(fid, "force an interior site into the aquifer class. Stage a SUMup ");
   fprintf(fid, "SE point as a standalone aquifer case to close it.\n");
end

function s = legSpan(leg)
   if isempty(leg.years)
      s = "absent on disk";
   else
      s = sprintf("%d-%d", leg.year_min, leg.year_max);
   end
end

function s = cell2str(v)
   if iscell(v), v = v{1}; end
   if isstring(v) || ischar(v)
      s = string(v);
   elseif islogical(v)
      s = ternary(v, "yes", "no");
   elseif isnumeric(v)
      if isnan(v), s = "-"; else, s = string(v); end
   else
      s = string(v);
   end
end

function out = ternary(cond, a, b)
   if cond, out = a; else, out = b; end
end
