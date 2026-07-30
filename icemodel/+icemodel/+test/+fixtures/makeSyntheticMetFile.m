function [met, filepath] = makeSyntheticMetFile(simyear, kwargs)
   %MAKESYNTHETICMETFILE Build a simple synthetic forcing timetable for tests.
   %
   %  met = icemodel.test.fixtures.makeSyntheticMetFile(simyear)
   %  [met, filepath] = icemodel.test.fixtures.makeSyntheticMetFile(simyear, ...
   %     sitename="kanm", forcings="kanm", metdir="/tmp")

   arguments
      simyear (1, 1) double {mustBeInteger, mustBePositive}
      kwargs.sitename (1, :) char = 'kanm'
      kwargs.forcings (1, :) char = 'kanm'
      kwargs.nsteps (1, 1) double {mustBeInteger, mustBePositive} = 48
      kwargs.dt_seconds (1, 1) double {mustBePositive} = 3600
      kwargs.include_modis (1, 1) logical = false
      kwargs.modis_varname (1, :) char = 'MODIS'
      kwargs.snow_depth = []
      % Optional precipitation overrides: full-length vectors for the
      % rainf/snowf split channels and a ppt replacement (e.g. all-NaN to
      % exercise loadmet's runtime total-precip derivation).
      kwargs.ppt = []
      kwargs.rainf = []
      kwargs.snowf = []
      kwargs.metdir (1, :) char = ''
   end

   % Build one simple, smooth forcing year with a stable diurnal cycle.
   dt = seconds(kwargs.dt_seconds);
   t0 = datetime(simyear, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
   Time = transpose(t0 + dt * (0:kwargs.nsteps-1));
   tod_hours = hour(Time) + minute(Time) / 60 + second(Time) / 3600;

   tair = 263.15 + 4.0 * sin(2 * pi * tod_hours / 24) + 0.25 * (simyear - 2015);
   swd = max(0, 250 * sin(pi * (tod_hours - 6) / 12));
   lwd = 220 + 20 * cos(2 * pi * tod_hours / 24);
   albedo = 0.62 + 0.01 * cos(2 * pi * tod_hours / 24);
   wspd = 4.0 + 0.5 * sin(2 * pi * tod_hours / 24);
   rh = 75 + 5 * cos(2 * pi * tod_hours / 24);
   psfc = 78000 + zeros(kwargs.nsteps, 1);
   ppt = zeros(kwargs.nsteps, 1);

   met = timetable(Time, tair, swd, lwd, albedo, wspd, rh, psfc, ppt);
   if ~isempty(kwargs.snow_depth)
      snow_depth = kwargs.snow_depth;
      if isscalar(snow_depth)
         snow_depth = repmat(snow_depth, kwargs.nsteps, 1);
      else
         snow_depth = snow_depth(:);
      end
      if numel(snow_depth) ~= kwargs.nsteps
         error('snow_depth must be scalar or length nsteps')
      end
      met.snowd = snow_depth;
   end
   % Precipitation overrides install after the base table so tests can
   % model a split-only source (rainf/snowf present, ppt placeholder).
   if ~isempty(kwargs.ppt)
      met.ppt = checkLength(kwargs.ppt, kwargs.nsteps, 'ppt');
   end
   if ~isempty(kwargs.rainf)
      met.rainf = checkLength(kwargs.rainf, kwargs.nsteps, 'rainf');
   end
   if ~isempty(kwargs.snowf)
      met.snowf = checkLength(kwargs.snowf, kwargs.nsteps, 'snowf');
   end
   if kwargs.include_modis
      % Mirror the albedo shape for optional synthetic userdata tests.
      modis = min(max(albedo - 0.05, 0.05), 0.95);
      met.(kwargs.modis_varname) = modis;
   end

   % Optionally persist the timetable using the normal met-file naming.
   filepath = '';
   if ~isempty(kwargs.metdir)
      tag = formatTimestepTag(kwargs.dt_seconds);
      filepath = fullfile(kwargs.metdir, ...
         ['met_' kwargs.sitename '_' kwargs.forcings '_' int2str(simyear) ...
         '_' tag '.mat']);
      save(filepath, 'met');
   end
end

function x = checkLength(x, nsteps, name)
   %CHECKLENGTH Validate an override channel as a full-length column.
   x = x(:);
   if numel(x) ~= nsteps
      error('%s must be length nsteps', name)
   end
end

function tag = formatTimestepTag(dt_seconds)
   %FORMATTIMESTEPTAG Convert a synthetic timestep into the met filename tag.

   if mod(dt_seconds, 3600) == 0
      tag = [int2str(dt_seconds / 3600) 'hr'];
   elseif mod(dt_seconds, 60) == 0
      tag = [int2str(dt_seconds / 60) 'm'];
   else
      tag = [int2str(dt_seconds) 's'];
   end
end
