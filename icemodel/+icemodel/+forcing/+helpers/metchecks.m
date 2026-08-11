function [met, checks] = metchecks(met, kwargs)
   %METCHECKS Gap-fill and clamp met variables to physically valid ranges.
   %
   %  [met, checks] = icemodel.forcing.helpers.metchecks(met)
   %  [met, checks] = ... metchecks(met, fillgaps=false, clamp=false)
   %
   % METCHECKS is the standard QA/QC pass that every forcing builder applies
   % before it writes a met or Data timetable:
   %
   %  1. Counts NaN and complex-valued samples per variable. The counts
   %     are returned in CHECKS so builders can record data provenance.
   %  2. Gap-fills each variable by linear interpolation with
   %     nearest-value end fill. Wind direction (wdir) is circular, so this
   %     function fills it through its unit-vector components, not
   %     linearly. A linear fill corrupts a gap that crosses the 360/0
   %     wrap. The +forcing README describes the component fill.
   %  3. Clamps recognized variables to the legacy physical ranges:
   %
   %        albedo   [0.05, 0.98]   [-]
   %        rh       [5, 99.99]     [%]
   %        wspd     >= 0.1         [m s-1]
   %        wdir     wrapped to (0, 360]
   %        tsfc     <= 273.16 when in kelvin, <= 0 when in celsius
   %
   %     This function gap-fills an unrecognized variable but never clamps it.
   %
   % Inputs
   %  met - timetable holding any subset of the met-contract variables
   %        (see icemodel.forcing.helpers.metvariables)
   %
   % Outputs
   %  met    - the QA/QC'd timetable
   %  checks - struct of pre-fill diagnostics:
   %           checks.numnan     table of NaN counts per variable
   %           checks.numcomplex table of complex-sample counts per variable
   %
   % Legacy: reimplements runoff/functions/metchecks.m (retained, unchanged,
   % as the legacy reference). Differences: timetable-only input (the legacy
   % version also accepted structs/matrices); circular wdir gap-fill (the
   % legacy version linear-filled it); no plot option (inspect the returned
   % checks instead).
   %
   % See also: icemodel.forcing.helpers.metvariables,
   %  icemodel.forcing.helpers.validatemet, fillmissing

   arguments
      met timetable
      kwargs.fillgaps (1, 1) logical = true
      kwargs.clamp (1, 1) logical = true
   end

   varnames = string(met.Properties.VariableNames);
   n_vars = numel(varnames);

   % Count missing and complex samples per variable before any filling.
   numnan = zeros(1, n_vars);
   numcomplex = zeros(1, n_vars);
   for n = 1:n_vars
      data = met.(varnames(n));
      numnan(n) = sum(isnan(data), 'all');
      numcomplex(n) = sum(imag(data) ~= 0, 'all');
   end
   checks.numnan = array2table(numnan, 'VariableNames', varnames);
   checks.numcomplex = array2table(numcomplex, 'VariableNames', varnames);

   % Gap-fill. wdir is filled through its unit-vector components so fills
   % crossing the 360/0 wrap interpolate along the short arc.
   if kwargs.fillgaps
      for n = 1:n_vars
         varname = varnames(n);
         if varname == "wdir"
            met.wdir = fillWindDirection(met.wdir);
         else
            met.(varname) = fillmissing(met.(varname), 'linear', ...
               'EndValues', 'nearest');
         end
      end
   end

   % Clamp recognized variables to physical ranges.
   if kwargs.clamp
      if ismember("albedo", varnames)
         met.albedo = clampPresent(met.albedo, 0.05, 0.98);
      end
      if ismember("rh", varnames)
         met.rh = clampPresent(met.rh, 5, 99.99);
      end
      if ismember("wspd", varnames)
         % Apply the wind floor so a reconstructed met file is accepted at
         % runtime.
         wspd_bounds = ...
            icemodel.forcing.reconstruct.physicalBounds("wspd");
         met.wspd = clampPresent(met.wspd, wspd_bounds(1), Inf);
      end
      if ismember("wdir", varnames)
         met.wdir = wrapWindDirection(met.wdir);
      end
      if ismember("tsfc", varnames)
         % Surface temperature cannot exceed the melting point. Detect
         % kelvin vs celsius from the series magnitude (legacy
         % convention; glacier surfaces are never near 100 C).
         if min(met.tsfc, [], 'omitnan') > 100
            met.tsfc = clampPresent(met.tsfc, -Inf, ...
               icemodel.physicalConstant('Tf'));
         else
            met.tsfc = clampPresent(met.tsfc, -Inf, 0);
         end
      end
   end
end

%% Local functions
function data = clampPresent(data, lower, upper)
   %CLAMPPRESENT Clamp non-NaN samples while preserving NaN placeholders.
   valid = ~isnan(data);
   data(valid) = min(max(data(valid), lower), upper);
end

function wdir = fillWindDirection(wdir)
   %FILLWINDDIRECTION Gap-fill a circular direction series [degrees].
   x = fillmissing(sind(wdir), 'linear', 'EndValues', 'nearest');
   y = fillmissing(cosd(wdir), 'linear', 'EndValues', 'nearest');
   filled = wrapWindDirection(atan2d(x, y));
   wdir(isnan(wdir)) = filled(isnan(wdir));
end

function wdir = wrapWindDirection(wdir)
   %WRAPWINDDIRECTION Wrap direction to the meteorological range (0, 360].
   wdir = mod(wdir, 360);
   wdir(wdir == 0) = 360;
end
