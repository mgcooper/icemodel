function [ice_exposed, snow_censored, unknown_snow] = ...
      classifySnowDepth(snow_depth_m, ice_exposure_threshold_m)
   %CLASSIFYSNOWDEPTH Classify snow support without accepting negative depth.
   %
   % [ice_exposed, snow_censored, unknown_snow] = ...
   %    icemodel.verification.helpers.classifySnowDepth( ...
   %    snow_depth_m, ice_exposure_threshold_m)
   %
   % Finite nonnegative depths at or below the threshold indicate exposed ice;
   % larger finite depths are snow-censored. Negative and nonfinite depths are
   % invalid observations and therefore have unknown snow support.

   arguments
      snow_depth_m double {mustBeReal}
      ice_exposure_threshold_m (1, 1) double ...
         {mustBeReal, mustBeFinite, mustBeNonnegative}
   end

   % The mutually exclusive masks retain the input shape for direct use with
   % timetable rows in readiness, comparison, runner, and report paths.
   valid_snow = isfinite(snow_depth_m) & snow_depth_m >= 0;
   ice_exposed = valid_snow ...
      & snow_depth_m <= ice_exposure_threshold_m;
   snow_censored = valid_snow ...
      & snow_depth_m > ice_exposure_threshold_m;
   unknown_snow = ~valid_snow;
end
