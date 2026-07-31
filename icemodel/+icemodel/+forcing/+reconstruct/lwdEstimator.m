function estimate = lwdEstimator(tair, rh)
   %LWDESTIMATOR Empirical downward-longwave candidate from temperature and RH.
   %
   %  estimate = icemodel.forcing.reconstruct.lwdEstimator(tair, rh)
   %
   % Role
   %  The policy's empirical lwd candidate (per-variable table, lwd proxy
   %  (c)): a clear-sky Brutsaert (1975) effective emissivity from vapor
   %  pressure and air temperature, epsilon = 1.24 (ea_hPa / T)^(1/7),
   %  applied as lwd = epsilon * sigma * T^4. This is a deliberately
   %  uncalibrated physical form: it competes in the empirical-estimator
   %  tier (POLICY B1) only AFTER passing
   %  through the same overlap calibration as the model proxies
   %  (fitProxyCalibration), which absorbs the cloud contribution the
   %  clear-sky form omits. Clean-room from the published formulation; no
   %  legacy coefficients are reused (none were recoverable — POLICY
   %  compatibility audit).
   %
   % Inputs
   %  tair : air temperature, K.
   %  rh : relative humidity with respect to water, percent.
   %
   % Returns
   %  estimate : downward longwave, W m-2 (NaN where inputs are missing).
   %
   % See also: icemodel.forcing.reconstruct.fitProxyCalibration

   arguments
      tair (:, 1) double
      rh (:, 1) double
   end

   % Saturation vapor pressure over water (Magnus form, hPa), scaled by RH
   % to the actual vapor pressure the emissivity keys on.
   es_hpa = 6.112 * exp(17.62 * (tair - 273.15) ./ (tair - 30.03));
   ea_hpa = max(0, rh / 100) .* es_hpa;

   % Brutsaert clear-sky effective emissivity and the graybody flux.
   sigma = 5.670374419e-8;
   emissivity = 1.24 * (ea_hpa ./ tair).^(1 / 7);
   estimate = emissivity * sigma .* tair.^4;
   estimate(~isfinite(tair) | ~isfinite(rh)) = NaN;
end
