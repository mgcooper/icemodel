function bands = solarElevationBands()
   %SOLARELEVATIONBANDS Single source of the swd solar-elevation banding.
   %
   %  bands = icemodel.forcing.reconstruct.solarElevationBands()
   %
   % Role
   %  Dedicated single-source function (like physicalBounds and
   %  bucketEdges) for every solar-elevation threshold the swd science
   %  path consumes, so the darkness pre-pass and the elevation-binned
   %  proxy calibration can never drift apart.
   %
   % Returns
   %  bands : struct with fields
   %     civil_twilight_deg : -6, the standard civil-twilight boundary.
   %        Below it the sky contributes no meaningful shortwave, so a
   %        missing swd sample is a KNOWN zero (POLICY B2 refined by
   %        D-28). Between -6 and 0 deg stations measure 8-28 W m^-2 of
   %        real diffuse twilight light, so those samples are never
   %        hard-zeroed and instead reach the normal fill tiers.
   %     calibration_bin_edges_deg : [-90 -6 0 5 15 30 90], the ratio-bin
   %        edges for the D-28 elevation-binned swd proxy calibration. A
   %        single per-season ratio mixes regimes with opposite biases —
   %        RCM proxies fill 4-20x LOW near solar midnight (kanu hour 1:
   %        observed median 17 W m^-2 vs filled 1) and 1.5-2x HIGH on the
   %        morning shoulder (kanl hour 7: observed 81 vs MAR 171) — so
   %        the bins isolate deep darkness [-90, -6), the twilight band
   %        [-6, 0), the low-sun shoulders [0, 5) and [5, 15), mid sun
   %        [15, 30), and high sun [30, 90]. The twilight edge reuses
   %        civil_twilight_deg's value so both passes share one boundary.
   %     min_bin_samples : 200, the overlap samples a season x elevation
   %        bin needs before its own ratio is trusted; thinner bins
   %        inherit the season's single ratio. 200 hourly samples is a
   %        few weeks of a season's ~2200 hours: enough to stabilize a
   %        median ratio while still letting the narrow twilight bins
   %        qualify.
   %     toa_ceiling_floor_wm2 : 5, the minimum candidate ceiling when
   %        1.05 times geometric TOA is smaller (sensor night noise).
   %     toa_ceiling_multiplier : 1.05, the candidate allowance above
   %        geometric top-of-atmosphere irradiance.
   %     solar_constant_wm2 : 1361, the extraterrestrial direct-normal
   %        ceiling used by toaIrradiance and the D-32 flux-anchor guard.
   %     twilight_ceiling_wm2 : 50, the candidate ceiling within civil
   %        twilight, where diffuse incident light is real despite zero
   %        geometric TOA (POLICY D-28).
   %     horizon_deg : 0, the upper edge of civil twilight.
   %     twilight_climatology_min_support : 30, the minimum local
   %        day-of-year/posting observations for the D-44 final-edge fill.
   %
   % See also: icemodel.forcing.reconstruct.fitProxyCalibration,
   %  icemodel.forcing.reconstruct.applyProxyCalibration,
   %  icemodel.forcing.reconstruct.reconstructSeries

   % One struct so callers destructure only the values they need; the
   % numbers live here and nowhere else (single-source rule).
   civil_twilight_deg = -6;
   horizon_deg = 0;
   bands = struct( ...
      'civil_twilight_deg', civil_twilight_deg, ...
      'horizon_deg', horizon_deg, ...
      'calibration_bin_edges_deg', ...
      [-90, civil_twilight_deg, horizon_deg, 5, 15, 30, 90], ...
      'min_bin_samples', 200, ...
      'toa_ceiling_floor_wm2', 5, ...
      'toa_ceiling_multiplier', 1.05, ...
      'solar_constant_wm2', 1361, ...
      'twilight_ceiling_wm2', 50, ...
      'twilight_climatology_min_support', 30);
end
