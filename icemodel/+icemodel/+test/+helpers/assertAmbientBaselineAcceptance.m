function assertAmbientBaselineAcceptance(ambient_stable, anchor_valid, ...
      anchor_ratio, accept_ambient_drift)
   %ASSERTAMBIENTBASELINEACCEPTANCE Validate the final baseline anchor.
   %
   %  icemodel.test.helpers.assertAmbientBaselineAcceptance( ...
   %     ambient_stable, anchor_valid, anchor_ratio, accept_ambient_drift)
   %
   % A release can accept a finite anchor drift when ACCEPT_AMBIENT_DRIFT is
   % true. An invalid anchor sample remains an error.
   %
   % Inputs
   %  ambient_stable      - final ambient-anchor verdict
   %  anchor_valid        - validity of the final anchor samples
   %  anchor_ratio        - final anchor median / first case median
   %  accept_ambient_drift - true to accept finite drifted release timings
   %
   % See also: icemodel.test.helpers.ambientAnchorVerdict,
   %  build_perf_baseline

   arguments
      ambient_stable (1, 1) logical
      anchor_valid (1, 1) logical
      anchor_ratio (1, 1) double
      accept_ambient_drift (1, 1) logical
   end

   % The override applies only to a valid, finite anchor measurement.
   if ambient_stable || (accept_ambient_drift && anchor_valid ...
         && isfinite(anchor_ratio))
      return
   end

   error('icemodel:test:perf:ambientDrift', ...
      ['Ambient conditions shifted during the baseline build, or the ' ...
      'anchor remeasurement was invalid (anchor ratio %.3f).'], anchor_ratio)
end
