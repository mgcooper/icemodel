function transfer = fitDonorTransfer(times, x_target, x_donor, channel, kwargs)
   %FITDONORTRANSFER Fit an overlap-calibrated donor-to-target transfer.
   %
   %  transfer = icemodel.forcing.reconstruct.fitDonorTransfer( ...
   %     times, x_target, x_donor, "tair", fit_years=split.years_selection)
   %
   % Role
   %  The donor-transfer model (POLICY B4): per-season linear
   %  regression target~donor on the concurrent overlap, with an optional
   %  monotone piecewise-linear adjustment whose knot count is the policy's
   %  validated spline hyperparameter (the clean-room stand-in for the
   %  paper's six-interval monotone spline — a recorded deviation), and a
   %  ±lag search applied only when it clearly improves the overlap
   %  correlation. Elevation adjustment happens BEFORE this fit via
   %  elevationAdjust. Fitting uses only the caller-supplied fit years so
   %  selection/evaluation hygiene is enforced by construction.
   %
   % Name-value
   %  fit_years : calendar years eligible for fitting (required — pass the
   %     split's selection years).
   %  knots : 0 for pure linear, or a knot count for the monotone
   %     piecewise adjustment (default 0; the candidate set is the
   %     Section-C knot_candidates parameter, currently [0 6]).
   %  lag_search : search donor lags within ±max_lag_hours (default true).
   %  max_lag_hours : lag search half-width (default 18, legacy precedent).
   %  min_lag_gain : correlation improvement required to adopt a nonzero
   %     lag (default 0.02); below it the lag stays 0 and is recorded so.
   %  min_overlap_hours : minimum concurrent finite overlap (default 8760,
   %     the policy's one-year requirement).
   %  target_location, donor_location : station points required for SWD;
   %     the transfer is fit in station-specific clear-sky-index space.
   %  toa_dark_wm2 : meaningful-sun threshold for the SWD CSI fit.
   %
   % Returns
   %  transfer : struct — channel, lag_hours, knots, per-season model
   %     (slope/intercept or breakpoint tables), fitted donor range,
   %     n_overlap, n_overlap_hours, overlap correlation before/after lag
   %     (diagnostics
   %     only — admission always comes from held-out gates).
   %
   % See also: icemodel.forcing.reconstruct.applyDonorTransfer,
   %  icemodel.forcing.reconstruct.elevationAdjust

   arguments
      times datetime
      x_target (:, 1) double
      x_donor (:, 1) double
      channel (1, 1) string
      kwargs.fit_years (1, :) double {mustBeInteger, mustBeNonempty}
      kwargs.knots (1, 1) double {mustBeInteger, mustBeNonnegative} = 0
      kwargs.lag_search (1, 1) logical = true
      kwargs.max_lag_hours (1, 1) double {mustBeNonnegative} = ...
         icemodel.forcing.reconstruct.setopts().max_lag_hours
      kwargs.min_lag_gain (1, 1) double {mustBeNonnegative} = ...
         icemodel.forcing.reconstruct.setopts().min_lag_gain
      kwargs.min_overlap_hours (1, 1) double {mustBePositive} = ...
         icemodel.forcing.reconstruct.setopts().min_overlap_hours
      kwargs.target_location (1, 1) struct = struct()
      kwargs.donor_location (1, 1) struct = struct()
      kwargs.toa_dark_wm2 (1, 1) double {mustBePositive} = ...
         icemodel.forcing.reconstruct.setopts().toa_dark_wm2
   end

   % A sampling interval cannot be inferred from one timestamp; report the
   % same policy failure used for every other under-constrained fit.
   if numel(times) < 2
      error('icemodel:reconstruct:fitDonorTransfer:insufficientOverlap', ...
         'at least two timestamps are required to fit %s', channel);
   end

   % SWD lag and regression operate on cloud transmissivity rather than raw
   % flux, so station solar geometry cannot masquerade as transfer skill.
   transfer_space = "native";
   if channel == "swd"
      x_target = icemodel.forcing.reconstruct.clearSkyIndex( ...
         times, x_target, kwargs.target_location, ...
         toa_dark_wm2=kwargs.toa_dark_wm2);
      x_donor = icemodel.forcing.reconstruct.clearSkyIndex( ...
         times, x_donor, kwargs.donor_location, ...
         toa_dark_wm2=kwargs.toa_dark_wm2);
      transfer_space = "clear_sky_index";
   end

   dt_target_hours = hours(median(diff(times)));
   in_years = ismember(year(times), kwargs.fit_years);

   % Lag search on the fit-year transfer space: shift the donor by whole steps
   % and keep the lag only when the correlation gain is material.
   max_lag_steps = round(kwargs.max_lag_hours / dt_target_hours);
   base_r = overlapCorrelation(x_target, x_donor, in_years, 0);
   best_r = base_r;
   best_lag = 0;
   if kwargs.lag_search
      for lag = [-max_lag_steps:-1, 1:max_lag_steps]
         r = overlapCorrelation(x_target, x_donor, in_years, lag);
         if r > best_r
            best_r = r;
            best_lag = lag;
         end
      end
      if best_r - base_r < kwargs.min_lag_gain
         best_lag = 0;
         best_r = base_r;
      end
   end
   donor = shiftSeries(x_donor, best_lag);

   % Integrate concurrent target intervals after the lag. Support-held coarse
   % donors repeat across target rows, so row count times TARGET cadence equals
   % elapsed support (e.g. four quarter-hours still count as one hour).
   overlap = in_years & isfinite(x_target) & isfinite(donor);
   n_overlap_hours = nnz(overlap) * dt_target_hours;
   if n_overlap_hours < kwargs.min_overlap_hours
      error('icemodel:reconstruct:fitDonorTransfer:insufficientOverlap', ...
         '%.0f overlap hours < required %.0f for %s', ...
         n_overlap_hours, kwargs.min_overlap_hours, channel);
   end

   % Per-season models; a season without enough overlap inherits the
   % all-season fit so application never silently mixes seasons.
   season = icemodel.forcing.reconstruct.seasonOf(times);
   all_model = fitOne(x_target(overlap), donor(overlap), kwargs.knots);
   models = struct();
   for name = ["DJF", "MAM", "JJA", "SON"]
      in_season = overlap & season == name;
      models.(char(name)) = all_model;
      if nnz(in_season) * dt_target_hours >= kwargs.min_overlap_hours / 4
         models.(char(name)) = fitOne(x_target(in_season), ...
            donor(in_season), kwargs.knots);
      end
   end

   transfer = struct('channel', channel, ...
      'lag_hours', best_lag * dt_target_hours, ...
      'knots', kwargs.knots, 'models', models, ...
      'fit_years', kwargs.fit_years, ...
      'donor_range', [min(donor(overlap)), max(donor(overlap))], ...
      'n_overlap', nnz(overlap), 'n_overlap_hours', n_overlap_hours, ...
      'overlap_r_base', base_r, 'overlap_r_lagged', best_r, ...
      'transfer_space', transfer_space, ...
      'target_location', kwargs.target_location, ...
      'donor_location', kwargs.donor_location, ...
      'toa_dark_wm2', kwargs.toa_dark_wm2);
end

function r = overlapCorrelation(x_target, x_donor, in_years, lag)
   %OVERLAPCORRELATION Pearson correlation of the lagged fit-year overlap.
   donor = shiftSeries(x_donor, lag);
   both = in_years & isfinite(x_target) & isfinite(donor);
   if nnz(both) < 3
      r = -Inf;
      return
   end
   c = corrcoef(x_target(both), donor(both));
   r = c(1, 2);
end

function x = shiftSeries(x, lag)
   %SHIFTSERIES Shift a series by whole steps, padding with NaN.
   if abs(lag) >= numel(x)
      % No source sample can remain on the requested axis.
      x(:) = NaN;
      return
   end
   if lag > 0
      x = [nan(lag, 1); x(1:end - lag)];
   elseif lag < 0
      x = [x(1 - lag:end); nan(-lag, 1)];
   end
end

function model = fitOne(y, d, knots)
   %FITONE Fit one linear or monotone piecewise-linear transfer.
   % The linear fit is ordinary least squares. The piecewise variant bins
   % the donor at equally spaced quantile breakpoints, takes the target
   % median per bin, and enforces monotonicity with a running maximum —
   % the clean-room stand-in for the paper's monotone spline whose
   % interval count the policy treats as a validated hyperparameter.
   if knots == 0
      coeffs = [d, ones(size(d))] \ y;
      model = struct('kind', "linear", 'slope', coeffs(1), ...
         'intercept', coeffs(2));
      return
   end
   edges = quantile(d, linspace(0, 1, knots + 1));
   edges = unique(edges);
   if numel(edges) < 2
      % A near-constant donor collapses every quantile to one value; a
      % scalar second argument would flip discretize into its N-bins form
      % and error, so fall back to the linear fit immediately.
      coeffs = [d, ones(size(d))] \ y;
      model = struct('kind', "linear", 'slope', coeffs(1), ...
         'intercept', coeffs(2));
      return
   end
   centers = (edges(1:end - 1) + edges(2:end)) / 2;
   node_values = nan(size(centers));
   bin = discretize(d, edges);
   for k = 1:numel(centers)
      node_values(k) = median(y(bin == k), 'omitnan');
   end
   keep = isfinite(node_values);
   if nnz(keep) < 2
      % Quantized or narrow donors can collapse the breakpoints; a
      % degenerate piecewise map falls back to the linear fit rather than
      % crash or fabricate a flat transfer.
      coeffs = [d, ones(size(d))] \ y;
      model = struct('kind', "linear", 'slope', coeffs(1), ...
         'intercept', coeffs(2));
      return
   end
   node_values = cummax(node_values(keep));
   model = struct('kind', "piecewise", ...
      'breakpoints', centers(keep), 'node_values', node_values);
end
