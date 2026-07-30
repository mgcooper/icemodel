function estimate = applyDonorTransfer(transfer, times, x_donor, kwargs)
   %APPLYDONORTRANSFER Apply a fitted donor transfer to donor samples.
   %
   %  estimate = icemodel.forcing.reconstruct.applyDonorTransfer( ...
   %     transfer, times, x_donor)
   %
   % Role
   %  Application half of the donor-transfer tier: shifts the donor by the
   %  fitted lag, evaluates the per-season model, and enforces the
   %  policy's extrapolation limit — donor samples beyond the fitted
   %  donor range plus the allowed fraction of its span produce NaN
   %  rather than an extrapolated value (POLICY B4), so the next tier or
   %  missingness handles them honestly.
   %
   % Name-value
   %  max_extrapolation_fraction : how far beyond the fitted donor range
   %     the transfer applies, as a fraction of that range (default from
   %     icemodel.forcing.reconstruct.setopts).
   %
   % Returns
   %  estimate : double column of transferred values (NaN where the donor
   %     is missing or out of the admitted range).
   %
   % See also: icemodel.forcing.reconstruct.fitDonorTransfer

   arguments
      transfer (1, 1) struct
      times datetime
      x_donor (:, 1) double
      kwargs.max_extrapolation_fraction (1, 1) double ...
         {mustBeNonnegative} = ...
         icemodel.forcing.reconstruct.setopts().max_extrapolation_fraction
   end

   donor = x_donor;
   target_toa = ones(size(donor));
   is_csi = isfield(transfer, 'transfer_space') ...
      && string(transfer.transfer_space) == "clear_sky_index";
   if is_csi
      donor = icemodel.forcing.reconstruct.clearSkyIndex( ...
         times, donor, transfer.donor_location, ...
         toa_dark_wm2=transfer.toa_dark_wm2);
      target_toa = icemodel.forcing.reconstruct.toaIrradiance(times, ...
         transfer.target_location.lat_wgs84, ...
         transfer.target_location.lon_wgs84);
   end
   if numel(times) < 2
      % A one-sample axis cannot represent a nonzero fitted time shift.
      if transfer.lag_hours ~= 0
         donor(:) = NaN;
      end
      lag_steps = 0;
   else
      dt_hours = hours(median(diff(times)));
      lag_steps = round(transfer.lag_hours / dt_hours);
   end
   if abs(lag_steps) >= numel(donor) && lag_steps ~= 0
      % A transfer fitted on a longer axis cannot retain any donor sample.
      donor(:) = NaN;
   elseif lag_steps > 0
      donor = [nan(lag_steps, 1); donor(1:end - lag_steps)];
   elseif lag_steps < 0
      donor = [donor(1 - lag_steps:end); nan(-lag_steps, 1)];
   end

   % The extrapolation window widens the fitted range by the allowed
   % fraction of its span on each side; beyond it the transfer refuses
   % rather than extrapolates.
   span = diff(transfer.donor_range);
   lo = transfer.donor_range(1) ...
      - kwargs.max_extrapolation_fraction * span;
   hi = transfer.donor_range(2) ...
      + kwargs.max_extrapolation_fraction * span;
   admissible = isfinite(donor) & donor >= lo & donor <= hi;

   season = icemodel.forcing.reconstruct.seasonOf(times);
   estimate = nan(numel(donor), 1);
   for name = ["DJF", "MAM", "JJA", "SON"]
      in_season = admissible & season == name;
      if ~any(in_season)
         continue
      end
      model = transfer.models.(char(name));
      if model.kind == "linear"
         estimate(in_season) = model.slope * donor(in_season) + ...
            model.intercept;
      else
         % Monotone piecewise evaluation with end-segment continuation
         % inside the admitted range only.
         estimate(in_season) = interp1(model.breakpoints, ...
            model.node_values, donor(in_season), 'linear', 'extrap');
      end
   end
   if is_csi
      estimate = estimate .* target_toa;
      estimate(target_toa < transfer.toa_dark_wm2) = NaN;
   end
end
