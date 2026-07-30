function rows = auditSegments(times, mask, channel, method, detail, kwargs)
   %AUDITSEGMENTS Build one audit row per contiguous selected segment.
   %
   %  rows = icemodel.forcing.reconstruct.auditSegments( ...
   %     times, mask, channel, method, detail)
   %
   % Returns a cell column whose rows match reconstructSeries.audit. A
   % disjoint selection is split so no audit row claims an unfilled span.
   % context_id joins a production segment to the exact fitted parameters
   % and held-out evidence in the persisted station plan.

   arguments
      times datetime
      mask (:, 1) logical
      channel (1, 1) string
      method (1, 1) string
      detail (1, 1) string
      kwargs.context_id (1, 1) string = method
   end

   if numel(times) ~= numel(mask)
      error('icemodel:reconstruct:auditSegments:sizeMismatch', ...
         'times and mask must share one sample axis');
   end

   % Empty selections need no cadence and produce the typed empty shape
   % callers can concatenate with other audit-cell columns.
   if ~any(mask)
      rows = cell(0, 1);
      return
   end
   if numel(times) < 2
      error('icemodel:reconstruct:auditSegments:missingCadence', ...
         'at least two timestamps are required to audit duration');
   end

   % Contiguous runs are the policy audit unit; duration uses the regular
   % posting support represented by each selected sample.
   dt_hours = hours(median(diff(times)));
   edges = diff([false; mask; false]);
   starts = find(edges == 1);
   stops = find(edges == -1) - 1;
   rows = cell(numel(starts), 1);
   for k = 1:numel(starts)
      n_samples = stops(k) - starts(k) + 1;
      rows{k} = {char(channel), times(starts(k)), times(stops(k)), ...
         n_samples * dt_hours, char(method), ...
         sprintf('%d samples; %s', n_samples, detail), ...
         char(kwargs.context_id)};
   end
end
