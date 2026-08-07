function [per_case, aggregate, diagnostics] = ablationPerformanceMetrics( ...
      site_results, kwargs)
   %ABLATIONPERFORMANCEMETRICS Score every modeled ablation diagnostic.
   %
   %  [per_case, aggregate, diagnostics] = ...
   %     icemodel.verification.ablationPerformanceMetrics(site_results)
   %
   % SITE_RESULTS is the saved per-site-year result array. Each completed case
   % supplies an ALIGNED timetable whose rows are the eligible observation
   % timestamps inside that case's selected evaluation window, so model and
   % observation values are already on identical timestamps and no metric
   % interpolates or reindexes anything.
   %
   % PER_CASE has one row per completed case and scored diagnostic. AGGREGATE
   % collapses those rows to one row per diagnostic. DIAGNOSTICS is the ordered
   % list of scored model channels with their display labels, which is the one
   % source consumers use so a channel cannot be scored under two names.
   %
   % Density is an explicit scoring axis, not a fixed choice. Every diagnostic
   % is scored once per density in policy.effective_density_kg_m3, against the
   % measured lowering converted at that density. Reporting one density would
   % hide that the sign of the model-observation difference reverses across the
   % band, so the caller compares the whole set and the report ranks at the
   % policy reference density.
   %
   % Metrics per case and diagnostic. Only the endpoint error is a statement
   % about the cumulative curves. Every distributional metric is computed on
   % per-step INCREMENTS, because two cumulative series that both rise all
   % season are trivially correlated: their residuals are strongly
   % autocorrelated and a Nash-Sutcliffe efficiency computed on them is
   % inflated toward 1 regardless of whether the model gets the ablation rate
   % right. Differencing tests the rate, which is the physical claim.
   %   endpoint_error_mwe   - cumulative model minus observation at the last row
   %   cumulative_rmse_mwe  - RMSE of the cumulative curves, retained only so
   %                          the curve-level view stays reportable
   %   bias_mwe             - mean signed per-step increment error
   %   mae_mwe              - mean absolute per-step increment error
   %   rmse_mwe             - RMSE of the per-step increment error
   %   nse                  - Nash-Sutcliffe efficiency on the increments
   %
   % Increments are formed only between consecutive aligned rows separated by
   % exactly one model output step. A pair spanning an observation gap would
   % otherwise lump several hours into one increment and understate the rate.
   %
   % Cases that are incomplete, or whose aligned payload cannot support a
   % metric, are retained with a stated reason rather than dropped silently.
   %
   % See also: icemodel.verification.compareAblation

   arguments
      site_results struct
      kwargs.output_step (1, 1) duration = hours(1)
      kwargs.densities (1, :) double = double.empty(1, 0)
   end

   % Converting observed lowering to water equivalent needs a density, and the
   % choice matters more for increments than for endpoints: an hour of lowering
   % may be intact ice or porous weathering crust, and the density scales the
   % observed increments while leaving the modeled ones untouched. Score every
   % density in the policy band so that sensitivity is visible instead of
   % buried in a single intact-ice assumption.
   densities = kwargs.densities;
   if isempty(densities)
      policy = icemodel.verification.namelists.promiceAblationPolicy();
      densities = policy.effective_density_kg_m3;
   end
   ro_liq = icemodel.physicalConstant('ro_liq');

   diagnostics = performanceDiagnostics();
   n_diagnostics = numel(diagnostics);
   n_results = numel(site_results);

   % Preallocate the full case-by-diagnostic grid so no row grows in a loop and
   % every excluded case still produces an auditable row.
   n_densities = numel(densities);
   n_rows = n_results * n_diagnostics * n_densities;
   case_id = strings(n_rows, 1);
   site_id = strings(n_rows, 1);
   year = zeros(n_rows, 1);
   diagnostic = strings(n_rows, 1);
   label = strings(n_rows, 1);
   density_kg_m3 = zeros(n_rows, 1);
   scored = false(n_rows, 1);
   reason = strings(n_rows, 1);
   n_samples = zeros(n_rows, 1);
   n_increments = zeros(n_rows, 1);
   window_start = NaT(n_rows, 1, 'TimeZone', 'UTC');
   window_end = NaT(n_rows, 1, 'TimeZone', 'UTC');
   observation_endpoint_mwe = NaN(n_rows, 1);
   model_endpoint_mwe = NaN(n_rows, 1);
   endpoint_error_mwe = NaN(n_rows, 1);
   cumulative_rmse_mwe = NaN(n_rows, 1);
   bias_mwe = NaN(n_rows, 1);
   mae_mwe = NaN(n_rows, 1);
   rmse_mwe = NaN(n_rows, 1);
   nse = NaN(n_rows, 1);

   row = 0;
   for k = 1:n_results
      result = site_results(k);
      [aligned, case_reason] = alignedPayload(result);
      for d = 1:n_diagnostics
       for q = 1:n_densities
         row = row + 1;
         case_id(row) = string(result.case_id);
         site_id(row) = string(result.site_id);
         year(row) = double(result.year);
         diagnostic(row) = diagnostics(d).name;
         label(row) = diagnostics(d).label;
         density_kg_m3(row) = densities(q);
         if case_reason ~= ""
            reason(row) = case_reason;
            continue
         end
         if ~ismember(diagnostics(d).name, ...
               string(aligned.Properties.VariableNames))
            reason(row) = "aligned payload lacks " + diagnostics(d).name;
            continue
         end
         observation = ...
            aligned.observation_lowering_m * densities(q) / ro_liq;
         model = aligned.(char(diagnostics(d).name));
         finite = isfinite(observation) & isfinite(model);
         if nnz(finite) < 2
            reason(row) = "fewer than two finite aligned pairs";
            continue
         end
         time = aligned.Time(finite);
         observation = observation(finite);
         model = model(finite);
         residual = model - observation;

         % Increment pairs must be exactly one output step apart so a gap in
         % the observation record cannot masquerade as a large hourly rate.
         adjacent = diff(time) == kwargs.output_step;
         d_observation = diff(observation);
         d_model = diff(model);
         d_observation = d_observation(adjacent);
         d_model = d_model(adjacent);
         d_residual = d_model - d_observation;

         scored(row) = true;
         n_samples(row) = numel(observation);
         n_increments(row) = numel(d_residual);
         window_start(row) = time(1);
         window_end(row) = time(end);
         observation_endpoint_mwe(row) = observation(end);
         model_endpoint_mwe(row) = model(end);
         endpoint_error_mwe(row) = model(end) - observation(end);
         cumulative_rmse_mwe(row) = sqrt(mean(residual .^ 2));
         if isempty(d_residual)
            continue
         end
         bias_mwe(row) = mean(d_residual);
         mae_mwe(row) = mean(abs(d_residual));
         rmse_mwe(row) = sqrt(mean(d_residual .^ 2));
         nse(row) = nashSutcliffe(d_observation, d_model);
       end
      end
   end

   per_case = table(case_id, site_id, year, diagnostic, label, ...
      density_kg_m3, scored, ...
      reason, n_samples, n_increments, window_start, window_end, ...
      observation_endpoint_mwe, model_endpoint_mwe, endpoint_error_mwe, ...
      cumulative_rmse_mwe, bias_mwe, mae_mwe, rmse_mwe, nse);
   aggregate = aggregateMetrics(per_case, diagnostics, densities);
end

function diagnostics = performanceDiagnostics()
   %PERFORMANCEDIAGNOSTICS Return the scored model channels in display order.

   % Every scored channel is a cumulative metres-water-equivalent series that a
   % reader could plausibly compare against measured lowering. They answer
   % different physical questions, which is exactly what the scoring exposes.
   % model_surface_mass_loss_mwe is deliberately NOT scored. It measures the
   % mass merge_layers destroys, not ablation: a merge keeps the mean of the
   % two cells in one cell, so removing an empty top cell still deletes mass
   % from the cell below. See icemodel-4nv.
   diagnostics = struct( ...
      'name', { ...
      "model_melt_mwe", ...
      "model_runoff_mwe", ...
      "model_ablation_proxy_mwe", ...
      "model_solid_loss_mwe"}, ...
      'label', { ...
      "gross melt", ...
      "runoff", ...
      "runoff + vapor loss", ...
      "signed net solid balance"});
end

function [aligned, reason] = alignedPayload(result)
   %ALIGNEDPAYLOAD Return one case's aligned comparison rows or a stated reason.

   aligned = timetable.empty(0, 0);
   reason = "";
   if ~isfield(result, 'status') || string(result.status) ~= "completed"
      reason = "case status is not completed";
      return
   end
   if ~isfield(result, 'aligned') || ~istimetable(result.aligned) ...
         || isempty(result.aligned)
      reason = "case has no aligned comparison payload";
      return
   end
   if ~ismember("observation_intact_mwe", ...
         string(result.aligned.Properties.VariableNames))
      reason = "aligned payload lacks observation_intact_mwe";
      return
   end
   aligned = result.aligned;
end

function value = nashSutcliffe(observation, model)
   %NASHSUTCLIFFE Return the efficiency of MODEL against OBSERVATION.

   % A constant observed series has no variance to explain, so the efficiency
   % is undefined rather than zero or one.
   denominator = sum((observation - mean(observation)) .^ 2);
   if denominator <= 0
      value = NaN;
      return
   end
   value = 1 - sum((model - observation) .^ 2) / denominator;
end

function aggregate = aggregateMetrics(per_case, diagnostics, densities)
   %AGGREGATEMETRICS Collapse the per-case rows to one row per diagnostic
   %and density.

   % The practical-agreement tolerance is a reporting policy, not a constant of
   % the method, so it comes from the namelist that the report also reads.
   policy = icemodel.verification.namelists.promiceAblationPolicy();
   tolerance = policy.scientific.endpoint_tolerance_fraction;

   n = numel(diagnostics) * numel(densities);
   diagnostic = strings(n, 1);
   label = strings(n, 1);
   density_kg_m3 = zeros(n, 1);
   n_scored = zeros(n, 1);
   n_excluded = zeros(n, 1);
   median_endpoint_error_mwe = NaN(n, 1);
   mean_endpoint_error_mwe = NaN(n, 1);
   mean_bias_mwe = NaN(n, 1);
   mean_mae_mwe = NaN(n, 1);
   mean_rmse_mwe = NaN(n, 1);
   pooled_rmse_mwe = NaN(n, 1);
   median_nse = NaN(n, 1);
   n_within_tolerance = zeros(n, 1);

   k = 0;
   for d = 1:numel(diagnostics)
    for q = 1:numel(densities)
      k = k + 1;
      rows = per_case.diagnostic == diagnostics(d).name ...
         & per_case.density_kg_m3 == densities(q);
      scored = rows & per_case.scored;
      diagnostic(k) = diagnostics(d).name;
      label(k) = diagnostics(d).label;
      density_kg_m3(k) = densities(q);
      n_scored(k) = nnz(scored);
      n_excluded(k) = nnz(rows & ~per_case.scored);
      if n_scored(k) == 0
         continue
      end
      endpoint = per_case.endpoint_error_mwe(scored);
      median_endpoint_error_mwe(k) = median(endpoint);
      mean_endpoint_error_mwe(k) = mean(endpoint);
      mean_bias_mwe(k) = mean(per_case.bias_mwe(scored), 'omitnan');
      mean_mae_mwe(k) = mean(per_case.mae_mwe(scored), 'omitnan');
      mean_rmse_mwe(k) = mean(per_case.rmse_mwe(scored), 'omitnan');

      % Pool the squared increment error by increment count so long site-years
      % are not weighted the same as short ones in the headline number.
      weights = per_case.n_increments(scored);
      usable = weights > 0 & isfinite(per_case.rmse_mwe(scored));
      case_rmse = per_case.rmse_mwe(scored);
      if any(usable)
         pooled_rmse_mwe(k) = sqrt(sum(weights(usable) ...
            .* case_rmse(usable) .^ 2) / sum(weights(usable)));
      end
      median_nse(k) = median(per_case.nse(scored), 'omitnan');

      % A tolerance count communicates practical agreement more directly than a
      % mean error that opposite-signed cases can cancel.
      observed = abs(per_case.observation_endpoint_mwe(scored));
      n_within_tolerance(k) = nnz(abs(endpoint) <= tolerance * observed);
    end
   end

   aggregate = table(diagnostic, label, density_kg_m3, n_scored, n_excluded, ...
      median_endpoint_error_mwe, mean_endpoint_error_mwe, mean_bias_mwe, ...
      mean_mae_mwe, mean_rmse_mwe, pooled_rmse_mwe, median_nse, ...
      n_within_tolerance);
end
