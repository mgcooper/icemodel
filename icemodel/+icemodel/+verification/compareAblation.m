function [summary, aligned, diagnostics, policy] = compareAblation( ...
      observations, model, kwargs)
   %COMPAREABLATION Compare PROMICE lowering with modeled solid-ice loss.
   %
   %  [summary, aligned, diagnostics, policy] = ...
   %     icemodel.verification.compareAblation(observations, model)
   %  [...] = icemodel.verification.compareAblation(..., ...,
   %     window_start=t0, window_end=t1,
   %     endpoint_deficit_kg_m2=[lower central upper])
   %
   % OBSERVATIONS is a staged PROMICE timetable or a bundle whose data field is
   % that timetable. MODEL is a postprocessed diagnostic IceModel timetable or
   % an equivalent data bundle. PROMICE ablation is geometric positive-down
   % lowering [m]. The modeled primary is positive solid-ice loss
   % -(mass_budget_phase_solid_mwe + mass_budget_vapor_solid_mwe) [m w.e.], the
   % solid_balance term from
   % icemodel.verification.helpers.ablationLedgerIncrements. Remeshing is never
   % added to it.
   % Each diagnostic model row is a forcing-step ledger stamped at the interval
   % start, so observations at t0 and t1 consume model rows on [t0,t1).
   % Paired values additionally require finite snow depth at or below the fixed
   % exposed-ice threshold; trace-snow rows remain inside the cumulative span
   % but are censored from alignment.
   % ENDPOINT_DEFICIT_KG_M2 supplies caller-provided sensitivity scenarios for
   % D(t1)-D(t0). This comparator cannot validate their observation provenance,
   % so they remain visible but never govern classification. This function
   % also converts geometric observations over the fixed sensitivity range
   % 600--900 kg m^-3, from a porous weathering crust to intact ice. The
   % 600 kg m^-3 value is the porous endpoint. It is not an intact-ice
   % density, and it is not always the numeric lower bound, because signed
   % lowering can be negative.

   arguments
      observations
      model
      kwargs.window_start (1, 1) datetime = NaT
      kwargs.window_end (1, 1) datetime = NaT
      kwargs.endpoint_deficit_kg_m2 double = double.empty(1, 0)
   end

   policy = icemodel.verification.namelists.promiceAblationPolicy();
   obs_tt = comparisonTimetable(observations, "observation");
   model_tt = comparisonTimetable(model, "model");

   % Validate the staged and diagnostic schemas before selecting a window so a
   % missing physical term cannot be mistaken for unavailable temporal support.
   requireVariables(obs_tt, policy.required_observation_fields, ...
      'icemodel:verification:compareAblation:missingObservationField');
   % Require only the channels this comparison reads.
   % required_model_fields is derived live from the output namelists, so
   % requiring all of it rejects a cohort saved before any later channel was
   % appended, over channels the comparison never touches. Compatibility, not
   % equality (see icemodel.verification.helpers.validateAblationModelSchema,
   % which applies the same rule to the report).
   requireVariables(model_tt, ...
      icemodel.verification.namelists.ablationReportChannels('ledger'), ...
      'icemodel:verification:compareAblation:missingModelField');

   % The finiteness and sign checks below run over the channels this cohort
   % carries. A channel it predates cannot be non-finite in it.
   checked_fields = intersect(policy.required_model_fields, ...
      string(model_tt.Properties.VariableNames), 'stable');
   obs_tt = normalizeTimetable(obs_tt, "observation");
   model_tt = normalizeTimetable(model_tt, "model");

   endpoint_deficit = kwargs.endpoint_deficit_kg_m2(:);
   if any(~isfinite(endpoint_deficit))
      error('icemodel:verification:compareAblation:badEndpointDeficit', ...
         'endpoint_deficit_kg_m2 must contain only finite scenarios')
   end

   % The requested bounds are explicit in every result; omitted bounds mean the
   % full overlap, not a hidden calendar or site-specific study window.
   [window_start, window_end] = comparisonWindow(obs_tt.Time, model_tt.Time, ...
      kwargs.window_start, kwargs.window_end);
   [common_time, obs_idx, model_idx] = intersect( ...
      obs_tt.Time, model_tt.Time, 'stable');
   in_window = common_time >= window_start & common_time <= window_end;

   obs_values = numericValues(obs_tt, obs_idx, ...
      policy.required_observation_fields, "observation");
   model_values = numericValues(model_tt, model_idx, ...
      checked_fields, "model");
   observation_fields = policy.required_observation_fields;
   snow_values = obs_values(:, ...
      observation_fields == policy.snow_variable);
   % classifyObservationSupport applies the PROMICE flag rules and returns the
   % row-shaped support masks used below.
   support = icemodel.verification.helpers.classifyObservationSupport( ...
      obs_values, observation_fields, policy.observation_field, policy);
   metadata_flagged = support.metadata_flagged;
   obs_finite = support.target_finite;
   model_finite = all(isfinite(model_values), 2);
   quality_finite = support.quality_finite;
   gap_bridged = support.gap_flagged;
   station_transition = support.station_transition;
   unresolved_step = support.unresolved_step;

   % Direct support here does not require a finite target: obs_finite is
   % tested separately in the eligibility mask below, and the exclusion
   % diagnostics need to tell a missing observation from a flagged one.
   direct = quality_finite & support.direct_flags_zero;
   [ice_exposed, snow_censored, unknown_snow] = ...
      icemodel.verification.helpers.classifySnowDepth( ...
      snow_values, policy.ice_exposure_threshold_m);
   eligible = in_window & obs_finite & model_finite & direct & ice_exposed;

   if nnz(eligible) < 2
      error('icemodel:verification:compareAblation:insufficientSupport', ...
         ['the requested window has %d direct finite exact-common samples; ' ...
         'at least two are required'], nnz(eligible))
   end

   eligible_time = common_time(eligible);
   actual_start = eligible_time(1);
   actual_end = eligible_time(end);

   % The comparator uses cumulative levels at exact direct endpoints, so an
   % ordinary surface-height gap may remain between plotted markers. A station
   % transition or unresolved/unknown step changes the cumulative datum and
   % therefore invalidates any window that crosses it.
   observation_window = obs_tt(obs_tt.Time >= actual_start ...
      & obs_tt.Time <= actual_end, :);
   observation_values = numericValues(observation_window, ...
      (1:height(observation_window))', ...
      policy.required_observation_fields, "observation");
   window_support = ...
      icemodel.verification.helpers.classifyObservationSupport( ...
      observation_values, observation_fields, policy.observation_field, ...
      policy);
   datum_intact = window_support.datum_intact;
   if any(~datum_intact)
      error('icemodel:verification:compareAblation:observationDatumBreak', ...
         ['comparison endpoints cannot cross a station transition or an ' ...
         'unresolved/unknown observation step'])
   end

   ledger = model_tt(model_tt.Time >= actual_start ...
      & model_tt.Time < actual_end, :);
   expected_step = seconds(policy.model_output_cadence_seconds);
   model_support_time = [ledger.Time; actual_end];
   model_day_start = dateshift(model_support_time, 'start', 'day');
   model_elapsed_seconds = seconds(model_support_time - model_day_start);
   if any(mod(model_elapsed_seconds, ...
         policy.model_output_cadence_seconds) ~= 0)
      error('icemodel:verification:compareAblation:offGridModelTime', ...
         'model diagnostic rows must lie on native UTC output boundaries')
   end
   if isempty(ledger) || ledger.Time(1) ~= actual_start ...
         || ledger.Time(end) + expected_step ~= actual_end ...
         || any(diff(ledger.Time) ~= expected_step)
      error('icemodel:verification:compareAblation:noncontiguousModelSupport', ...
         ['model diagnostics must contain every %g-second interval-start ' ...
         'row on [t0,t1)'], policy.model_output_cadence_seconds)
   end
   ledger_values = numericValues(ledger, (1:height(ledger))', ...
      checked_fields, "model");
   if any(~isfinite(ledger_values), 'all')
      error('icemodel:verification:compareAblation:nonfiniteModelWindow', ...
         'every required diagnostic must be finite throughout the model window')
   end

   % Independently rebase observation geometry and model cumulative loss at the
   % first direct common point. Prefix zero is the state at t0; each subsequent
   % value includes rows strictly before its state-observation timestamp.
   increments = ...
      icemodel.verification.helpers.ablationLedgerIncrements(ledger);
   physical_increment = increments.solid_balance;
   model_prefix = [0; cumsum(physical_increment)];
   [~, eligible_model_idx] = ismember(eligible_time, model_tt.Time);
   start_model_idx = eligible_model_idx(1);
   model_aligned = model_prefix(eligible_model_idx - start_model_idx + 1);
   cumulative_idx = eligible_model_idx - 1;
   if any(cumulative_idx < 1) || any( ...
         model_tt.Time(cumulative_idx) + expected_step ~= eligible_time)
      error('icemodel:verification:compareAblation:missingBoundaryState', ...
         ['cumulative diagnostics require the preceding hourly endpoint ' ...
         'state at every eligible observation time'])
   end
   cumulative_fields = string( ...
      icemodel.namelists.cumulativeoutputs());
   cumulative_values = numericValues(model_tt, cumulative_idx, ...
      cumulative_fields, "model cumulative boundary");
   if any(~isfinite(cumulative_values), 'all')
      error('icemodel:verification:compareAblation:nonfiniteBoundaryState', ...
         'every cumulative diagnostic boundary state must be finite')
   end
   model_melt = cumulative_values(:, cumulative_fields == "melt");
   model_runoff = cumulative_values(:, cumulative_fields == "runoff");
   model_freeze = cumulative_values(:, cumulative_fields == "freeze");
   % dlayer is the mass remeshing removed, in metres water equivalent. A merge
   % keeps exactly half the pair's water and the fixed-depth grid clones a
   % bottom cell, so this is a remeshing term rather than an ablation flux.
   % See icemodel-4nv.
   model_layer_change = cumulative_values(:, cumulative_fields == "dlayer");
   model_melt = model_melt - model_melt(1);
   model_runoff = model_runoff - model_runoff(1);
   model_freeze = model_freeze - model_freeze(1);
   model_layer_change = model_layer_change - model_layer_change(1);
   model_surface_loss = [0; cumsum(increments.surface_loss)];
   model_surface_loss = model_surface_loss( ...
      eligible_model_idx - start_model_idx + 1);

   model_vapor_loss = [0; cumsum(increments.solid_vapor_loss)];
   model_vapor_loss = model_vapor_loss( ...
      eligible_model_idx - start_model_idx + 1);
   model_ablation_proxy = model_runoff + model_vapor_loss;
   obs_raw = obs_tt.(policy.observation_field)(obs_idx(eligible));
   obs_lowering = obs_raw - obs_raw(1);
   [Ls, Lv, ro_ice, ro_liq] = icemodel.physicalConstant( ...
      'Ls', 'Lv', 'ro_ice', 'ro_liq');
   obs_intact_mwe = obs_lowering * ro_ice / ro_liq;
   density = policy.effective_density_kg_m3(:);
   obs_density_mwe = obs_lowering .* density.' / ro_liq;
   % Pointwise extrema keep the lower and upper field names correct when
   % signed lowering is negative and the two density endpoints swap order.
   obs_density_lower_mwe = min(obs_density_mwe, [], 2);
   obs_density_upper_mwe = max(obs_density_mwe, [], 2);

   % The reference density is the representative bubbly near-surface ice value
   % drawn inside the band, so a reader sees one best estimate rather than only
   % an envelope.
   obs_reference_mwe = obs_lowering ...
      * policy.effective_density_reference_kg_m3 / ro_liq;

   aligned = timetable(obs_raw, obs_lowering, obs_intact_mwe, ...
      obs_density_lower_mwe, obs_density_upper_mwe, obs_reference_mwe, ...
      model_melt, ...
      model_runoff, model_freeze, model_aligned, model_layer_change, ...
      model_surface_loss, model_ablation_proxy, ...
      model_aligned - obs_intact_mwe, ...
      'RowTimes', eligible_time, ...
      'VariableNames', {'observation_ablation_m', ...
      'observation_lowering_m', 'observation_intact_mwe', ...
      'observation_lower_mwe', 'observation_upper_mwe', ...
      'observation_reference_mwe', ...
      'model_melt_mwe', 'model_runoff_mwe', 'model_freeze_mwe', ...
      'model_solid_loss_mwe', 'model_layer_change_mwe', ...
      'model_surface_mass_loss_mwe', 'model_ablation_proxy_mwe', ...
      'model_minus_observation_mwe'});

   % Evaluate the checkpoint, phase, mixed-latent-heat, and B/O remesh
   % identities over exactly the interval represented by the rebased series.
   diagnostics.identities = closureIdentities( ...
      ledger, policy, Ls, Lv, ro_liq);
   diagnostics.physical_comparable = all(diagnostics.identities.passed);

   a_model = model_aligned(end);
   h_obs = obs_lowering(end);
   a_obs = obs_intact_mwe(end);
   signal = max([abs(a_model), abs(a_obs), ...
      policy.scientific.signal_floor_mwe]);

   % Retain signed and non-cancelling ratios for each accounting channel. A
   % large positive and negative exchange must not disappear by cancellation.
   diagnostics.materiality = materialityDiagnostics( ...
      ledger, endpoint_deficit, signal, policy, Ls, ro_liq);
   [diagnostics.scenarios, non_identifiable_reasons] = ...
      scenarioDiagnostics(a_model, a_obs, endpoint_deficit, ...
      diagnostics.materiality, ledger, policy, Ls, ro_liq);

   % Effective-density conversions are always labeled sensitivities and never
   % promoted to rigorous endpoint bounds or identifiability evidence.
   density_mwe = h_obs .* density ./ ro_liq;
   diagnostics.effective_density = table(density, density_mwe, ...
      false(numel(density), 1), ...
      'VariableNames', {'effective_density_kg_m3', ...
      'observation_sensitivity_mwe', 'rigorous_bound'});
   diagnostics.excluded = exclusionCounts(in_window, eligible, obs_finite, ...
      model_finite, quality_finite, gap_bridged, station_transition, ...
      unresolved_step, metadata_flagged, unknown_snow, snow_censored);
   diagnostics.top_deletion = struct( ...
      'count', sum(ledger.mass_budget_top_deletion_count), ...
      'height_m', sum(ledger.mass_budget_top_deletion_height_m), ...
      'export_solid_mwe', sum(ledger.mass_budget_top_export_solid_mwe), ...
      'export_liquid_mwe', sum(ledger.mass_budget_top_export_liquid_mwe), ...
      'interior_merge_count', sum(ledger.mass_budget_interior_merge_count), ...
      'role', "count and height are quantized grid geometry and are never " + ...
      "mass; export_* is the mass a surface removal actually removed");

   % Closure failure takes precedence over scientific classification. Otherwise
   % a credible material-accounting alternative that changes the result makes
   % the comparison non-identifiable instead of forcing a directional claim.
   difference = a_model - a_obs;
   base_classification = classifyDifference(difference, a_model, a_obs, policy);
   if ~diagnostics.physical_comparable
      classification = "not_physically_comparable";
   elseif ~isempty(non_identifiable_reasons)
      classification = "non_identifiable";
   else
      classification = base_classification;
   end

   summary = struct( ...
      'classification', classification, ...
      'base_classification', base_classification, ...
      'physical_comparable', diagnostics.physical_comparable, ...
      'requested_window_start', window_start, ...
      'requested_window_end', window_end, ...
      'window_start', actual_start, ...
      'window_end', actual_end, ...
      'common_sample_count', nnz(in_window), ...
      'eligible_sample_count', nnz(eligible), ...
      'observation_lowering_m', h_obs, ...
      'observation_intact_mwe', a_obs, ...
      'model_solid_loss_mwe', a_model, ...
      'model_minus_observation_mwe', difference, ...
      'relative_difference', abs(difference) / signal, ...
      'endpoint_deficit_scenario_count', numel(endpoint_deficit), ...
      'non_identifiable_reasons', non_identifiable_reasons, ...
      'policy_version', policy.version);
end

function tt = comparisonTimetable(value, role)
   %COMPARISONTIMETABLE Unwrap a timetable or verification data bundle.
   if istimetable(value)
      tt = value;
   elseif isstruct(value) && isfield(value, 'data') ...
         && istimetable(value.data)
      tt = value.data;
   else
      error('icemodel:verification:compareAblation:badInput', ...
         '%s input must be a timetable or a struct with timetable data', role)
   end
end


function requireVariables(tt, required, error_id)
   %REQUIREVARIABLES Error if a required variable is missing.
   missing = setdiff(required, string(tt.Properties.VariableNames), 'stable');
   if ~isempty(missing)
      error(error_id, 'missing required field(s): %s', strjoin(missing, ', '))
   end
end

function tt = normalizeTimetable(tt, role)
   %NORMALIZETIMETABLE Normalize timestamps and reject ambiguous duplicates.
   if ~isdatetime(tt.Properties.RowTimes) || any(isnat(tt.Properties.RowTimes))
      error('icemodel:verification:compareAblation:badTime', ...
         '%s row times must be finite datetimes', role)
   end
   tt.Properties.RowTimes = icemodel.verification.setup.ensureUtc( ...
      tt.Properties.RowTimes);
   tt = sortrows(tt);
   if numel(unique(tt.Time)) ~= height(tt)
      error('icemodel:verification:compareAblation:duplicateTime', ...
         '%s timetable must have unique timestamps', role)
   end
end

function [window_start, window_end] = comparisonWindow( ...
      obs_time, model_time, requested_start, requested_end)
   %COMPARISONWINDOW Resolve explicit bounds against the available overlap.
   overlap_start = max(obs_time(1), model_time(1));
   overlap_end = min(obs_time(end), model_time(end));
   if isnat(requested_start)
      window_start = overlap_start;
   else
      window_start = icemodel.verification.setup.ensureUtc(requested_start);
   end
   if isnat(requested_end)
      window_end = overlap_end;
   else
      window_end = icemodel.verification.setup.ensureUtc(requested_end);
   end
   if window_start > window_end || window_end < overlap_start ...
         || window_start > overlap_end
      error('icemodel:verification:compareAblation:invalidWindow', ...
         'the requested comparison window does not overlap both inputs')
   end
end

function values = numericValues(tt, rows, names, role)
   %NUMERICVALUES Extract numeric comparison fields with one stable failure.
   try
      values = double(tt{rows, cellstr(names)});
   catch
      error('icemodel:verification:compareAblation:nonNumericField', ...
         '%s comparison fields must be numeric or logical', role)
   end
end

function identities = closureIdentities(ledger, policy, Ls, Lv, ro_liq)
   %CLOSUREIDENTITIES Evaluate mass, vapor-energy, and remesh identities.
   p_s = ledger.mass_budget_phase_solid_mwe;
   p_l = ledger.mass_budget_phase_liquid_mwe;
   v_s = ledger.mass_budget_vapor_solid_mwe;
   v_l = ledger.mass_budget_vapor_liquid_mwe;
   r_s = ledger.mass_budget_remesh_solid_mwe;
   r_l = ledger.mass_budget_remesh_liquid_mwe;
   d_s = ledger.mass_budget_solid_end_mwe(end) ...
      - ledger.mass_budget_solid_start_mwe(1);
   d_l = ledger.mass_budget_liquid_end_mwe(end) ...
      - ledger.mass_budget_liquid_start_mwe(1);
   step_d_s = ledger.mass_budget_solid_end_mwe ...
      - ledger.mass_budget_solid_start_mwe;
   step_d_l = ledger.mass_budget_liquid_end_mwe ...
      - ledger.mass_budget_liquid_start_mwe;

   % Coupled interior transport moves storage between the phases without
   % changing their sum. The per-phase closures need its increments.
   x_s = ledger.mass_budget_vapor_redistribution_solid_mwe;
   x_l = ledger.mass_budget_vapor_redistribution_liquid_mwe;

   % Two checkpoint mass identities retain the solid and liquid storage
   % views. Their sum reproduces the total-storage identity, so no separate
   % total row is needed. The magnitude term sums the component channels
   % separately: opposing components inside one row must widen the
   % relative tolerance, not cancel it down to the absolute floor.
   rows(1) = identityRow("solid_storage", "mwe", ...
      d_s, step_d_s, p_s + v_s + r_s + x_s, ...
      abs(p_s) + abs(v_s) + abs(r_s) + abs(x_s), policy, 1);
   rows(2) = identityRow("liquid_storage", "mwe", ...
      d_l, step_d_l, p_l + v_l + r_l + x_l, ...
      abs(p_l) + abs(v_l) + abs(r_l) + abs(x_l), policy, 1);

   % Treat solid phase change as the increment and its liquid counterpart as
   % the opposing transfer. Their combined per-row residual must close; equal
   % and opposite residuals in different forcing rows cannot cancel.
   rows(3) = identityRow("phase_mass", "mwe", sum(p_s), p_s, -p_l, ...
      abs(p_l), policy, 1);

   % Vapor energy keeps solid sublimation/deposition on Ls and liquid vapor or
   % overflow on Lv; signed unapplied energy remains an explicit remainder.
   %
   % Coupled interior transport runs after the vapor budget takes its
   % storage baseline, so no interior term reaches these channels and the
   % identity needs no redistribution correction. The redistribution
   % channels stand alone as the transport's own accounting.
   potential = ledger.mass_budget_vapor_potential_j_m2;
   overflow = ledger.mass_budget_condensation_overflow_mwe;
   unapplied = ledger.mass_budget_unapplied_vapor_j_m2;
   accepted = ro_liq .* (Ls .* v_s + Lv .* v_l + Lv .* overflow) ...
      + unapplied;
   accepted_mag = ro_liq .* (Ls .* abs(v_s) + Lv .* abs(v_l) ...
      + Lv .* abs(overflow)) + abs(unapplied);
   energy_scale = ro_liq * Lv;
   rows(4) = identityRow("vapor_energy", "J m-2", ...
      sum(potential), potential, accepted, accepted_mag, policy, ...
      energy_scale);

   % B/O decomposition is checked for solid mass, the phase the ablation
   % comparator scores. The liquid clone and export carry no standing
   % channels; remesh_liquid stays visible through the storage identity
   % and the materiality diagnostics.
   b_s = ledger.mass_budget_cloned_bottom_solid_mwe;
   o_s = ledger.mass_budget_merge_export_solid_mwe;
   rows(5) = identityRow("remesh_solid", "mwe", ...
      sum(r_s), r_s, b_s - o_s, abs(b_s) + abs(o_s), policy, 1);
   identities = struct2table(rows);
end

function row = identityRow(name, units, window_delta, step_delta, step_flux, ...
      step_flux_mag, policy, scale)
   %IDENTITYROW Apply exact window and per-forcing-step closure acceptance.
   %
   % STEP_FLUX is the per-row signed net the identity closes against.
   % STEP_FLUX_MAG sums the component channels' magnitudes for the same
   % rows, so opposing components inside one row widen the relative
   % tolerance instead of cancelling it down to the absolute floor.
   residual = window_delta - sum(step_flux);
   normalization = max(policy.numerical.q_floor_mwe * scale, ...
      abs(window_delta) + sum(step_flux_mag));
   tolerance = policy.numerical.atol_mwe * scale ...
      + policy.numerical.rtol * normalization;

   % A signed window residual can vanish even when individual forcing rows do
   % not close. Test every row with the same absolute-plus-relative policy.
   step_residual = step_delta - step_flux;
   step_normalization = max(policy.numerical.q_floor_mwe * scale, ...
      abs(step_delta) + step_flux_mag);
   step_tolerance = policy.numerical.atol_mwe * scale ...
      + policy.numerical.rtol .* step_normalization;
   window_passed = abs(residual) <= tolerance;
   step_passed = all(abs(step_residual) <= step_tolerance);
   row = struct('identity', name, 'units', units, 'residual', residual, ...
      'normalization', normalization, 'tolerance', tolerance, ...
      'window_passed', window_passed, 'step_passed', step_passed, ...
      'failed_step_count', nnz(abs(step_residual) > step_tolerance), ...
      'passed', window_passed && step_passed);
end

function materiality = materialityDiagnostics(ledger, endpoint_deficit, ...
      signal, policy, Ls, ro_liq)
   %MATERIALITYDIAGNOSTICS Preserve signed and non-cancelling channel ratios.
   step_d_l = ledger.mass_budget_liquid_end_mwe ...
      - ledger.mass_budget_liquid_start_mwe;
   rows = [ ...
      materialityRow("endpoint_liquid_storage", step_d_l, signal, policy); ...
      materialityRow("remesh_solid", ...
      ledger.mass_budget_remesh_solid_mwe, signal, policy); ...
      materialityRow("remesh_liquid", ...
      ledger.mass_budget_remesh_liquid_mwe, signal, policy); ...
      materialityRow("merge_delete_solid", ...
      ledger.mass_budget_merge_export_solid_mwe, signal, policy); ...
      materialityRow("cloned_bottom_solid", ...
      ledger.mass_budget_cloned_bottom_solid_mwe, signal, policy); ...
      materialityRow("condensation_overflow", ...
      ledger.mass_budget_condensation_overflow_mwe, signal, policy); ...
      materialityRow("unapplied_vapor_solid_equivalent", ...
      ledger.mass_budget_unapplied_vapor_j_m2 ./ (ro_liq * Ls), ...
      signal, policy)];
   % The caller's deficit list has a known length, so size the array once
   % rather than growing it per iteration.
   n_deficit = numel(endpoint_deficit);
   if n_deficit > 0
      n_fixed = numel(rows);
      rows(n_fixed + n_deficit) = rows(n_fixed);
      for n = 1:n_deficit
         rows(n_fixed + n) = materialityRow( ...
            "endpoint_deficit_" + n, endpoint_deficit(n) / ro_liq, ...
            signal, policy);
      end
   end
   materiality = struct2table(rows);
end

function row = materialityRow(channel, values, signal, policy)
   %MATERIALITYROW Score one channel without cross-channel cancellation.
   %
   % The gross is the sum of the per-row magnitudes, so a large positive and
   % negative exchange in different rows cannot disappear by cancellation.
   signed_net = sum(values);
   gross = sum(abs(values));
   row = struct('channel', channel, 'signed_net_mwe', signed_net, ...
      'gross_mwe', gross, 'signed_ratio', signed_net / signal, ...
      'gross_ratio', gross / signal, ...
      'material', gross / signal > policy.scientific.materiality_limit);
end

function [scenarios, reasons] = scenarioDiagnostics(a_model, a_obs, ...
      endpoint_deficit, materiality, ledger, policy, Ls, ro_liq)
   %SCENARIODIAGNOSTICS Test endpoint sensitivities and accounting alternatives.
   base_difference = a_model - a_obs;
   base_class = classifyDifference(base_difference, a_model, a_obs, policy);
   rows = scenarioRow("central_intact_ice", "central", a_model, a_obs, ...
      base_difference, base_class, policy, true, true);

   % The public numeric input has no validated observation-artifact provenance.
   % Preserve each sensitivity value, but mark it unverified and non-governing.
   n_deficit = numel(endpoint_deficit);
   if n_deficit > 0
      n_fixed = numel(rows);
      rows(n_fixed + n_deficit) = rows(n_fixed);
      for n = 1:n_deficit
         scenario_obs = a_obs + endpoint_deficit(n) / ro_liq;
         rows(n_fixed + n) = scenarioRow("endpoint_deficit_" + n, ...
            "endpoint", a_model, scenario_obs, base_difference, ...
            base_class, policy, ...
            isMaterial(materiality, "endpoint_deficit_" + n), false);
      end
   end

   % Only demonstrated solid-equivalent alternatives can adjust the solid-loss
   % comparator. Liquid storage, liquid domain exchange, and overflow remain
   % explicit materiality diagnostics because no observation operator maps
   % them one-for-one onto current-window solid loss.
   scenario_names = ["remesh_solid", ...
      "merge_delete_solid", "cloned_bottom_solid", ...
      "unapplied_vapor_solid_equivalent"];
   adjustments = [-sum(ledger.mass_budget_remesh_solid_mwe), ...
      sum(ledger.mass_budget_merge_export_solid_mwe), ...
      -sum(ledger.mass_budget_cloned_bottom_solid_mwe), ...
      -sum(ledger.mass_budget_unapplied_vapor_j_m2) / (ro_liq * Ls)];
   n_fixed = numel(rows);
   rows(n_fixed + numel(scenario_names)) = rows(n_fixed);
   for n = 1:numel(scenario_names)
      rows(n_fixed + n) = scenarioRow(scenario_names(n), "accounting", ...
         a_model + adjustments(n), a_obs, base_difference, base_class, ...
         policy, isMaterial(materiality, scenario_names(n)), true);
   end
   scenarios = struct2table(rows);

   % Only credible material accounting can govern classification. Endpoint
   % sensitivities remain diagnostic until a validated runner-derived evidence
   % path exists. A caller-supplied value carries no provenance of its own.
   changed = scenarios.changes_sign | scenarios.changes_classification;
   governs = scenarios.role == "accounting" & scenarios.material;
   reasons = scenarios.scenario(changed & governs & scenarios.credible);
end

function row = scenarioRow(name, role, a_model, a_obs, base_difference, ...
      base_class, policy, material, credible)
   %SCENARIOROW Classify one observation or accounting alternative.
   difference = a_model - a_obs;
   classification = classifyDifference(difference, a_model, a_obs, policy);
   signal = max([abs(a_model), abs(a_obs), ...
      policy.scientific.signal_floor_mwe]);
   row = struct('scenario', name, 'role', role, ...
      'model_mwe', a_model, 'observation_mwe', a_obs, ...
      'difference_mwe', difference, ...
      'relative_difference', abs(difference) / signal, ...
      'classification', classification, ...
      'changes_sign', sign(difference) ~= sign(base_difference), ...
      'changes_classification', classification ~= base_class, ...
      'material', material, 'credible', credible);
end

function classification = classifyDifference(difference, a_model, a_obs, policy)
   %CLASSIFYDIFFERENCE Assign a neutral bias class at the fixed threshold.
   signal = max([abs(a_model), abs(a_obs), ...
      policy.scientific.signal_floor_mwe]);
   if abs(difference) / signal <= policy.scientific.materiality_limit
      classification = "within_materiality";
   elseif difference > 0
      classification = "model_high";
   else
      classification = "model_low";
   end
end

function tf = isMaterial(materiality, channel)
   %ISMATERIAL Return the non-cancelling gate for one named channel.
   idx = materiality.channel == channel;
   tf = any(materiality.material(idx));
end

function counts = exclusionCounts(in_window, eligible, obs_finite, ...
      model_finite, quality_finite, gap_bridged, station_transition, ...
      unresolved_step, metadata_flagged, unknown_snow, snow_censored)
   %EXCLUSIONCOUNTS Report every primary-support exclusion independently.
   %
   % metadata_flagged arrives as a logical mask from
   % classifyObservationSupport, so it is already finite-and-nonzero and must
   % not be re-tested with > 0. A malformed negative posting therefore counts
   % as flagged.
   counts = struct( ...
      'gap_bridged', nnz(in_window & gap_bridged), ...
      'station_transition', nnz(in_window & station_transition), ...
      'unresolved_step', nnz(in_window & unresolved_step), ...
      'step_correctable_but_unresolved', ...
      nnz(in_window & unresolved_step & metadata_flagged), ...
      'nonfinite_observation', nnz(in_window & ~obs_finite), ...
      'unknown_quality_flag', nnz(in_window & ~quality_finite), ...
      'unknown_snow_depth', nnz(in_window & unknown_snow), ...
      'snow_censored', nnz(in_window & snow_censored), ...
      'nonfinite_model', nnz(in_window & ~model_finite), ...
      'total_unique_excluded', nnz(in_window & ~eligible));
end
