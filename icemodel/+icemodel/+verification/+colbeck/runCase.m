function candidate = runCase(case_manifest, kwargs)
   %RUNCASE Build a Colbeck candidate bundle for the verification suite.
   %
   %  candidate = icemodel.verification.colbeck.runCase(case_manifest)
   %  candidate = icemodel.verification.colbeck.runCase(case_manifest, ...
   %     kind="numerical", experiment_names=["exp1", "exp2", "exp3"])
   %
   %  Single candidate-provider entry point for the Colbeck verification case.
   %  Dispatches internally on `kind`: the "numerical" branch runs a thin
   %  time-loop calling icemodel.column.infiltration once per step (all physics
   %  lives in the kernel); the "analytical" branch wraps
   %  icemodel.verification.colbeck.analyticalSolution. Both branches return the
   %  same experiment_bundle schema so the comparison driver treats them
   %  uniformly.
   %
   %  Inputs
   %    case_manifest      Resolved Colbeck case manifest (currently
   %                       informational; the canonical case definition is
   %                       loaded directly from caseDefinition).
   %  Name-value
   %    kind               "numerical" (default) | "analytical"
   %    experiment_names   String row vector. Default
   %                       ["exp1", "exp2", "exp3"].
   %
   %  Outputs
   %    candidate   Struct with fields:
   %       format        "experiment_bundle"
   %       experiments   Struct keyed exp1/exp2/exp3 with timetables of
   %                     snow_liquid_water_storage_m and bottom_outflow_mps.
   %       metadata      solver_kind, kernel provenance, options used.
   %
   %  See also: icemodel.column.infiltration,
   %    icemodel.verification.colbeck.analyticalSolution,
   %    icemodel.verification.colbeck.compareSolutions

   arguments
      case_manifest %#ok<INUSA>
      kwargs.kind (1, 1) string {mustBeMember(kwargs.kind, ...
         ["numerical", "analytical"])} = "numerical"
      kwargs.experiment_names (1, :) string = ["exp1", "exp2", "exp3"]
   end

   def = icemodel.verification.colbeck.caseDefinition();

   experiments = struct();
   solver_meta = struct();
   for i = 1:numel(kwargs.experiment_names)
      name = kwargs.experiment_names(i);
      switch kwargs.kind
         case "numerical"
            [tt, meta] = run_numerical(name, def);
         case "analytical"
            [tt, meta] = run_analytical(name, def);
      end
      experiments.(char(name)) = tt;
      solver_meta.(char(name)) = meta;
   end

   switch kwargs.kind
      case "numerical"
         metadata = struct( ...
            'solver_kind', "numerical_icemodel", ...
            'kernel',      "icemodel.column.infiltration", ...
            'note',        ['IceModel Colbeck numerical candidate. ' ...
            'Uses icemodel.column.infiltration for the ' ...
            'full liquid-mass + cold-content + ' ...
            'conduction step.'], ...
            'per_experiment', solver_meta);
      case "analytical"
         metadata = struct( ...
            'solver_kind', "analytical_icemodel", ...
            'kernel',      "icemodel.verification.colbeck.analyticalSolution", ...
            'note',        ['IceModel Colbeck analytical candidate. ' ...
            'Closed-form Clark 2017 wetting-front / ' ...
            'kinematic-wave solution.'], ...
            'per_experiment', solver_meta);
   end

   candidate = struct( ...
      'format',      'experiment_bundle', ...
      'experiments', experiments, ...
      'metadata',    metadata);
end

% =====================================================================
% Numerical branch: time loop calling icemodel.column.infiltration.
% =====================================================================
function [tt, meta] = run_numerical(experiment_name, def)

   % Resolve the named experiment from the case definition.
   names = string({def.experiments.name});
   idx = find(names == experiment_name, 1);
   if isempty(idx)
      error('icemodel:colbeck:unknownExperiment', ...
         'unknown Colbeck experiment: %s', experiment_name);
   end
   exp = def.experiments(idx);

   % Grid, time, surface forcing, and k_sat dispatch values.
   %
   % The Colbeck verification pins kappa to the SUMMA Laugh-Tests parameter
   % trial value (matches Shimizu evaluated at initial state), so we use method
   % = "darcy" with the case-specific permeability.
   N           = def.grid.n_layers;
   dz          = def.grid.dz_m;
   dt          = def.time.dt_s;
   n_steps     = def.time.n_steps;
   q_top_full  = def.boundary.q_top_m_per_s;
   rain_window = def.time.rain_window_s;
   permeability = exp.permeability_m2;

   % Initial column state (uniform across layers per Clark 2017 Table 1).
   T     = exp.T_K   * ones(N, 1);
   f_ice = exp.f_ice * ones(N, 1);
   f_liq = exp.f_liq * ones(N, 1);

   % Per-output-step diagnostics.
   storage       = zeros(n_steps, 1);
   q_bot         = zeros(n_steps, 1);
   inflow_total  = 0;
   outflow_total = 0;

   % Time loop: one infiltration call per output step. The kernel handles its
   % own CFL substepping internally, so this loop body is deliberately thin —
   % all physics lives in icemodel.column.infiltration.
   for step = 1:n_steps

      % Top boundary follows the Colbeck rain window. q_top is constant during
      % rain (steps 1..rain_window/dt) and zero afterwards.
      t_end = step * dt;
      q_top = q_top_full * (t_end <= rain_window + eps);

      % Advance the column one main step. Returns updated state plus a diag
      % struct with inflow / outflow totals for this step.
      [f_liq, f_ice, T, diag] = icemodel.column.infiltration( ...
         f_liq, f_ice, T, dz, dt, q_top, ...
         k_sat_method="darcy", permeability=permeability);

      % Step diagnostics: column-integrated liquid storage [m] and the
      % step-averaged bottom outflow [m s-1].
      storage(step) = sum(f_liq) * dz;
      q_bot(step)   = diag.outflow_total / dt;
      inflow_total  = inflow_total  + diag.inflow_total;
      outflow_total = outflow_total + diag.outflow_total;
   end

   % Package the per-step diagnostics as a timetable matching the staged target
   % schema.
   time_seconds  = (1:n_steps).' * dt;
   time_datetime = datetime(1990, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') ...
      + seconds(time_seconds);
   tt = timetable(time_datetime, storage, q_bot, ...
      'VariableNames', ...
      {'snow_liquid_water_storage_m', 'bottom_outflow_mps'});
   tt.Properties.DimensionNames{1} = 'Time';

   % Per-experiment provenance metadata for the candidate bundle.
   meta = struct( ...
      'inflow_total_m',  inflow_total, ...
      'outflow_total_m', outflow_total, ...
      'final_T_min',     min(T), ...
      'final_T_max',     max(T), ...
      'k_sat_method',    "darcy", ...
      'permeability_m2', permeability);
end

% =====================================================================
% Analytical branch: wrap analyticalSolution into the bundle schema.
% =====================================================================
function [tt, meta] = run_analytical(experiment_name, def)
   sol = icemodel.verification.colbeck.analyticalSolution( ...
      experiment_name, def);
   tt = timetable(sol.time_datetime, ...
      sol.snow_liquid_water_storage_m, sol.bottom_outflow_mps, ...
      'VariableNames', {'snow_liquid_water_storage_m', 'bottom_outflow_mps'});
   tt.Properties.DimensionNames{1} = 'Time';
   meta = sol.parameters;
end
