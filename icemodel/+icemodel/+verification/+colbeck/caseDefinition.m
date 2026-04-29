function def = caseDefinition()
   %CASEDEFINITION Return the canonical Colbeck 1976 verification case.
   %
   %  def = icemodel.verification.colbeck.caseDefinition()
   %
   % Outputs
   %  def   Struct with the canonical Colbeck 1976 infiltration experiment
   %        parameters, grid, time window, surface boundary, and per-experiment
   %        initial state needed by analytical solutions, the IceModel numerical
   %        candidate, and verification target packaging.
   %
   % Role
   %  Single source of truth for Colbeck 1976 case constants. Analytical and
   %  numerical solvers must consume this definition rather than duplicate
   %  values. The struct intentionally keeps fields flat and named with units so
   %  downstream code does not need to translate values.
   %
   % References
   %  Clark, M. P., et al. (2021): The numerical implementation of land
   %  models: Problem formulation and laugh tests. JHM, Eq. 60-64 and Table 5.
   %  Colbeck, S. C. (1976): An analysis of water flow in dry snow.
   %  Laugh-Tests bundle: KyleKlenk/Laugh-Tests, validation_data/m2_mac_Sept23
   %  for SUMMA reference state and parameter trial files.

   % Three experiments share the same grid, timestep, and surface forcing. They
   % differ in initial T, f_ice, f_liq, and snow grain size. Values are the
   % SUMMA initial conditions and parameter trials used by Laugh-Tests; see
   % provenance below.
   %
   % permeability_m2 is the intrinsic permeability kappa [m^2] used by the
   % "darcy" k_sat method: k_sat = kappa * ro_liq * g / mu_water_ref. This is
   % SUMMA's effective `k_snow` from the Laugh-Tests Colbeck parameter trial
   % file, expressed as a permeability so the units are physically transparent.
   % The values match Shimizu 1970 evaluated at the initial state for each
   % experiment:
   %
   %   kappa = shi_k0 * grainsz^2 * exp(-shi_a_wat * f_wat_initial)
   %
   % with shi_k0 = 0.077, shi_a_wat = 7.8 from parameterLookup. For all three
   % experiments f_wat_initial = 0.300 (exp1: 0.0510 + 0.917 * 0.272; exp2/exp3:
   % 0.917 * 0.327), so the only state difference between experiments is
   % grainsz: 2 mm for exp1/exp2 and 0.2 mm for exp3. The 100x grainsz^2 factor
   % explains the 100x permeability drop for exp3. SUMMA freezes kappa at the
   % initial state and reuses it for the run, hence the constant value here.
   experiments = struct( ...
      'name',          {'exp1';        'exp2';        'exp3'}, ...
      'description',   {'ripe';        'refrozen';    'fresh'}, ...
      'T_K',           {273.116;       268.16;        268.16}, ...
      'f_ice',         {0.271547;      0.327154;      0.327154}, ...
      'f_liq',         {0.0509917;     0.0;           0.0}, ...
      'grainsz_m',     {2.0e-3;        2.0e-3;        2.0e-4}, ...
      'permeability_m2', {2.967e-8;    2.967e-8;      2.967e-10});
   % mgc note: fcp = (sqrt(1/0.0509917 - 1)) / (273.16 - 273.116) = 98.0466

   def.experiments = experiments;

   % Uniform vertical grid: 100 layers x 0.01 m = 1.0 m snow column.
   def.grid = struct( ...
      'n_layers',      100, ...
      'dz_m',          0.01, ...
      'total_depth_m', 1.00);

   % Time window: 60 s steps for 600 steps = 10 hours total simulation time. The
   % surface liquid flux is active for the first 3 hours (10800 s) then shuts
   % off, matching Colbeck 1976 / Clark et al. 2021.
   def.time = struct( ...
      'dt_s',           60, ...
      'n_steps',        600, ...
      't_end_s',        36000, ...
      'rain_window_s',  10800);

   % Surface and bottom boundary fluxes. The SUMMA forcing supplies pptrate as
   % 0.01 kg m-2 s-1 (mass flux of liquid water at the surface). Convert to a
   % volumetric water flux using ro_liq so analytical and numerical solvers
   % consume one boundary value with consistent units.
   ro_liq = icemodel.physicalConstant('ro_liq');
   q_top_kg_per_m2_per_s = 0.01;
   def.boundary = struct( ...
      'q_top_kg_per_m2_per_s', q_top_kg_per_m2_per_s, ...
      'q_top_m_per_s',         q_top_kg_per_m2_per_s / ro_liq, ...
      'q_bottom',              "free");

   % Snow water-flow case parameters. f_res_pore matches SUMMA's Fcapil
   % parameter for the Colbeck verification target. The constitutive-law
   % exponent (m_exp) and other empirical coefficients are loaded from
   % icemodel.parameterLookup at the kernel boundary; this struct records only
   % case-specific values that vary with the experiment.
   def.params = struct( ...
      'f_res_pore', 0.07);

   % Provenance: pin every value to a paper or upstream file so developers
   % can audit numbers without re-deriving them.
   def.provenance = struct( ...
      'clark_2021_eqs',    "Eq. 60-64", ...
      'clark_2021_table5', "Table 5 / pp. 11-13", ...
      'colbeck_1976',      "Colbeck (1976), An analysis of water flow in dry snow", ...
      'summa_initcond',    fullfile("test_cases", "settings", ...
      "syntheticTestCases", "colbeck1976", "common", ...
      "summa_zInitialCond_colbeck1976-exp{1,2,3}.nc"), ...
      'summa_paramtrial',  fullfile("test_cases", "settings", ...
      "syntheticTestCases", "colbeck1976", "common", ...
      "summa_zParamTrial_colbeck1976-exp{1,2,3}.nc"), ...
      'summa_forcing',     fullfile("test_cases", "input_data", ...
      "colbeck1976", "colbeck1976_forcing.nc"), ...
      'laughtests_url',    "https://github.com/KyleKlenk/Laugh-Tests");
end
