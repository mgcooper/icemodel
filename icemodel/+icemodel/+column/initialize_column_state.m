function [ice1, ice2, T_sfc, T_ice, f_ice, f_liq, Sc, Sp, r_eff, k_eff, fn, ...
      dz, delz, z_nodes, f_liq_res] = initialize_column_state(opts, tair, r_eff)
   %INITIALIZE_COLUMN_STATE Initialize the 1-d ice column state.
   %
   %  [ice1, ice2, T_sfc, T_ice, f_ice, f_liq, Sc, Sp, r_eff, k_eff, fn, ...
   %   dz, delz, z_nodes, f_liq_res] = ...
   %     icemodel.column.initialize_column_state(opts, tair, r_eff)
   %
   %  Returns the smallest set of column-state variables that the model kernel
   %  needs. This function reads physical constants (TL, TH, Lf, …) from
   %  icemodel.physicalConstant. It reads thermodynamic parameters
   %  (f_ell_min, f_ell_max, ro_iwe, …) from icemodel.parameterLookup.
   %
   %  opts.use_ro_glc changes only the densities used to construct initial
   %  phase fractions. Solvers use persistent physical constants. See setopts.
   %
   % See also: icemodel, skinmodel
   %
   %#codegen

   debug = false;

   % LOAD PHYSICAL CONSTANTS AND PARAMETERS
   [cp_ice, ro_ice, ro_liq, k_ice, Tf] ...
      = icemodel.physicalConstant('cp_ice', 'ro_ice', 'ro_liq', 'k_ice', 'Tf');

   fcp = icemodel.parameterLookup('fcp');

   % GENERATE A THERMAL MESH
   [dz, delz, z_nodes, ~, fn] = ...
      icemodel.column.control_volume_mesh(opts.z0_thermal, opts.dz_thermal);

   % NUMBER OF TIMESTEPS TO INITIALIZE OUTPUTS
   maxiter = numel(tair) / opts.numyears;

   % NUMBER OF NODES IN THE THERMAL MESH
   JJ = numel(dz);

   % DERIVED INTRINSIC DENSITIES
   ro_glc = (ro_ice + ro_liq) / 2; % [kg/m3]
   if opts.use_ro_glc == true
      ro_ice = ro_glc;
      ro_liq = ro_glc;
   end

   % Lower melt-zone temperature bound (needed for initial T_ice profile only)
   TL = icemodel.parameterLookup('TL');

   % INITIALIZE CORE STATE VARIABLES
   if opts.use_restart
      restart = icemodel.loadRestartState(opts);
      validateRestartOpts(opts, restart);
      T_ice = restart.T_ice;
      f_ice = restart.f_ice;
      f_liq = restart.f_liq;
      T_sfc = restart.T_sfc;
      r_eff = restart.r_eff;
      validateRestartState(T_ice, f_ice, f_liq, T_sfc, r_eff, JJ);
   else
      % Initialize with a physically scaled exponential profile anchored to the
      % top node so T_ice(1) = T_ref exactly. Use a slight cold offset from air
      % to avoid a zero-gradient startup.
      T_bot = min(TL - 1, Tf + opts.T_ice_init);
      T_ref = min(TL - 1, min(tair(1), Tf) - 1.0);
      omega_yr = 2.0 * pi / (365.0 * 86400.0);
      kappa_ice = k_ice / (ro_ice * cp_ice);
      z_scale = sqrt(2.0 * kappa_ice / omega_yr);
      T_ice = T_bot + (T_ref - T_bot) ...
         .* exp(-(z_nodes - z_nodes(1)) ./ z_scale);

      % INITIALIZE LIQUID/ICE WATER FRACTION (f) AND BULK DENSITIES (g)
      T_dep = Tf - T_ice;                          % [K]
      f_ell = 1 ./ (1 + (fcp * T_dep) .^ 2);       % [-], liquid mass fraction
      g_ice = opts.ro_ice_init;                    % [kg m-3]
      g_liq = g_ice .* (f_ell ./ (1 - f_ell));     % [kg m-3]
      f_liq = g_liq ./ ro_liq;                     % [1]
      f_ice = g_ice ./ ro_ice .* ones(JJ, 1);      % [1]
      T_sfc = (min(tair(1), Tf) + T_ice(1)) / 2;   % [K]
      r_eff = r_eff / 1000 * ones(JJ, 1);          % [m] mm->m for grain growth
   end

   % INITIAL SOLVER CONDUCTIVITY
   % The column solvers add vapor transport on faces. Initialize the
   % checkpoint with the same vapor-free node conductivity that subsequent
   % accepted solves return.
   k_eff = icemodel.column.bulk_thermal_conductivity(T_ice, f_ice, f_liq, 0);

   % SOURCE TERM LINEARIZATION VECTORS
   Sc = zeros(JJ, 1);
   Sp = zeros(JJ, 1);

   % RESIDUAL LIQUID-WATER FRACTION (set to the snow value; each substep
   % overrides it)
   f_liq_res = opts.f_res_pore_snow;

   % INITIALIZE THE OUTPUT STRUCTURES
   for n = 1:numel(opts.vars1) % ice1 = 1-d data
      ice1.(opts.vars1{n}) = nan(maxiter, 1);
   end
   for n = 1:numel(opts.vars2) % ice2 = 2-d data
      if strcmp(opts.vars2{n}, 'lcflag')
         ice2.(opts.vars2{n}) = false(JJ, maxiter);
      else
         ice2.(opts.vars2{n}) = nan(JJ, maxiter);
      end
   end

   function plot_T_init(T_ice, Z, Ta, T_sfc)
      figure; hold on
      plot(T_ice, Z)
      scatter(T_sfc, 0, 'filled')
      scatter(Ta, 0, 'filled')
      set(gca, 'YDir', 'reverse')
      legend('T_ice', 'T_sfc', 'Ta')
   end
   if debug == true
      plot_T_init(T_ice, z_nodes, tair(1), T_sfc);
   end
end

function validateRestartState(T_ice, f_ice, f_liq, T_sfc, r_eff, JJ)

   state_vars = {T_ice, f_ice, f_liq, r_eff};
   for n = 1:numel(state_vars)
      x = state_vars{n};
      if ~isvector(x) || numel(x) ~= JJ
         error('restart state size does not match the current thermal mesh')
      end
   end

   if ~isscalar(T_sfc) || ~isfinite(T_sfc)
      error('restart state must provide a finite scalar T_sfc')
   end
end

function validateRestartOpts(opts, restart)

   if ~isfield(restart, 'opts') || ~isstruct(restart.opts)
      return
   end

   saved = restart.opts;
   requireMatchingField(opts, saved, 'smbmodel');
   requireMatchingField(opts, saved, 'dz_thermal');
   requireMatchingField(opts, saved, 'z0_thermal');
   requireMatchingField(opts, saved, 'use_ro_glc');
end

function requireMatchingField(current, saved, name)

   if ~isfield(current, name) || ~isfield(saved, name)
      return
   end

   x = current.(name);
   y = saved.(name);
   if ischar(x) || isstring(x)
      ok = strcmpi(string(x), string(y));
   else
      ok = isequaln(x, y);
   end

   if ~ok
      error('restart opts mismatch for "%s"', name)
   end
end
