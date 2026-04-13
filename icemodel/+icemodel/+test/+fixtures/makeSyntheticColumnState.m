function state = makeSyntheticColumnState(workspace, smbmodel, kwargs)
   %MAKESYNTHETICCOLUMNSTATE Build a resolved synthetic column kernel state.
   %
   % state = icemodel.test.fixtures.makeSyntheticColumnState(...
   %    workspace, "icemodel")
   %
   % The returned struct contains a benign initialized column, resolved opts,
   % forcing data for one met step, and the physical constants needed by the
   % core solver kernels.

   arguments
      workspace struct
      smbmodel (1, :) char {mustBeMember(smbmodel, {'icemodel', 'skinmodel'})}
      kwargs.simyears = []
      kwargs.solver (1, 1) double = NaN
      kwargs.seb_solver (1, 1) double = NaN
      kwargs.turbulent_flux_scheme (1, :) char = ''
      kwargs.z0_bulk (1, 1) double = NaN
      kwargs.z0_ice (1, 1) double = NaN
      kwargs.metstep (1, 1) double {mustBeInteger, mustBePositive} = 1
      kwargs.testname (1, :) char = 'column_state'
      kwargs.include_spectral (1, 1) logical = false
   end

   simyears = kwargs.simyears;
   if isempty(simyears)
      simyears = workspace.simyears(1);
   end

   % Resolve one benign OPTS struct and load the matching forcing slice.
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, smbmodel, simyears, ...
      solver=kwargs.solver, seb_solver=kwargs.seb_solver, ...
      turbulent_flux_scheme=kwargs.turbulent_flux_scheme, ...
      z0_bulk=kwargs.z0_bulk, ...
      z0_ice=kwargs.z0_ice, testname=kwargs.testname, ...
      output_profile='standard');

   met = icemodel.loadmet(opts);

   % Seed the thermal grain-radius state only for icemodel.
   r_eff = 0;
   if strcmp(smbmodel, 'icemodel')
      [~, ~, ~, ~, ~, ~, ~, ~, r_eff] = ...
         icemodel.radiation.initialize_spectral_model(opts);
   end

   % Initialize the column state exactly as the model kernel would.
   [ice1, ice2, Ts, T, f_ice, f_liq, r_eff, k_eff, fn, dz, delz, ...
      z_nodes, ro_air_Lv, liqflag, JJ, Sc, Sp, Fc, Fp, TL, TH, f_ell_min, ...
      f_ell_max, f_ice_min, f_liq_res, ro_iwe, ro_wie] = ...
      icemodel.column.initialize_column_state(opts, met.tair, r_eff);

   metstep = kwargs.metstep;

   % Extract one forcing step and the corresponding transfer coefficients.
   [tair, swd, lwd, albedo, wspd, psfc, ea, De, forcing_snow_depth] ...
      = icemodel.timestepping.getforcings(met, metstep, liqflag, opts);

   [~, scoef] = ...
      icemodel.surface.turbulence.bulk_richardson.exchange_coefficients( ...
      wspd, opts.z0_bulk, opts.z_tair, opts.z_wind);

   % Carry precipitation fields when the synthetic met fixture defines them.
   if ismember('ppt', met.Properties.VariableNames)
      ppt = met.ppt(metstep);
   else
      ppt = 0.0;
   end
   if ismember('tppt', met.Properties.VariableNames)
      tppt = met.tppt(metstep);
   else
      tppt = tair;
   end

   [cv_air, cv_liq, cv_ice, SB, k_liq, ro_ice, ro_liq, ro_air, ...
      Ls, Lf, roLf, Rv, Tf, epsilon, roLs, roLv] = ...
      icemodel.physicalConstant('cv_air', 'cv_liq', 'cv_ice', ...
      'SB', 'k_liq', 'ro_ice', 'ro_liq', 'ro_air', 'Ls', 'Lf', 'roLf', ...
      'Rv', 'Tf', 'epsilon', 'roLs', 'roLv');
   [emiss, fcp] = icemodel.parameterLookup('emiss', 'fcp');

   % Package the resolved state so kernel tests can reuse it directly.
   state = struct();
   state.workspace = workspace;
   state.opts = opts;
   state.met = met;
   state.ice1 = ice1;
   state.ice2 = ice2;
   state.metstep = metstep;

   state.T = T;
   state.f_ice = f_ice;
   state.f_liq = f_liq;
   state.r_eff = r_eff;
   state.f_wat = icemodel.column.water_fraction(f_ice, f_liq);
   state.k_eff = k_eff;
   state.fn = fn;
   state.dz = dz;
   state.delz = delz;
   state.z_nodes = z_nodes;
   state.Ts = Ts;
   state.ro_sfc = icemodel.surface.surface_bulk_density(f_ice(1), f_liq(1));
   state.snow_depth = icemodel.surface.resolve_forcing_snow_depth( ...
      forcing_snow_depth, 1, opts.use_forcing_snow_depth_for_thf);
   state.JJ = JJ;
   state.Sc = Sc;
   state.Sp = Sp;
   state.Fc = Fc;
   state.Fp = Fp;
   state.TL = TL;
   state.TH = TH;
   state.f_ell_min = f_ell_min;
   state.f_ell_max = f_ell_max;
   state.f_ice_min = f_ice_min;
   state.f_liq_res = f_liq_res;
   state.ro_iwe = ro_iwe;
   state.ro_wie = ro_wie;
   state.liqflag = liqflag;
   state.ro_air_Lv = ro_air_Lv;

   state.tair = tair;
   state.swd = swd;
   state.lwd = lwd;
   state.albedo = albedo;
   state.wspd = wspd;
   state.ppt = ppt;
   state.tppt = tppt;
   state.psfc = psfc;
   state.De = De;
   state.ea_atm = ea;
   state.br_coefs = scoef;
   state.chi = 1.0;

   state.cv_air = cv_air;
   state.cv_liq = cv_liq;
   state.cv_ice = cv_ice;
   state.emiss = emiss;
   state.SB = SB;
   state.k_liq = k_liq;
   state.ro_ice = ro_ice;
   state.ro_liq = ro_liq;
   state.ro_air = ro_air;
   state.Ls = Ls;
   state.Lf = Lf;
   state.roLf = roLf;
   state.Rv = Rv;
   state.Tf = Tf;
   state.epsilon = epsilon;
   state.roLs = roLs;
   state.roLv = roLv;
   state.fcp = fcp;

   state.tol = opts.tol;
   state.maxiter = opts.maxiter;
   state.alpha = opts.alpha;
   state.use_aitken = opts.use_aitken;
   state.jumpmax = opts.jumpmax;
   state.seb_solver = opts.seb_solver;
   state.cpl_maxiter = opts.cpl_maxiter;
   state.cpl_Ts_tol = opts.cpl_Ts_tol;
   state.cpl_seb_tol = opts.cpl_seb_tol;
   state.cpl_alpha = opts.cpl_alpha;
   state.cpl_aitken = opts.cpl_aitken;
   state.cpl_jumpmax = opts.cpl_jumpmax;

   % Attach the spectral grid only for tests that exercise that path.
   if kwargs.include_spectral
      [I0, dz_spect, z_nodes_spect, z_edges_spect, tau_N, tau_S, ...
         solar_dwavel, ~, ~, qext, g, coalbedo, kabs, kice, wavel, radii] ...
         = icemodel.radiation.initialize_spectral_model(opts);

      % Recover k_ext for tests that need the initialized coefficients.
      [~, ~, ~, k_ext] = icemodel.radiation.update_extinction_coefficients( ...
         qext, g, coalbedo, kabs, kice, wavel, radii, opts.i_grainradius, ...
         z_edges_spect, dz_spect, solar_dwavel, ro_ice, false);

      state.I0 = I0;
      state.dz_spect = dz_spect;
      state.z_nodes_spect = z_nodes_spect;
      state.z_edges_spect = z_edges_spect;
      state.tau_N = tau_N;
      state.tau_S = tau_S;
      state.solar_dwavel = solar_dwavel;
      state.qext = qext;
      state.g = g;
      state.coalbedo = coalbedo;
      state.radii = radii;
      state.wavel = wavel;
      state.k_ext = k_ext;
      state.kice = kice;
      state.kabs = kabs;
   end
end
