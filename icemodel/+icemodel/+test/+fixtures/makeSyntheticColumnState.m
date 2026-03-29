function state = makeSyntheticColumnState(workspace, smbmodel, kwargs)
   %MAKESYNTHETICCOLUMNSTATE Build a resolved synthetic column kernel state.
   %
   %  state = icemodel.test.fixtures.makeSyntheticColumnState(workspace, "icemodel")
   %
   % The returned struct contains a benign initialized column, resolved opts,
   % forcing data for one met step, and the physical constants needed by the
   % core solver kernels.

   arguments
      workspace struct
      smbmodel (1, :) char {mustBeMember(smbmodel, {'icemodel', 'skinmodel'})}
      kwargs.simyears = []
      kwargs.solver (1, 1) double = NaN
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
      solver=kwargs.solver, testname=kwargs.testname, output_profile='standard');
   met = icemodel.loadmet(opts);

   % Initialize the column state exactly as the model kernel would.
   [ice1, ice2, T, f_ice, f_liq, k_eff, fn, dz, delz, z_nodes, roL, ...
      liqflag, Ts, JJ, Sc, Sp, Fc, Fp, TL, TH, f_ell_min, f_ell_max, ...
      f_ice_min, f_liq_res, ro_iwe, ro_wie] = ICEINIT(opts, met.tair);

   metstep = kwargs.metstep;
   % Extract one forcing step and the corresponding transfer coefficients.
   [tair, swd, lwd, albedo, wspd, psfc, De, ea] = LOADMETDATA(met, metstep, ...
      liqflag);
   [~, scoef] = WINDCOEF(wspd, opts.z_0, opts.z_tair, opts.z_wind);

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
   state.f_wat = f_liq + f_ice * ro_ice / ro_liq;
   state.k_eff = k_eff;
   state.fn = fn;
   state.dz = dz;
   state.delz = delz;
   state.z_nodes = z_nodes;
   state.Ts = Ts;
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
   state.roL = roL;

   state.tair = tair;
   state.swd = swd;
   state.lwd = lwd;
   state.albedo = albedo;
   state.wspd = wspd;
   state.ppt = ppt;
   state.tppt = tppt;
   state.psfc = psfc;
   state.De = De;
   state.ea = ea;
   state.scoef = scoef;
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
      [I0, dz_spect, z_nodes_spect, z_edges_spect, tau_N, ...
         tau_S, solar_dwavel, ~, ~, qext, g, coalbedo, kabs, kice, ...
         wavel, radii] = EXTCOEFSINIT(opts, ro_ice);

      % Recover k_ext for tests that need the initialized coefficients.
      [~, ~, ~, k_ext] = UPDATEEXTCOEFS(qext, g, coalbedo, kabs, kice, ...
         wavel, radii, opts.i_grainradius, z_edges_spect, dz_spect, ...
         solar_dwavel, ro_ice, false);

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
