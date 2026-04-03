function [T_sfc, ok] = SEBSOLVE(T_sfc, tair, Qsi, Qli, albedo, wspd, ppt, tppt, psfc, ...
      De, ea_atm, chi, roL, br_coefs, liqflag, T, k_eff, dz, ro_sfc, snow_depth, ...
      solver, debug, opts)
   %SEBSOLVE solve the surface energy balance for the skin temperature
   %
   % Solver options:
   %  0 = derivative-free, slower but more robust, used as a fall back.
   %  1 = newton-rhapson, fast, requires analytic derivative, used as a default.
   %  2 = complex-step, fast and does not require an analytical derivative.
   % -1 = no outer iterations, use this if phase change is not represented
   %      explicitly in the model.
   % >2 = experimental.
   %
   % Important programming notes:
   %  - The outer iterations control the convergence of the Ts calculation wrt
   %  the conduction term. Thus old is initialized to Ts outside the outer loop.
   %
   %  - In SFCTEMP or complexstep or any derivative-based method, old must be
   %  initialized to Ta to avoid divergence during spinup, keeping in mind that
   %  the conduction passed into SFCTEMP is computed with old = Ts.
   %
   %  - For a "skinmodel", Ts never exceeds Tf when passed into functions, but
   %  within the iterations of SFCTEMP and when it comes out of SFCTEMP it can
   %  exceed Tf.
   %
   %#codegen

   tol = 1e-3;
   seb_tol = 1.0;
   maxiter = 100;
   iterflag = true;
   if solver < 0
      maxiter = 1;
      solver = -solver;
      iterflag = false;
   end

   ok_cpl = false;
   T_sfc_old = T_sfc;
   switch solver

      case 1 % Newton-Rhapson - analytical derivative

         for iter = 1:maxiter
            old = T_sfc;
            [T_sfc, ok] = SFCTEMP(tair, Qsi, Qli, albedo, wspd, ppt, tppt, psfc, ...
               De, ea_atm, chi, roL, br_coefs, liqflag, CONDUCT(k_eff, T, dz, old), ...
               k_eff, T, dz, iterflag);

            if ~ok
               break
            end
            if abs(old - T_sfc) < tol && abs(fSEB(T_sfc)) < seb_tol
               ok_cpl = true;
               break
            end
         end

      case 2 % Complex-step - numerical derivative

         for iter = 1:maxiter
            old = T_sfc;
            [T_sfc, ok] = complexstep(@fSEB, old);
            if ~ok
               break
            end
            if abs(old - T_sfc) < tol && abs(fSEB(T_sfc)) < seb_tol
               ok_cpl = true;
               break
            end
         end

      otherwise % Derivative-free
         solver = 0;
   end

   % Derivative-free solver
   if solver == 0 || ~ok_cpl || abs(T_sfc - tair) > 20

      T_sfc = T_sfc_old;

      for iter = 1:maxiter
         old = T_sfc;
         [T_sfc, ~, ok] = fsearchzero(@fSEB, old, tair-50, tair+50, tair, tol);
         if ~ok
            break
         end
         if abs(old - T_sfc) < tol && abs(fSEB(T_sfc)) < seb_tol
            ok_cpl = true;
            break
         end
      end
   end

   % Hitting max coupling iterations without ok_cpl is a substep fail
   ok = ok && ok_cpl;

   if ~ok && debug
      dumpSebSolveFailure(solver, iter, T_sfc_old, T_sfc, tair, Qsi, Qli, albedo, ...
         wspd, ppt, tppt, psfc, De, ea_atm, chi, roL, br_coefs, liqflag, k_eff, T, ...
         dz, ro_sfc, snow_depth, ok_cpl, opts);
   end

   % Note: nested function captures updated Qe, Qh, and Qc on each iteration.
   function f = fSEB(T_sfc_local)
      f = icemodel.surface.surface_energy_balance_residual(T_sfc_local, tair, Qsi, ...
         Qli, albedo, wspd, ppt, tppt, psfc, De, ea_atm, chi, roL, br_coefs, ...
         CONDUCT(k_eff, T, dz, T_sfc_local), liqflag, ro_sfc, snow_depth, opts);
   end
end

function [x, ok, iter] = complexstep(f, x0)
   %COMPLEXSTEP Find root of nonlinear function using the Newton Rhapson method.
   %
   % Note: this uses the complex-step numerical derivative.
   %
   % The Newton iterate itself must remain on the real axis. The
   % underlying SEB residual contains regime switches (for example the
   % stable/unstable bulk-Richardson branch selection), so the complex-step
   % perturbation is only used to obtain dfdx at a real state. The
   % unperturbed residual f(old) is required to be real-valued on the real
   % axis; if not, fail loudly because that indicates a contract bug in the
   % residual implementation rather than a valid Newton state.

   persistent h dh tol maxiter imag_factor_tol
   if isempty(tol)
      h = 1e-10;
      dh = 1i * h;
      tol = 1e-3;
      maxiter = 100;
      imag_factor_tol = 100;
   end

   ok = false;
   old = real(x0);
   for iter = 1:maxiter

      f_old = f(old);
      imag_tol = imag_factor_tol * eps(max(1.0, abs(real(f_old))));
      if abs(imag(f_old)) > imag_tol
         error('icemodel:ComplexStepResidualNotReal', ...
            ['complexstep residual must stay real on the real axis; ', ...
            'got imag(f(x)) = %.3e at x = %.15g.'], imag(f_old), old);
      end
      dfdx = imag(f(old + dh)) / h;
      x = real(old - real(f_old) / dfdx);

      if abs(x - old) < tol
         ok = true;
         return

      elseif ~isfinite(x)
         x = x0;
         return
      end
      old = x;
   end
   x = x0; % ok = false
end

function dumpSebSolveFailure(solver, iter, T_sfc_old, T_sfc, tair, Qsi, Qli, ...
      albedo, wspd, ppt, tppt, psfc, De, ea_atm, chi, roL, br_coefs, liqflag, ...
      k_eff, T, dz, ro_sfc, snow_depth, ok_cpl, opts)
   %DUMPSEBSOLVEFAILURE Save SEB root-find failure diagnostics on demand.

   debug_file = getenv('ICEMODEL_DEBUG_SEBSOLVE_FILE');
   if isempty(debug_file)
      return
   end

   debug_state = struct();
   debug_state.timestamp_utc = datetime('now', 'TimeZone', 'UTC');
   debug_state.solver = solver;
   debug_state.iter = iter;
   debug_state.Ts_old = T_sfc_old;
   debug_state.Ts = T_sfc;
   debug_state.tair = tair;
   debug_state.Qsi = Qsi;
   debug_state.Qli = Qli;
   debug_state.albedo = albedo;
   debug_state.wspd = wspd;
   debug_state.ppt = ppt;
   debug_state.tppt = tppt;
   debug_state.psfc = psfc;
   debug_state.De = De;
   debug_state.ea_atm = ea_atm;
   debug_state.chi = chi;
   debug_state.roL = roL;
   debug_state.br_coefs = br_coefs;
   debug_state.liqflag = liqflag;
   debug_state.k_eff = k_eff;
   debug_state.T = T;
   debug_state.dz = dz;
   debug_state.ro_sfc = ro_sfc;
   debug_state.snow_depth = snow_depth;
   debug_state.turbulent_flux_scheme = char(opts.turbulent_flux_scheme);
   debug_state.ok_cpl = ok_cpl;

   save(debug_file, 'debug_state');

   icemodel.surface.dump_turbulent_heat_flux_debug_state( ...
      "sebsolve_nonconvergence", T_sfc_old, T_sfc, tair, Qsi, Qli, albedo, wspd, ...
      ppt, tppt, psfc, De, ea_atm, chi, roL, br_coefs, liqflag, T, k_eff, dz, ...
      ro_sfc, snow_depth, opts);
end
