function [T_sfc, ok] = solve_surface_energy_balance(T_sfc, tair, Qsi, Qli, ...
      albedo, wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, ro_air_Lv, liqflag, ...
      chi, T_ice, k_eff, dz, ro_sfc, snow_depth, opts)
   %SOLVE_SURFACE_ENERGY_BALANCE Solve the nonlinear surface energy balance.
   %
   %  [T_sfc, ok] = icemodel.surface.solve_surface_energy_balance(...)
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
   %  - The outer iterations provide a shared convergence / fallback wrapper
   %  across the analytical Newton, complex-step, and derivative-free paths.
   %  old is initialized to T_sfc outside the outer loop and is passed as
   %  the initial guess to the selected inner root finder.
   %
   %  - solve_surface_temperature starts its Newton-Raphson from the outer
   %  iterate (old = T_sfc on each call), and evaluates Qc and its derivative
   %  at each inner step.
   %
   %  - For a "skinmodel", Ts never exceeds Tf when passed into functions, but
   %  within the iterations of solve_surface_temperature and when it comes out
   %  of solve_surface_temperature it can exceed Tf.
   %
   %#codegen

   debug = opts.debug;
   seb_solver = opts.seb_solver;

   tol = 1e-3;
   seb_tol = 1.0;
   maxiter = 100;

   % Diagnostic option: seb_solver < 0 limits to one outer iteration.
   % Useful for testing with opts.seb_solver = -1, -2, etc.
   if seb_solver < 0
      maxiter = 1;
      seb_solver = -seb_solver;
   end

   ok_cpl = false;
   T_sfc_old = T_sfc;
   switch seb_solver

      case 1 % Newton-Raphson - analytical derivative

         for iter = 1:maxiter
            old = T_sfc;
            [T_sfc, ok] = icemodel.surface.solve_surface_temperature( ...
               old, tair, Qsi, Qli, albedo, wspd, ppt, tppt, psfc, De, ...
               ea_atm, br_coefs, ro_air_Lv, liqflag, chi, k_eff, T_ice, dz);

            if ~ok
               break
            end
            if abs(old - T_sfc) < tol && abs(surface_residual(T_sfc)) < seb_tol
               ok_cpl = true;
               break
            end
         end

      case 2 % Complex-step - numerical derivative

         for iter = 1:maxiter
            old = T_sfc;
            [T_sfc, ok] = icemodel.numerics.complexstep(@surface_residual, old);
            if ~ok
               break
            end
            if abs(old - T_sfc) < tol && abs(surface_residual(T_sfc)) < seb_tol
               ok_cpl = true;
               break
            end
         end

      otherwise % Derivative-free
         seb_solver = 0;
   end

   % Derivative-free solver
   if seb_solver == 0 || ~ok_cpl || abs(T_sfc - tair) > 20

      T_sfc = T_sfc_old;

      for iter = 1:maxiter
         old = T_sfc;
         [T_sfc, ~, ok] = icemodel.numerics.fsearchzero(@surface_residual, ...
            old, tair-50, tair+50, tair, tol);
         if ~ok
            break
         end
         if abs(old - T_sfc) < tol && abs(surface_residual(T_sfc)) < seb_tol
            ok_cpl = true;
            break
         end
      end
   end

   % Hitting max coupling iterations without ok_cpl is a substep fail
   ok = ok && ok_cpl;

   if ~ok && debug
      dumpSebSolveFailure(seb_solver, iter, T_sfc_old, T_sfc, tair, Qsi, ...
         Qli, albedo, wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, ro_air_Lv, ...
         liqflag, chi, k_eff, T_ice, dz, ro_sfc, snow_depth, ok_cpl, opts);
   end

   % Nested surface residual function for convenience
   function f = surface_residual(T_sfc_local)
      f = icemodel.surface.surface_energy_balance_residual(T_sfc_local, tair, ...
         Qsi, Qli, albedo, wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, ...
         ro_air_Lv, liqflag, chi, T_ice, k_eff, dz, ro_sfc, snow_depth, opts);
   end
end

function dumpSebSolveFailure(solver, iter, T_sfc_old, T_sfc, tair, Qsi, Qli, ...
      albedo, wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, ro_air_Lv, liqflag, ...
      chi, k_eff, T_ice, dz, ro_sfc, snow_depth, ok_cpl, opts)
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
   debug_state.br_coefs = br_coefs;
   debug_state.ro_air_Lv = ro_air_Lv;
   debug_state.liqflag = liqflag;
   debug_state.chi = chi;
   debug_state.k_eff = k_eff;
   debug_state.T_ice = T_ice;
   debug_state.dz = dz;
   debug_state.ro_sfc = ro_sfc;
   debug_state.snow_depth = snow_depth;
   debug_state.turbulent_flux_scheme = char(opts.turbulent_flux_scheme);
   debug_state.ok_cpl = ok_cpl;

   save(debug_file, 'debug_state');

   icemodel.surface.dump_turbulent_heat_flux_debug_state( ...
      "sebsolve_nonconvergence", T_sfc_old, T_sfc, tair, Qsi, Qli, albedo, ...
      wspd, ppt, tppt, psfc, De, ea_atm, chi, ro_air_Lv, br_coefs, liqflag, ...
      T_ice, k_eff, dz, ro_sfc, snow_depth, opts);
end
