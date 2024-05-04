function [Ts, ok] = SEBSOLVE(Ta, Qsi, Qli, albedo, wspd, ppt, tppt, Pa, De, ...
      ea, cv_air, cv_liq, emiss, SB, Tf, chi, roL, scoef, liqflag, Ts, T, ...
      k_eff, dz, solver)
   %SEBSOLVE solve the surface energy balance for the skin temperature
   %
   % Solver options:
   % 0 = derivative-free, slower but more robust, thus used as a fall back.
   % 1 = newton-rhapson, fast, requires analytic derivative, used as a default.
   % 2 = complex-step, fast and does not require an analytical derivative.
   % -1 = no outer iterations, use this if phase change is not represented
   % explicitly in the model.
   % >2 = experimental.
   %
   % Important programming notes:
   %  - The outer iterations control the convergence of the Ts calculation wrt
   %  the conduction term. Thus old is initialized to Ts outside the outer loop.
   %  - In SFCTEMP or complexstep or any derivative-based method, old must be
   %  initialized to Ta to avoid divergence during spinup, keeping in mind that
   %  the conduction passed into SFCTEMP is computed with old = Ts.
   %  - For a "skinmodel", Ts never exceeds Tf when passed into functions, but
   %  within the iterations of SFCTEMP and when it comes out of SFCTEMP it can
   %  exceed Tf.
   %
   %#codegen

   tol = 1e-3;
   maxiter = 100;
   iterflag = true;
   if solver < 0
      maxiter = 1;
      solver = -solver;
      iterflag = false;
   end

   Ts_old = Ts;
   switch solver

      case 1 % Newton-Rhapson - analytical derivative

         for iter = 1:maxiter
            old = Ts;
            [Ts, ok] = SFCTEMP(Ta, Qsi, Qli, albedo, wspd, Pa, De, ...
               ea, cv_air, emiss, SB, Tf, chi, roL, scoef, liqflag, ...
               CONDUCT(k_eff, T, dz, old), k_eff, T, dz, iterflag);

            if not(ok) || abs(old - Ts) < tol
               break
            end
         end

      case 2 % Complex-step - numerical derivative

         for iter = 1:maxiter
            old = Ts;
            [Ts, ok] = complexstep(@fSEB, Ta);
            if not(ok) || abs(old - Ts) < tol
               break
            end
         end

      otherwise % Derivative-free
         solver = 0;
   end

   % Derivative-free solver
   if solver == 0 || not(ok) || abs(Ts - Ta) > 20

      Ts = Ts_old;

      for iter = 1:maxiter
         old = Ts;
         [Ts, ~, ok] = fsearchzero(@fSEB, old, Ta-50, Ta+50, Ta, tol);
         if not(ok) || abs(old - Ts) < tol
            break
         end
      end
   end

   % Note: nested function captures updated Qc on each iteration.
   function f = fSEB(Ts)
      f = chi * (1.0 - albedo) * Qsi + emiss * (Qli - SB * Ts ^ 4) ...
         + CONDUCT(k_eff, T, dz, Ts) ...
         + cv_air * De * (Ta - Ts) * STABLEFN(Ta, Ts, wspd, scoef) ...
         + roL * De * 0.622 / Pa * (ea - VAPPRESS(Ts, Tf, liqflag)) ...
         * STABLEFN(Ta, Ts, wspd, scoef) ...
         + QADVECT(ppt, tppt, cv_liq);
   end
end

function [x, ok, iter] = complexstep(f, x0)
   %COMPLEXSTEP Find root of nonlinear function using the Newton Rhapson method.
   %
   % Note: this uses the complex step numerical derivative.

   persistent h dh tol maxiter
   if isempty(tol)
      h = 1e-10;
      dh = 1i * h;
      tol = 1e-3;
      maxiter = 100;
   end

   ok = false;
   old = x0;
   for iter = 1:maxiter

      x = old - f(old) / (imag(f(old + dh)) / h);

      if abs(x - old) < tol
         ok = true;
         return

      elseif isnan(x)
         x = x0;
         return
      end
      old = x;
   end
   x = x0; % ok = false
end
