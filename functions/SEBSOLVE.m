function Ts = SEBSOLVE(Ta, Qsi, Qli, albedo, wspd, Pa, De, ea, cv_air, ...
      emiss, SB, Tf, chi, roL, scoef, liqflag, Ts, T, k_eff, dz, solver)
   %SEBSOLVE solve the surface energy balance for the skin temperature
   %#codegen

   persistent tol maxiter
   if isempty(tol)
      tol = 1e-3;
      maxiter = 100;
   end

   switch solver

      case 1 % Call SFCTEMP directly

         for iter = 1:maxiter

            old = Ts;
            [Ts, ok] = SFCTEMP(Ta, Qsi, Qli, albedo, wspd, Pa, De, ...
               ea, cv_air, emiss, SB, Tf, chi, roL, scoef, ...
               CONDUCT(k_eff, T, dz, old), liqflag);

            if not(ok) || abs(old - Ts) < tol
               break
            end
         end

      case 2 % Use nested function to create the function handle

         for iter = 1:maxiter

            old = Ts;
            Qc = CONDUCT(k_eff, T, dz, old);
            % [Ts, ok] = newtrhap(@fSEB, Ta);
            [Ts, ~, ok] = fsearchzero(@fSEB, old, Ta-50, Ta+50, Ta, tol);

            if not(ok) || abs(old - Ts) < tol
               break
            end
         end

      case 3 % Pass the function handle directly to root finder
         for iter = 1:maxiter
            old = Ts;
            [Ts, ~, ok] = fsearchzero(@(Ts) CONDUCT(k_eff, T, dz, old) ...
               + chi * (1.0 - albedo) * Qsi + emiss * (Qli - SB * Ts ^ 4) ...
               + cv_air * De * (Ta - Ts) * STABLEFN(Ta, Ts, wspd, scoef) ...
               + roL * De * 0.622 / Pa * (ea - VAPPRESS(Ts, Tf, liqflag)) ...
               * STABLEFN(Ta, Ts, wspd, scoef), ...
               old, Ta-50, Ta+50, Ta, tol);

            if not(ok) || abs(old - Ts) < tol
               break
            end
         end
   end

   % Note: must use nested function to capture updated Qc on each iteration.
   function f = fSEB(Ts)
      f = chi * (1.0 - albedo) * Qsi + emiss * (Qli - SB * Ts ^ 4) ...
         + Qc + cv_air * De * (Ta - Ts) * STABLEFN(Ta, Ts, wspd, scoef) ...
         + roL * De * 0.622 / Pa * (ea - VAPPRESS(Ts, Tf, liqflag)) ...
         * STABLEFN(Ta, Ts, wspd, scoef);
      % + cp_liq * ppt * Tppt; % ppt in kg/m2/s
   end
end

%% Newton-rhapson
function [x, ok, iter] = newtrhap(f, x0) %#ok<*DEFNU>
   %NEWTRHAP Find root of nonlinear function using the Newton Rhapson method.
   %
   % Note: this uses the complex step method.

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
