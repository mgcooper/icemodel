function [x, ok, iter] = complexstep(f, x0)
   %COMPLEXSTEP Find a scalar root with a complex-step Newton iteration.
   %
   % The Newton iterate itself remains on the real axis. The complex-step
   % perturbation is used only to estimate the derivative at a real state.
   %
   %#codegen

   persistent h dh tol maxiter
   if isempty(tol)
      h = 1e-10;
      dh = 1i * h;
      tol = 1e-3;
      maxiter = 100;
   end

   ok = false;
   old = real(x0);
   for iter = 1:maxiter

      % The residual contract must stay real on the real axis.
      f_old = f(old);
      if abs(imag(f_old)) > imag_tol(f_old)
         error('icemodel:ComplexStepResidualNotReal', ...
            ['complexstep residual must stay real on the real axis; ', ...
            'got imag(f(x)) = %.3e at x = %.15g.'], imag(f_old), old);
      end

      % Estimate the derivative with a complex perturbation and take one
      % real-valued Newton step.
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

function tol = imag_tol(f_old)

   persistent tol_factor
   if isempty(tol_factor)
      tol_factor = 100;
   end

   tol = tol_factor * eps(max(1.0, abs(real(f_old))));
end
