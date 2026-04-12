function dfdx = complex_step_derivative(f, x, h)
   %COMPLEX_STEP_DERIVATIVE Estimate a scalar derivative with a complex step.
   %
   %  dfdx = icemodel.numerics.complex_step_derivative(f, x)
   %  dfdx = icemodel.numerics.complex_step_derivative(f, x, h)
   %
   % Evaluates the derivative of a scalar-valued function handle `f` at the
   % real state `x` using the complex-step identity
   %
   %   df/dx ≈ imag(f(x + 1i*h)) / h
   %
   % which avoids the subtractive cancellation present in a forward-difference
   % quotient while keeping the Newton iterate itself on the real axis.
   %
   %#codegen

   if nargin < 3
      h = 1e-10;
   end

   % Perturb only the derivative evaluation, not the iterate itself.
   dfdx = imag(f(x + 1i * h)) / h;
end
