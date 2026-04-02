function s = sign_or_one(x)
   %SIGN_OR_ONE Return the sign of X, treating zero as positive one.
   %
   % This helper preserves the sign of a regularization term while avoiding a
   % zero multiplier when the real part of x vanishes exactly. The real part is
   % used so complex-step perturbations keep the same branch choice as the
   % underlying real-valued state.

   s = sign(real(x));
   if s == 0
      s = 1;
   end
end
