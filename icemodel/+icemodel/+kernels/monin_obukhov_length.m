function L = monin_obukhov_length(u_star, theta_air, q_air, theta_star, q_star)
   %MONIN_OBUKHOV_LENGTH Return the Monin-Obukhov stability length.
   %
   % L = u_*^2 theta_v / (kappa g theta_v*)
   %
   % The virtual-potential-temperature corrections are written here in the
   % humidity form used by the bulk-MO scheme:
   %   theta_v  = theta * (1 + c_q q)
   %   theta_v* = theta* * (1 + c_q q*)

   %#codegen

   persistent epsilon kappa gravity
   if isempty(epsilon)
      [epsilon, kappa, gravity] = icemodel.physicalConstant( ...
         'epsilon', 'kappa', 'gravity');
   end

   moisture_coeff = (1 - epsilon) / epsilon;
   numer = u_star ^ 2 * theta_air * (1 + moisture_coeff * q_air);
   denom = gravity * kappa * theta_star * (1 + moisture_coeff * q_star);

   if abs(denom) < 1e-12
      denom = denom + sign_or_one(denom) * 1e-12;
   end

   L = numer / denom;
end
