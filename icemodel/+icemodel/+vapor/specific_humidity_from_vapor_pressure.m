function q = specific_humidity_from_vapor_pressure(e, p)
   %SPECIFIC_HUMIDITY_FROM_VAPOR_PRESSURE Convert vapor pressure to q.
   %
   % q = epsilon * e / (p - (1 - epsilon) * e)
   %
   % where e is vapor pressure, p is total pressure, and epsilon is the
   % molecular-weight ratio Rd / Rv.

   %#codegen

   persistent epsilon
   if isempty(epsilon)
      epsilon = icemodel.physicalConstant('epsilon');
   end

   q = epsilon * e / (p - (1 - epsilon) * e);
end
