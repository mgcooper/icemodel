function e = vapor_pressure_from_specific_humidity(q, p)
   %SPECIFIC_HUMIDITY_FROM_VAPOR_PRESSURE Convert specific humidity to vapor pressure
   %
   % e = q * p / (epsilon + (1 - epsilon) * q)
   %
   % where e is vapor pressure, q is specific humidity, p is total pressure, and
   % epsilon is the molecular-weight ratio Rd / Rv.

   %#codegen

   persistent epsilon
   if isempty(epsilon)
      epsilon = icemodel.physicalConstant('epsilon');
   end

   e = q .* p ./ (epsilon + (1 - epsilon) .* q);
end
