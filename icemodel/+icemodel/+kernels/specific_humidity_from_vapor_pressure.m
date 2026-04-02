function q = specific_humidity_from_vapor_pressure(vap_press, Pa)
   %SPECIFIC_HUMIDITY_FROM_VAPOR_PRESSURE Convert vapor pressure to q.
   %
   % q = epsilon * e / (p - (1 - epsilon) * e)
   %
   % where e is vapor pressure, p is total pressure, and epsilon is the
   % molecular-weight ratio R_d / R_v.

   %#codegen

   persistent epsilon
   if isempty(epsilon)
      epsilon = icemodel.physicalConstant('epsilon');
   end

   q = epsilon * vap_press / (Pa - (1 - epsilon) * vap_press);
end
