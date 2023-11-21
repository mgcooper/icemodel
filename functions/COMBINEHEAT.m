function [f_liq_C, f_ice_C, T_C, Sc_C, Sp_C, d_drn, d_liq] = COMBINEHEAT( ...
      f_liq, f_ice, T, Sc, Sp, Tf, TL, cv_ice, cv_liq, Lf, ro_ice, ro_liq, ...
      fcp, j1, j2, d_drn, d_liq, dz, fopts)
   %COMBINEHEAT Combine layers conserving mass, enthalpy, and absorbed radiation
   %
   % j1 is the thin layer that is removed, the combined values are put into j2.
   %
   % T = control volume temperature.
   % Tf = freezing point temperature.
   % Td = depression temperature, Td = Tf - T.
   %
   % f_liq = liquid fraction, volume of liquid water per cv volume
   % f_ice = frozen fraction, volume of frozen water per cv volume
   %
   % f_wat = water fraction, volume of melted ice + liquid water per cv volume:
   % f_wat = f_liq + f_ice * ro_ice / ro_liq
   %
   % ro_liq = liquid water intrinsic density, 1000 kg m-3
   % ro_ice = frozen water intrinsic density, 917 kg m-3
   %
   % The temperature is a function of the liquid fraction:
   %
   % T = Tf - sqrt((f_wat / f_liq - 1)) / fcp;
   %
   % where fcp is the "freezing curve parameter" that controls the slope of
   % the f_liq = f(T) relationship in the mushy zone.
   %
   % Rearrange that for depression temp as a function of f_liq:
   %
   % Td = sqrt((f_wat / f_liq - 1)) / fcp;
   %
   % This function solves for Td by equating the sum of the enthalpies of two
   % cv's given their respective Td's to their combined enthalpy given a common Td.
   %
   % The Td solution is then used to compute the combined layer f_liq:
   %
   % f_liq = f_wat / (1 + (fcp * Td)^2)
   %
   % Subscript _C is used for the combined layer quantities.
   % Subscript _12 is used for the arithmetic sum of the separate layer quantities
   %
   %#codegen

   % Combine absorbed solar radiation [W m-3]
   Sc_C = Sc(j1) + Sc(j2);
   Sp_C = Sp(j1) + Sp(j2);

   % Depression temperature of each cv
   Td_1 = Tf-T(j1);
   Td_2 = Tf-T(j2);

   % Compute water mass (ice + liquid water)
   m_wat_1 = (ro_ice*f_ice(j1) + ro_liq*f_liq(j1)) * dz;
   m_wat_2 = (ro_ice*f_ice(j2) + ro_liq*f_liq(j2)) * dz;
   m_wat_C = (ro_ice*(f_ice(j1)+f_ice(j2)) + ro_liq*(f_liq(j1)+f_liq(j2))) * dz;

   % Compute combined temperature assuming T1 != T2 but neither CV is melting.
   % This temperature also applies to the case where T1 == T2, melting or not.
   T_C = (T(j1)*m_wat_1 + T(j2)*m_wat_2) / m_wat_C;
   Td_C = Tf - T_C;

   % If either CV is melting, compute combined temperature of dry and wet ice.
   if T(j1) ~= T(j2) && (T(j1) > TL || T(j2) > TL)

      % Equate the sum of the enthalpies of the separate cv's to that of the
      % combined cv. Use m_liq_C to express the equation as a third-order
      % polynomial of the combined Td. Define the coefficients:

      f = -cv_ice * ( f_ice(j1) * Td_1 + f_ice(j2) * Td_2 ) + ...
         -cv_liq * ( f_liq(j1) * Td_1 + f_liq(j2) * Td_2 ) + ...
         Lf * ro_liq * ( f_liq(j1) + f_liq(j2) );
      g = cv_ice * (f_ice(j1) + f_ice(j2)) + cv_liq * (f_liq(j1) + f_liq(j2));

      a = 1;
      b = f / g;
      c = 1 / fcp^2;
      d = c * (b - Lf * m_wat_C / g / dz);

      F = @(Td_C) a*Td_C^3 + b*Td_C^2 + c*Td_C + d;

      % solve using an endpoint constrained range
      [Td_C, ok] = fsearchzero(F, Td_C, min(Td_1,Td_2), max(Td_1,Td_2), Td_C, 1e-4);

      %[Td_C, ~, ok] = fzero(F, [min(Td_1,Td_2) max(Td_1,Td_2)], fopts);

      % Check if the function successfully found a zero
      if ok == 1
         T_C = Tf - Td_C;

      elseif ok <= 0
         % fall back to the mass-weighted average temperature

         % activate this to have an error message
         % [errStatus, errMsg] = errorcheck(ok);
      end

      % Set T=Tf for a fully-melted node and ICEMF will remove it
      if Td_C < 0
         T_C  = (T(j1) * m_wat_1 + T(j2) * m_wat_2) / m_wat_C;
         Td_C = Tf - T_C;
      end
   end

   % Invert T_dC to obtain f_liq from the f_liq = f(f_wat, Tdc) relationship.
   f_wat_C = m_wat_C / (ro_liq * 2 * dz);
   f_liq_C = f_wat_C / (1.0 + (fcp * Td_C) ^ 2.0); % eq 67, Jordan
   f_ice_C = (f_wat_C - f_liq_C) * ro_liq / ro_ice;
   
   % positive d_liq means ice melted, negative means liquid refroze.
   d_liq(j1) = d_liq(j1) + 2 * f_liq_C - (f_liq(j1) + f_liq(j2));

   % if there is less water in the combined layer then water drained
   % d_drn(j1) = d_drn(j1) + max((f_liq(j1)+f_liq(j2))-f_liq_C-d_liq(j1), 0);
   d_drn(j1) = d_drn(j1) + max(f_liq_C / 2, 0);
end
