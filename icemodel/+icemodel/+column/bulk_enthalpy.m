function [H, dHdT, dFdT, dVdT] = bulk_enthalpy(T, f_ice, f_liq, f_wat, ...
      ro_vap, dro_vapdT)
   %BULK_ENTHALPY Compute the solver-state bulk enthalpy [J m^-3].
   %
   %  H = icemodel.column.bulk_enthalpy(T, f_ice, f_liq, f_wat)
   %  H = icemodel.column.bulk_enthalpy(T, f_ice, f_liq, f_wat, ro_vap)
   %  [H, dHdT, dFdT] = icemodel.column.bulk_enthalpy(...)
   %  [H, dHdT, dFdT, dVdT] = icemodel.column.bulk_enthalpy(...)
   %
   % This helper returns the enthalpy-like bulk state measure used by the
   % column solver. The exact mixture enthalpy relative to Tf can be written
   % schematically as
   %
   %   H(T) = integral_Tf^T [dH_sensible/dT ...
   %      + ro_liq * Lf * df_liq/dT ...
   %      + L_v/s * f_air * d(ro_vap)/dT] dT.
   %
   % The strict primitive of that expression carries reference terms such as
   % f_liq(Tf) and ro_vap(Tf). The solver, however, only uses this helper
   % through H - H_old. Under the current Picard solve assumptions, the
   % omitted reference terms are fixed for a node over the nonlinear update
   % and therefore cancel in the subtraction.
   %
   % In that sense, this helper should be read as the solver-consistent bulk
   % enthalpy measure, not as a claim that the returned H is the uniquely
   % normalized thermodynamic enthalpy of the mixture in every reference
   % convention.
   %
   % The current implementation assumes:
   %  - the liquid-fraction curve is evaluated at fixed f_wat within the solve
   %  - vapor heat is referenced to saturation over the active phase
   %  - the vapor contribution is accurate to the degree that the saturation
   %    reference assumption is accurate for the modeled cell state
   %
   % The returned derivatives are:
   %  dHdT - bulk sensible-heat coefficient [J m^-3 K^-1]
   %  dFdT - liquid-fraction derivative df_liq/dT [K^-1]
   %  dVdT - optional vapor-enthalpy derivative term [J m^-3 K^-1]
   %
   %#codegen

   persistent Tf cv_ice cv_liq roLf Ls Lv f_liq_phase_switch_threshold
   if isempty(Tf)
      [Tf, cv_ice, cv_liq, roLf, Ls, Lv] = icemodel.physicalConstant( ...
         'Tf', 'cv_ice', 'cv_liq', 'roLf', 'Ls', 'Lv');
      f_liq_phase_switch_threshold = icemodel.parameterLookup( ...
         'f_liq_phase_switch_threshold');
   end

   % Option to ignore vapor heat.
   if nargin < 5
      ro_vap = 0;
      dro_vapdT = 0;
   elseif nargin < 6
      dro_vapdT = 0;
   end

   % Air fraction
   f_air = 1.0 - f_liq - f_ice;

   % Vapor heat wrt sublimation [J m^-3]
   roLv = ro_vap * Ls;

   % Vapor heat wrt evaporation (switch to Lv for wet cells)
   wet = f_liq > f_liq_phase_switch_threshold;
   if any(wet)
      roLv(wet) = ro_vap(wet) * Lv;
   end

   % Temperature derivative of the sensible term [J m^-3 K^-1].
   dHdT = cv_ice * f_ice + cv_liq * f_liq;

   % Temperature derivative of the liquid fraction [K^-1].
   if nargout > 2
      dFdT = icemodel.column.liquid_fraction_derivative(T, f_ice, f_liq, f_wat);
   else
      dFdT = [];
   end

   % Total (bulk) enthalpy [J m^-3].
   H = ( ...
      dHdT .* (T - Tf) ...         % sensible heat relative to Tf
      + roLf .* f_liq ...          % latent heat wrt the dry reference state
      + roLv .* f_air ...          % vapor heat at saturation
      );

   % Optional vapor-enthalpy derivative [J m^-3 K^-1].
   if nargout > 3
      dVdT = Ls * dro_vapdT .* f_air;
      if any(wet)
         dVdT(wet) = Lv * dro_vapdT(wet) .* f_air(wet);
      end
   else
      dVdT = [];
   end
end
