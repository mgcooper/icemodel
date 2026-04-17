function [H, dHdT, dFdT, dLdT, dVdT] = bulk_enthalpy(T, f_ice, f_liq, f_wat, ...
      ro_vap, dro_vapdT)
   %BULK_ENTHALPY Compute the solver-state bulk enthalpy [J m-3].
   %
   %  H = icemodel.column.bulk_enthalpy(T, f_ice, f_liq, f_wat)
   %  H = icemodel.column.bulk_enthalpy(T, f_ice, f_liq, f_wat, ro_vap)
   %  [H, dHdT, dFdT] = icemodel.column.bulk_enthalpy(...)
   %  [H, dHdT, dFdT, dLdT, dVdT] = icemodel.column.bulk_enthalpy(...)
   %
   % This helper returns the enthalpy-like bulk state measure used by the
   % column solver. The exact mixture enthalpy relative to Tf is defined as:
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
   %  - vapor heat is referenced to saturation vapor pressure
   %  - the vapor contribution is accurate to the degree that the saturation
   %    reference assumption is accurate for the modeled cell state
   %
   % The returned derivatives are:
   %  dHdT - bulk sensible-heat coefficient [J m-3 K-1]
   %  dFdT - liquid-fraction derivative df_liq/dT [K-1]
   %  dLdT - liquid-enthalpy derivative [J m-3 K-1]
   %  dVdT - optional vapor-enthalpy derivative term [J m-3 K-1]
   %
   % The return signature is dHdT, dFdT because downstream call sites use these
   % two terms rather than dLdT and dVdT directly.
   %
   %#codegen

   persistent Tf cv_ice cv_liq roLf
   if isempty(Tf)
      [Tf, cv_ice, cv_liq, roLf] = icemodel.physicalConstant( ...
         'Tf', 'cv_ice', 'cv_liq', 'roLf');
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

   % Phase-aware latent heat: Ls for dry/cold cells, Lv for wet cells.
   Lv = icemodel.vapor.latent_enthalpy_switch(f_liq);

   % Vapor heat [J m-3], using the phase-appropriate latent heat.
   roLv = ro_vap .* Lv;

   % Temperature derivative of the sensible term [J m-3 K-1].
   dHdT = cv_ice * f_ice + cv_liq * f_liq;

   % Total (bulk) enthalpy [J m-3].
   H = ( ...
      dHdT .* (T - Tf) ...         % sensible heat relative to Tf
      + roLf .* f_liq ...          % latent heat wrt the dry reference state
      + roLv .* f_air ...          % vapor heat at saturation
      );

   % Optional liquid-fraction derivative [K-1].
   if nargout > 2
      dFdT = icemodel.column.liquid_fraction_derivative(T, f_ice, f_liq, f_wat);
   else
      dFdT = [];
   end

   % Optional liquid-enthalpy derivative [J m-3 K-1].
   if nargout > 3
      dLdT = roLf * dFdT;
   else
      dLdT = [];
   end

   % Optional vapor-enthalpy derivative [J m-3 K-1].
   % Lv already encodes the phase-aware latent heat (Ls or Lv per cell).
   if nargout > 4
      dVdT = Lv .* dro_vapdT .* f_air;
   else
      dVdT = [];
   end
end
