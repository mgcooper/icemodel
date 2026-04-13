function Lv_or_Ls = latent_enthalpy_switch(f_liq, N)
   %LATENT_ENTHALPY_SWITCH Return Ls or Lv based on surface phase state.
   %
   %  Lv_or_Ls = icemodel.vapor.latent_enthalpy_switch(f_liq)
   %  Lv_or_Ls = icemodel.vapor.latent_enthalpy_switch(f_liq, N)
   %
   % Returns a column vector of length N containing the specific latent heat
   % of sublimation (Ls) for dry/cold cells and the latent heat of
   % vaporization (Lv) for wet cells.  The phase switch threshold is the same
   % f_liq_phase_switch_threshold used elsewhere in the SEB and column solver
   % stacks so that the latent-heat choice is consistent across the model.
   %
   % The calling functions assign the result to a variable named `Lv` because
   % from the caller's perspective the value is the active latent heat for vapor
   % exchange, which may be Ls or Lv depending on the local f_liq state.
   %
   % Inputs
   %   f_liq - Liquid fraction vector or scalar [-].
   %   N     - Optional output length.  Defaults to numel(f_liq).  Pass an
   %           explicit N when the caller already knows the column length to
   %           avoid a redundant numel call.
   %
   % Output
   %   Lv_or_Ls - Latent heat vector [J kg^-1], length N.
   %              Ls for cells where f_liq <= f_liq_phase_switch_threshold,
   %              Lv for cells where f_liq >  f_liq_phase_switch_threshold.
   %
   % See also: icemodel.column.bulk_enthalpy,
   %           icemodel.column.assemble_enthalpy_system,
   %           icemodel.timestepping.updatesubstep
   %
   %#codegen

   persistent Ls Lv f_liq_phase_switch_threshold
   if isempty(Ls)
      [Ls, Lv] = icemodel.physicalConstant('Ls', 'Lv');
      f_liq_phase_switch_threshold = icemodel.parameterLookup( ...
         'f_liq_phase_switch_threshold');
   end

   % Option to pass in N or derive it from f_liq.
   if nargin < 2
      N = numel(f_liq);
   end

   % Default to latent heat of sublimation (dry/cold ice).
   Lv_or_Ls = Ls * ones(N, 1);

   % Switch to latent heat of vaporization for wet cells.
   wet = f_liq > f_liq_phase_switch_threshold;
   if any(wet)
      Lv_or_Ls(wet) = Lv;
   end
end
