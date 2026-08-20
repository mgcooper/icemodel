function Lv_or_Ls = latent_enthalpy_switch(f_liq, N)
   %LATENT_ENTHALPY_SWITCH Return Ls or Lv based on surface phase state.
   %
   %  Lv_or_Ls = icemodel.vapor.latent_enthalpy_switch(f_liq)
   %  Lv_or_Ls = icemodel.vapor.latent_enthalpy_switch(f_liq, N)
   %
   % Returns an array the same shape as f_liq (or a column vector of length N
   % when N is given) containing the specific latent heat of sublimation (Ls)
   % for dry/cold cells and the latent heat of vaporization (Lv) for wet cells.
   % The phase-switch threshold is f_liq_phase_switch_threshold, the same
   % threshold the SEB and column solver stacks use, so the latent-heat choice
   % is consistent across the model.
   %
   % The calling functions assign the result to a variable named `Lv`, because
   % to the caller the value is the active latent heat for vapor exchange.
   % That value is Ls or Lv, depending on the local f_liq state.
   %
   % Shape: when you omit N, the output has the same size as f_liq. It then
   % works for column-vector inputs (the primary production path) and for 2-D
   % inputs such as the [JJ × numsteps] arrays that icemodel.postprocess and
   % the diagnostic routines use. When you supply N, the output is a column
   % vector of length N, for callers that compute JJ before the call.
   %
   % Inputs
   %   f_liq - Liquid fraction array of any shape [-].
   %   N     - Optional output column-vector length.  When omitted the output
   %           matches size(f_liq).
   %
   % Output
   %   Lv_or_Ls - Latent heat array [J kg^-1], same size as f_liq (or [N × 1]).
   %              Ls for cells where f_liq <= f_liq_phase_switch_threshold,
   %              Lv for cells where f_liq >  f_liq_phase_switch_threshold.
   %
   % See also: icemodel.column.bulk_enthalpy,
   %           icemodel.column.assemble_enthalpy_system,
   %           icemodel.timestepping.acceptsubstep
   %
   %#codegen

   persistent Ls Lv f_liq_phase_switch_threshold
   if isempty(Ls)
      [Ls, Lv] = icemodel.physicalConstant('Ls', 'Lv');
      f_liq_phase_switch_threshold = icemodel.parameterLookup( ...
         'f_liq_phase_switch_threshold');
   end

   % Default to latent heat of sublimation (dry/cold ice).
   % Match the shape of f_liq unless the caller explicitly requests a column
   % vector of length N.
   if nargin < 2
      Lv_or_Ls = Ls * ones(size(f_liq));
   else
      Lv_or_Ls = Ls * ones(N, 1);
   end

   % Switch to latent heat of vaporization for wet cells.
   wet = f_liq > f_liq_phase_switch_threshold;
   if any(wet(:))
      Lv_or_Ls(wet) = Lv;
   end
end
