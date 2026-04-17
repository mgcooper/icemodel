function Lv_or_Ls = latent_enthalpy_switch(f_liq, N)
   %LATENT_ENTHALPY_SWITCH Return Ls or Lv based on surface phase state.
   %
   %  Lv_or_Ls = icemodel.vapor.latent_enthalpy_switch(f_liq)
   %  Lv_or_Ls = icemodel.vapor.latent_enthalpy_switch(f_liq, N)
   %
   % Returns an array the same shape as f_liq (or a column vector of length N
   % when N is given) containing the specific latent heat of sublimation (Ls)
   % for dry/cold cells and the latent heat of vaporization (Lv) for wet cells.
   % The phase switch threshold is the same f_liq_phase_switch_threshold used
   % elsewhere in the SEB and column solver stacks so that the latent-heat
   % choice is consistent across the model.
   %
   % The calling functions assign the result to a variable named `Lv` because
   % from the caller's perspective the value is the active latent heat for vapor
   % exchange, which may be Ls or Lv depending on the local f_liq state.
   %
   % Shape: when N is omitted the output has the same size as f_liq and
   % therefore works correctly for both column-vector inputs (the primary
   % production path) and 2-D inputs such as [JJ × numsteps] arrays used by
   % icemodel.postprocess and diagnostic routines.  When N is supplied the
   % output is a column vector of length N (retained for callers that pre-compute
   % JJ before calling).
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
   %           icemodel.timestepping.updatesubstep
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
