function fields = surfaceoutputs(kind)
   %SURFACEOUTPUTS Return surface (ice1) output channel names.
   %
   %  fields = icemodel.namelists.surfaceoutputs()
   %  fields = icemodel.namelists.surfaceoutputs(kind)
   %
   % KIND is one of:
   %   'standard'          the channels every standard and diagnostic run
   %                       writes (the default)
   %   'diagnostic_suffix' the channels the diagnostic profile appends
   %   'icemodel_diagnostic_suffix' the diagnostic channels that only the
   %                       full column model writes
   %   'additive_diagnostic' the subset of the solver diagnostics that
   %                       retiming sums over a bin; isIncrementChannel
   %                       derives additivity from it
   %
   % configureRun composes vars1 from these lists. It adds the model-specific
   % channels, such as df_rof for icemodel but not skinmodel, and sets where
   % they sit in the order.
   %
   % A channel added to the standard list appears in both the standard and
   % diagnostic profiles.

   if nargin == 0
      kind = 'standard';
   end

   % The additive subset of the solver diagnostics. Retiming sums these
   % channels over a bin; the residual and iteration diagnostics
   % (cpl_iters, cpl_res, seb_res, Tice_numiter) are per-step records,
   % not sums. icemodel.isIncrementChannel derives additivity from this
   % list.
   additive_diagnostic = {'cpl_recovery_count'};

   switch lower(char(kind))
      case 'standard'
         fields = {'Tsfc', 'Qm', 'Qe', 'Qh', 'Qc', 'chi', 'balance', ...
            'dt_sum', 'Tsfc_converged', 'Tice_converged', 'Tice_numiter'};
      case 'diagnostic_suffix'
         fields = {'n_subfail', 'ea_atm', 'br_coefs_gamma', ...
            'br_coefs_b1_num', 'br_coefs_b2_num', 'hv_atm', 'ro_sfc', ...
            'thf_es_sfc', 'thf_stability_factor', 'thf_z0m', 'thf_z0h', ...
            'thf_z0q', 'thf_u_star', 'thf_L', 'thf_Re', 'thf_numiter', ...
             'thf_scalar_exchange_Qh', 'thf_scalar_exchange_Qe'};
      case 'icemodel_diagnostic_suffix'
         % Solver-health observability from the coupler diag struct,
         % reported from the ACCEPTED solve (sentinels after a step of
         % only forced advances). cpl_recovery_count counts accepted
         % conservative recoveries per forcing step.
         fields = [{'cpl_iters', 'cpl_res', 'seb_res'}, ...
            additive_diagnostic];
      case 'additive_diagnostic'
         fields = additive_diagnostic;
      otherwise
         error('icemodel:namelists:surfaceoutputs:kind', ...
            'unsupported surface-output kind: %s', kind)
   end
end
