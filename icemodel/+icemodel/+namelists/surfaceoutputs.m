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
   %
   % configureRun composes vars1 from these. It keeps the composition because
   % the run, not the profile, decides model-specific additions such as
   % df_rof for icemodel but not skinmodel, and where they sit in the order.
   %
   % Adding a standard channel here reaches both profiles, so the diagnostic
   % profile cannot silently fall short of the standard one.

   if nargin == 0
      kind = 'standard';
   end

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
      otherwise
         error('icemodel:namelists:surfaceoutputs:kind', ...
            'unsupported surface-output kind: %s', kind)
   end
end
