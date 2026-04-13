function results = validate_thf_cases(options)
   %VALIDATE_THF_CASES Run and optionally plot KANM/KANL THF validation cases.
   %
   %  results = validate_thf_cases()
   %  results = validate_thf_cases(make_plots=false)
   %  results = validate_thf_cases(case_idx=1:4)
   %  results = validate_thf_cases(group_idx=1:2)
   %  results = validate_thf_cases(case_idx=1:4, group_idx=1:2)
   %  results = validate_thf_cases(scalar_comparison=true)
   %
   % This wrapper calls run_thf_cases to execute the validation suite and
   % save the results, then optionally calls plot_thf_cases to produce and
   % save the figures.
   %
   % Options:
   %   make_plots         - produce and save figures (default: true)
   %   scalar_comparison  - also run with output_profile='diagnostic' and call
   %                        plot_thf_scalar_comparison (default: false)
   %   case_idx           - indices into the flat cases list for run_thf_cases
   %   group_idx          - indices into the sorted group list for plot_thf_cases
   %
   % See also: run_thf_cases, plot_thf_cases, plot_thf_scalar_comparison

   arguments
      options.make_plots (1,1) logical = true
      options.scalar_comparison (1,1) logical = false
      options.case_idx (1,:) double = []
      options.group_idx (1,:) double = []
   end

   results = run_thf_cases(case_idx=options.case_idx);

   if options.make_plots
      plot_thf_cases(group_idx=options.group_idx);
   end

   if options.scalar_comparison
      run_thf_cases(case_idx=options.case_idx, output_profile='diagnostic');
      plot_thf_scalar_comparison(group_idx=options.group_idx);
   end
end
