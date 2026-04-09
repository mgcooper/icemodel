function results = validate_thf_cases(options)
   %VALIDATE_THF_CASES Run and optionally plot KANM/KANL THF validation cases.
   %
   %  results = validate_thf_cases()
   %  results = validate_thf_cases(make_plots=false)
   %  results = validate_thf_cases(case_idx=1:4)
   %  results = validate_thf_cases(group_idx=1:2)
   %  results = validate_thf_cases(case_idx=1:4, group_idx=1:2)
   %
   % This wrapper calls run_thf_cases to execute the validation suite and
   % save the results, then optionally calls plot_thf_cases to produce and
   % save the figures.
   %
   % make_plots  - produce and save figures (default: true)
   % case_idx    - indices into the flat cases list for run_thf_cases (default: all)
   % group_idx   - indices into the sorted group list for plot_thf_cases (default: all)
   %
   % See also: run_thf_cases, plot_thf_cases

   arguments
      options.make_plots (1,1) logical = true
      options.case_idx (1,:) double = []
      options.group_idx (1,:) double = []
   end

   results = run_thf_cases(case_idx=options.case_idx);

   if options.make_plots
      plot_thf_cases(group_idx=options.group_idx);
   end
end
