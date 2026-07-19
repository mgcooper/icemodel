function summary = plotFirnArtifacts(kwargs)
   %PLOTFIRNARTIFACTS Plot staged firn-family artifacts for visual QA.
   %
   %  summary = icemodel.verification.plotFirnArtifacts(...)
   %
   % Use this entry point for firn verification families when creating staged
   % artifact diagnostics for visual inspection.
   %
   % See also: icemodel.verification.plotVerificationArtifacts

   arguments
      kwargs.dataset_family (1, :) string ...
         {icemodel.verification.validators.mustBeFirnDatasetFamilySelection} = "all"
      kwargs.case_ids (1, :) string = strings(1, 0)
      kwargs.output_root (1, 1) string = ""
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.input_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = ""
      kwargs.figure_root (1, 1) string = ""
      kwargs.save_figs (1, 1) logical = true
      kwargs.overwrite (1, 1) logical = false
      kwargs.visible (1, 1) logical = false
      kwargs.startdate = ""
      kwargs.enddate = ""
   end

   families = expandFamilies(kwargs.dataset_family);
   summary = icemodel.verification.plotVerificationArtifacts( ...
      dataset_family=families, ...
      case_ids=kwargs.case_ids, ...
      output_root=kwargs.output_root, ...
      evaluation_data_root=kwargs.evaluation_data_root, ...
      input_data_root=kwargs.input_data_root, ...
      icemodel_config_casename=kwargs.icemodel_config_casename, ...
      figure_root=kwargs.figure_root, ...
      save_figs=kwargs.save_figs, ...
      overwrite=kwargs.overwrite, ...
      visible=kwargs.visible, ...
      startdate=kwargs.startdate, ...
      enddate=kwargs.enddate);
end

function families = expandFamilies(families)
   %EXPANDFAMILIES Limit the firn wrapper's "all" selector to firn families.
   if any(families == "all")
      families = icemodel.verification.namelists.firndatasetfamily();
   end
   families = unique(families, 'stable');
end
