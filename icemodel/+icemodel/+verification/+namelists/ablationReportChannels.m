function fields = ablationReportChannels(kind)
   %ABLATIONREPORTCHANNELS Model channels read by the ablation report.
   %
   %  fields = icemodel.verification.namelists.ablationReportChannels()
   %  fields = icemodel.verification.namelists.ablationReportChannels(kind)
   %
   % KIND is 'all', 'components', 'grid', or 'ledger'. These channels are the
   % subset of model diagnostics read by
   % icemodel.verification.report.buildAblationEvaluationReport. They are
   % grouped by the report table that reads them.
   %
   % policy.required_model_fields defines the full cohort run schema from
   % icemodel.namelists.budgetoutputs('all') and
   % icemodel.namelists.cumulativeoutputs().
   %
   % It is not the only place the names appear. componentTable and
   % gridTranslationRows also read the same columns by name. Renaming a
   % channel therefore means editing this list, those two readers, and
   % icemodel.namelists.budgetoutputs in one change.
   %
   % Adding a channel here makes every saved cohort that lacks it invalid for
   % report building. Add one only together with the report code that reads it.
   %
   % See also: icemodel.namelists.budgetoutputs,
   %  icemodel.verification.namelists.promiceAblationPolicy,
   %  icemodel.verification.report.buildAblationEvaluationReport

   if nargin == 0
      kind = 'all';
   end

   % componentTable sums the signed physical solid-mass components.
   component_fields = [ ...
      "mass_budget_phase_solid_mwe", ...
      "mass_budget_vapor_solid_mwe"];

   % gridTranslationRows reads the quantized remesh event counters.
   grid_fields = [ ...
      "mass_budget_top_deletion_count", ...
      "mass_budget_top_deletion_height_m", ...
      "mass_budget_interior_merge_count"];

   % icemodel.verification.helpers.ablationLedgerIncrements reads these to
   % calculate the ablation terms for each interval. compareAblation requires
   % every budget and cumulative output in policy.required_model_fields before
   % it compares the data.
   %
   % Keep these out of 'all'. 'all' is the report's schema gate, and adding a
   % channel to it makes every cohort that lacks it invalid for report
   % building.
   ledger_fields = [ ...
      "mass_budget_phase_solid_mwe", ...
      "mass_budget_vapor_solid_mwe", ...
      "mass_budget_top_export_solid_mwe", ...
      "mass_budget_top_export_liquid_mwe"];

   switch char(kind)
      case 'all'
         fields = [component_fields, grid_fields];
      case 'components'
         fields = component_fields;
      case 'grid'
         fields = grid_fields;
      case 'ledger'
         fields = ledger_fields;
      otherwise
         error( ...
            'icemodel:verification:namelists:ablationReportChannels:kind', ...
            'KIND must be all, components, grid, or ledger')
   end
end
