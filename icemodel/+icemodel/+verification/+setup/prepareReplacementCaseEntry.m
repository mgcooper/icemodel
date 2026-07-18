function entry = prepareReplacementCaseEntry(entry, dataset_family)
   %PREPAREREPLACEMENTCASEENTRY Retain observations but clear prior runtime legs.
   %
   %  entry = icemodel.verification.setup.prepareReplacementCaseEntry( ...
   %     entry, dataset_family)
   %
   % A build_observations=false whole-family replacement needs the requested
   % observation contract, but its prior native/RCM runtime state is not part of
   % the replacement. This helper removes every RCM leg and clears the combined
   % native/evaluation fields used by PROMICE, IMAU, and RetMIP before requested
   % forcing is attached.
   %
   % See also: icemodel.verification.setup.reuseDatasetFamilyCases,
   %  icemodel.verification.setup.colocationSourceLists

   arguments
      entry (1, 1) struct
      dataset_family (1, 1) string
   end

   if ~isfield(entry, 'colocation') || ~isstruct(entry.colocation)
      return
   end

   % Whole-family replacement discards every prior RCM leg; selected sources are
   % rediscovered or rebuilt after the observation-only checkpoint.
   colocation = entry.colocation;
   for source = icemodel.verification.namelists.rcmsources()
      name = char(source);
      if isfield(colocation, name)
         colocation = rmfield(colocation, name);
      end
   end

   % PROMICE, IMAU, and RetMIP store native runtime and evaluation provenance
   % together, so clear only the runtime portion while retaining observations.
   switch dataset_family
      case "promice"
         if isfield(colocation, 'promice') && isstruct(colocation.promice)
            colocation.promice.kind = 'station_met_and_eval';
            colocation.promice.staged = false;
            colocation.promice.eval_staged = true;
            colocation.promice.met_files = strings(1, 0);
            colocation.promice.data_files = strings(1, 0);
            runtime_fields = [ ...
               "forcing_ready", "forcing_ready_reason", ...
               "forcing_complete_windows", "met_skipped_reason"];
            present = runtime_fields(ismember( ...
               runtime_fields, string(fieldnames(colocation.promice))));
            if ~isempty(present)
               colocation.promice = rmfield( ...
                  colocation.promice, cellstr(present));
            end
         end

      case "imau"
         if isfield(colocation, 'imau') && isstruct(colocation.imau)
            colocation.imau.kind = 'hourly_aws_eval';
            colocation.imau.staged = false;
            colocation.imau.eval_staged = true;
            colocation.imau.met_files = strings(1, 0);
            colocation.imau.data_files = strings(1, 0);
            colocation.imau.forcing_ready = false;
            colocation.imau.forcing_ready_reason = ...
               'imau was not requested in forcing_sources';
            if isfield(colocation.imau, 'forcing_complete_windows')
               colocation.imau.forcing_complete_windows = ...
                  colocation.imau.forcing_complete_windows([]);
            end
            if isfield(colocation.imau, 'artifact_metadata')
               colocation.imau = rmfield( ...
                  colocation.imau, 'artifact_metadata');
            end
         end

      case "retmip"
         if isfield(colocation, 'retmip') && isstruct(colocation.retmip)
            colocation.retmip.kind = 'protocol_userdata';
            runtime_fields = [ ...
               "native_source", "native_met_status", ...
               "native_met_skipped_reason", "met_files", "data_files", ...
               "forcing_ready", "forcing_ready_reason", ...
               "forcing_complete_windows", "window"];
            present = runtime_fields(ismember( ...
               runtime_fields, string(fieldnames(colocation.retmip))));
            if ~isempty(present)
               colocation.retmip = rmfield( ...
                  colocation.retmip, cellstr(present));
            end
         end
         if isfield(colocation, 'native_met')
            colocation = rmfield(colocation, 'native_met');
         end
         if isfield(entry, 'observation_variables') ...
               && isstruct(entry.observation_variables) ...
               && isfield(entry.observation_variables, 'native_source')
            entry.observation_variables = rmfield( ...
               entry.observation_variables, 'native_source');
         end
         if isfield(colocation, 'retmip') ...
               && isfield(colocation.retmip, 'protocol_timestep')
            entry.native_timestep = colocation.retmip.protocol_timestep;
         end
   end

   % Public source lists must describe the observation-only replacement state.
   entry.colocation = colocation;
   [forcing_sources, eval_sources] = ...
      icemodel.verification.setup.colocationSourceLists(colocation);
   entry.forcing_sources = cellstr(forcing_sources);
   entry.eval_sources = cellstr(eval_sources);
end
