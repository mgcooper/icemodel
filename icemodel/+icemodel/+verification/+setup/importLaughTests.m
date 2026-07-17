function manifest = importLaughTests(laugh_tests_source_dir, kwargs)
   %IMPORTLAUGHTESTS Stage selected Laugh-Tests synthetic snow benchmarks.
   %
   %  manifest = icemodel.verification.setup.importLaughTests(source_dir)
   %  manifest = icemodel.verification.setup.importLaughTests(source_dir, ...
   %     overwrite=true)
   %
   % Laugh cases intentionally expose no forcing_sources/build_observations
   % split: each case is an atomic evaluation/reference bundle and has no staged
   % runtime forcing source that can be attached independently.
   %
   % Inputs
   %  laugh_tests_source_dir     Root of a local Laugh-Tests checkout.
   %  output_root                Base output root; eval artifacts go to
   %                             output_root/eval.
   %  evaluation_data_root       Base evaluation-data root to stage into.
   %  icemodel_config_casename   Config casename used when evaluation_data_root
   %                             is blank. Blank uses repo-local data/eval.
   %  case_id                    Laugh-Tests case to stage. Currently only
   %                             "colbeck1976" is supported.
   %  overwrite                  Refresh setup artifacts when true; protect
   %                             existing staged data when false.
   %  overwrite_family           Replace the family manifest instead of merging.
   %  skip_missing               Record missing source data as skipped cases.
   %  dry_run                    Return the manifest shape without writes.
   %
   % Outputs
   %  manifest   Family manifest struct also written to manifest.json.
   %
   % Role
   %  Setup/update tooling. This function creates or refreshes staged data
   %  under the resolved data/eval/laugh_tests tree and is not part of normal verification
   %  runs.
   %
   % See also icemodel.verification.setup.importEsmSnowmip

   arguments
      laugh_tests_source_dir (1, 1) string = ""
      kwargs.output_root (1, 1) string = ""
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = ""
      kwargs.case_id (1, 1) string ...
         {icemodel.verification.validators.mustBeLaughTestCase} = ...
         icemodel.verification.namelists.caseid("laugh_tests")
      kwargs.overwrite (1, 1) logical = false
      kwargs.overwrite_family (1, 1) logical = false
      kwargs.skip_missing (1, 1) logical = false
      kwargs.dry_run (1, 1) logical = false
   end

   % Name the source family and runnable case once. dataset_family is the staged
   % source folder/manifest family; case_id is the benchmark case inside that
   % family (currently only "colbeck1976").
   dataset_family = "laugh_tests";
   case_id = kwargs.case_id;
   source_status = [];
   if ~kwargs.dry_run
      [laugh_tests_source_dir, source_status] = ...
         icemodel.verification.setup.fetchLaughTests( ...
         cache_dir=laugh_tests_source_dir, strict=~kwargs.skip_missing, ...
         silent=kwargs.skip_missing);
   end

   % Resolve the path to the evaluation data folder through the shared staging
   % root helper. Laugh-Tests has no staged input artifacts, so input_root is
   % intentionally ignored.
   [evaluation_data_root, ~] = ...
      icemodel.verification.setup.resolveStagingRoots( ...
      output_root=kwargs.output_root, ...
      evaluation_data_root=kwargs.evaluation_data_root, ...
      icemodel_config_casename=kwargs.icemodel_config_casename);
   family_root = fullfile(evaluation_data_root, dataset_family);
   if ~kwargs.dry_run
      icemodel.helpers.ensureDirExists(family_root);
   end
   manifest_file = fullfile(family_root, "manifest.json");

   % Stage each requested case into importer state.
   [state, alive, skipped] = ...
      icemodel.verification.setup.stageDatasetFamilyCases( ...
      case_id, emptyState(), ...
      @(id, ~) stageLaughCase(id, laugh_tests_source_dir, family_root, ...
      source_status, kwargs), ...
      skip_missing=kwargs.skip_missing, ...
      warning_id="icemodel:verification:importLaughTests:caseSkipped", ...
      label_callback=@(id, ~) id);

   % Record source provenance once so dry-run and persisted manifests match.
   source_doi = "";
   source_url = "https://github.com/KyleKlenk/Laugh-Tests";
   source_version = "m2_mac_Sept23 validation bundle";
   retrieval_date = string(datetime('today'));

   if kwargs.dry_run
      manifest = icemodel.verification.setup.runDatasetFamilyDryRun( ...
         state, alive, dataset_family=dataset_family, ...
         requested_ids=case_id, skipped=skipped, source_doi=source_doi, ...
         source_url=source_url, source_version=source_version, ...
         retrieval_date=retrieval_date, entry_callback=@caseEntry);
      return
   end

   [manifest, ~] = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family=dataset_family, ...
      manifest_file=manifest_file, requested_ids=case_id, skipped=skipped, ...
      source_doi=source_doi, source_url=source_url, ...
      source_version=source_version, retrieval_date=retrieval_date, ...
      overwrite_family=kwargs.overwrite_family, overwrite=kwargs.overwrite, ...
      entry_callback=@caseEntry);
end

function s = emptyState()
   %EMPTYSTATE Prototype Laugh-Tests staging state.
   s = struct('case_id', "", 'entry', struct());
end

function s = stageLaughCase(case_id, source_dir, family_root, source_status, ...
      kwargs)
   %STAGELAUGHCASE Stage one Laugh-Tests case and return importer state.
   case_root = fullfile(family_root, case_id);
   evaluation_output_file = fullfile(case_root, "evaluation.mat");
   reference_output_file = fullfile(case_root, "reference.mat");

   % The manifest definition validates the requested case without reading raw
   % files; the build helper owns every source-specific read/normalization step.
   case_values = dryRunLaughCaseValues(case_id);
   if ~kwargs.dry_run
      % Real staging delegates the complete raw build after cache validation.
      assertLaughSourceAvailable(source_status);
      [~, targets, reference] = ...
         icemodel.verification.setup.buildLaughTestsArtifacts( ...
         source_dir, case_id);
      targets = icemodel.verification.setup.stampArtifactMetadata(targets);
      reference = icemodel.verification.setup.stampArtifactMetadata(reference);
      % Create or clear the output root only after source reads succeeded, so
      % skip-missing imports do not delete existing cases or leave empty roots.
      write_artifacts = icemodel.verification.setup.prepareCaseRoot( ...
         case_root, kwargs.overwrite, ["evaluation.mat", "reference.mat"]);
      % Keep each already-current half of the atomic source pair byte-stable;
      % a missing sibling may still be added by an ordinary repeated import.
      if write_artifacts(1)
         save(evaluation_output_file, 'targets');
      end
      if write_artifacts(2)
         save(reference_output_file, 'reference');
      end
   end

   s = struct('case_id', case_id, ...
      'entry', icemodel.verification.setup.makeCaseManifestEntry(case_values));
end

function assertLaughSourceAvailable(status)
   %ASSERTLAUGHSOURCEAVAILABLE Throw a stable per-case missing-source error.
   if isempty(status) || all([status.present])
      return
   end
   error('icemodel:verification:importLaughTests:missingSource', ...
      'Laugh-Tests source checkout is incomplete')
end

function entry = caseEntry(s)
   %CASEENTRY Return the staged Laugh-Tests case entry.
   entry = s.entry;
end

function case_values = dryRunLaughCaseValues(case_id)
   %DRYRUNLAUGHCASEVALUES Return source-light manifest values for a case.
   switch case_id
      case "colbeck1976"
         case_values = colbeckCaseValues();
      otherwise
         valid_cases = icemodel.verification.namelists.caseid("laugh_tests");
         error('unsupported Laugh-Tests verification case %s. Valid cases: %s', ...
            case_id, strjoin(valid_cases, ', '))
   end
end

function case_values = colbeckCaseValues()
   %COLBECKCASEVALUES Canonical Colbeck manifest fields.
   note = ['Single canonical Colbeck 1976 case. The staged evaluation.mat ' ...
      'carries two target sources keyed numerical_summa and ' ...
      'analytical_clark2017. The numerical_summa bundle is derived from the ' ...
      'frozen SUMMA validation outputs in KyleKlenk/Laugh-Tests. The ' ...
      'analytical_clark2017 bundle is the closed-form Clark 2017 wetting-' ...
      'front / kinematic-wave solution from ' ...
      'icemodel.verification.colbeck.analyticalSolution.'];

   observation_variables = icemodel.verification.setup.metadataStruct({ ...
      'snow_liquid_water_storage_m', ...
      'derived from mLayerVolFracLiq * mLayerDepth in snow layers'
      'bottom_outflow_mps', 'scalarRainPlusMelt'});
   case_values = { ...
      'colbeck1976'
      'synthetic_process'
      'colbeck1976'
      'Colbeck 1976 synthetic snow infiltration benchmark'
      ''
      {}
      ''
      fullfile('colbeck1976', 'evaluation.mat')
      fullfile('colbeck1976', 'reference.mat')
      '1 minute output / sub-hour forcing'
      struct('start', '1990-01-01 00:01:00', 'end', '1990-01-01 10:00:00')
      {'snow_liquid_water_storage_m', 'bottom_outflow_mps'}
      observation_variables
      note};
end
