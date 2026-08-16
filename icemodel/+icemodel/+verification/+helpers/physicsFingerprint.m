function fingerprint = physicsFingerprint(opts)
   %PHYSICSFINGERPRINT Fingerprint the model's default physics configuration.
   %
   %  fingerprint = icemodel.verification.helpers.physicsFingerprint()
   %  fingerprint = icemodel.verification.helpers.physicsFingerprint(opts)
   %
   % Returns a digest of the resolved default model option VALUES together
   % with the IceModel version from CITATION.cff. A cohort run stamps this
   % value. When a report is later built from that cohort, the report
   % recomputes it and compares.
   %
   % What it covers, and what it does not. The digest moves when a default
   % option value changes, and the version moves on a release bump. Neither
   % moves for a source-level physics change that touches no option and no
   % release. A matching fingerprint is therefore not proof that the saved
   % numbers are reproducible. Bead icemodel-j2w holds the decision on
   % whether to widen the scope to the solver source.
   %
   % The gate that consumes this fingerprint warns; it never errors. A
   % flag-off option change does not make saved results stale, and the owner
   % decides when a warning warrants a rerun.
   %
   % Fields
   %   opts_sha256      - SHA-256 of the default option values [string].
   %   icemodel_version - icemodel.internal.version() [string].
   %   excluded_fields  - Fields left out of the digest [string row], so a
   %                      later reader can tell what it does not cover.
   %
   % Every default option enters the digest unless EXCLUDEDOPTIONFIELDS names
   % it. Excluding by name, rather than listing the physics options, means a
   % new physics option joins the digest without a second edit here. The cost
   % is a warning when a non-physics default changes. That is the safe
   % direction for a warn-only gate.
   %
   % OPTS replaces the reference resolution below. The runner and the report
   % both call the no-argument form, so they compare one reference. The
   % argument exists so a test can show which options move the digest.
   %
   % See also: icemodel.setopts,
   %  icemodel.verification.report.buildAblationEvaluationReport,
   %  icemodel.verification.setup.textSha256, icemodel.internal.version

   arguments
      % One fixed reference case resolves the defaults. Keep this reference
      % fixed. setopts calls configureRun, which resolves several options
      % from the forcing: z_tair, z_wind, and z_relh come from the station's
      % instrument heights. The exclusion list drops those three, so changing
      % the reference station does not move the digest. A reference forcing
      % that resolves other options differently still can.
      %
      % setopts also touches the filesystem through configureRun, which
      % asserts that ICEMODEL_INPUT_PATH exists. Callers that must not fail
      % on a missing workspace have to guard this call.
      opts (1, 1) struct = icemodel.setopts("icemodel", "kanm", 2016, "kanm")
   end

   excluded = excludedOptionFields();
   for k = 1:numel(excluded)
      field = char(excluded(k));
      if isfield(opts, field)
         opts = rmfield(opts, field);
      end
   end

   % Sort the field order so the digest depends on the values alone.
   % jsonencode writes fields in struct order, so the sort must precede it.
   opts = orderfields(opts);
   fingerprint = struct( ...
      'opts_sha256', ...
      icemodel.verification.setup.textSha256(string(jsonencode(opts))), ...
      'icemodel_version', string(icemodel.internal.version()), ...
      'excluded_fields', excluded);
end

function fields = excludedOptionFields()
   %EXCLUDEDOPTIONFIELDS Name the options that are run identity, not physics.
   %
   % Four groups. The workspace paths and file names resolve against the
   % machine's IceModel configuration. Hashing them would make the same code
   % fingerprint differently on two computers. The case identity and the
   % output controls name which run this is and what it may write, not how
   % the model computes. The instrument heights come from the reference
   % station, so hashing them would tie the digest to that choice. The output
   % channel lists are the resolved expansion of output_profile. Appending a
   % diagnostic channel changes them and changes no computed value. Warning
   % about that is the false alarm this gate exists to remove.

   fields = [ ...
      "pathdata", "pathinput", "patheval", "pathuserdata", ...
      "pathoutput", "pathrestart", ...
      "metfname", "userdatafname", "initfile", "restartfile", ...
      "debug_path", "readiness_file", "report_inputs_file", ...
      "sitename", "forcings", "userdata", "uservars", ...
      "simyears", "numyears", "output_years", "startdate", "enddate", ...
      "casename", "testname", ...
      "saveflag", "saveopts", "backupflag", "saverestart", ...
      "output_profile", "vars1", "vars2", ...
      "z_tair", "z_wind", "z_relh"];
end
