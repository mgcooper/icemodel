function report = audit_formal_substep_failures(kwargs)
   %AUDIT_FORMAL_SUBSTEP_FAILURES Probe formal cases for dt-min/maxsubstep fallback.
   %
   %  report = audit_formal_substep_failures()
   %  report = audit_formal_substep_failures(tier="full", smbmodel="icemodel")
   %  report = audit_formal_substep_failures(output_file="/tmp/substep_audit.mat")
   %
   % Use this when you want to identify which formal cases emit the
   % check_substep dt-min/maxsubstep diagnostics without rebuilding baselines.

   arguments
      kwargs.tier (1, :) string ...
         {icemodel.validators.mustBeTestTierName(kwargs.tier)} = "full"
      kwargs.smbmodel (1, :) string ...
         {icemodel.validators.mustBeTestSmbmodelSelector(kwargs.smbmodel)} = "all"
      kwargs.solver {icemodel.validators.mustBeSolverFilter(kwargs.solver)} = []
      kwargs.simyear (1, 1) double {mustBeInteger, mustBePositive} = 2016
      kwargs.smoke_sites string = "kanm"
      kwargs.full_sites string = ["kanm"; "kanl"]
      kwargs.stop_on_first_issue (1, 1) logical = false
      kwargs.output_file (1, :) string = ""
   end

   % Bootstrap the canonical test config once for the whole audit.
   [~, ~, ~, ~, suite_cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment(); %#ok<ASGLU>

   % Resolve the formal case matrix up front so the audit runs one canonical
   % single-case workflow at a time.
   cases = icemodel.test.helpers.getRegressionCaseMatrix( ...
      tier=kwargs.tier, smbmodel=kwargs.smbmodel, solver=kwargs.solver, ...
      simyear=kwargs.simyear, smoke_sites=kwargs.smoke_sites, ...
      full_sites=kwargs.full_sites);
   if isempty(cases)
      error('no formal regression cases matched the requested filters')
   end

   % Probe each case independently and record whether checksubstep emitted
   % either the dt-min warning or the forced-advance maxsubstep fallback.
   rows = repmat(struct( ...
      'case_id', "", ...
      'smbmodel', "", ...
      'sitename', "", ...
      'solver', NaN, ...
      'simyears', [], ...
      'n_spinup_years', NaN, ...
      'elapsed_s', NaN, ...
      'has_dt_min_warning', false, ...
      'has_maxsubstep', false, ...
      'issue_excerpt', ""), height(cases), 1);

   for icase = 1:height(cases)
      c = cases(icase, :);
      fprintf('Substep audit case %d/%d: %s\n', ...
         icase, height(cases), c.case_id)

      opts = icemodel.test.helpers.setModelOptsForCase(c);

      t0 = tic;
      txt = evalc('[ice1, ice2] = icemodel.test.helpers.runSmbModel(opts);'); %#ok<NASGU,ASGLU>
      elapsed_s = toc(t0);

      has_dt_min_warning = contains(string(txt), '(dt_min)');
      has_maxsubstep = contains(string(txt), 'n_subfail == maxsubstep');

      rows(icase).case_id = string(c.case_id);
      rows(icase).smbmodel = string(c.smbmodel);
      rows(icase).sitename = string(c.sitename);
      rows(icase).solver = c.solver;
      rows(icase).simyears = opts.simyears;
      rows(icase).n_spinup_years = opts.n_spinup_years;
      rows(icase).elapsed_s = elapsed_s;
      rows(icase).has_dt_min_warning = has_dt_min_warning;
      rows(icase).has_maxsubstep = has_maxsubstep;
      rows(icase).issue_excerpt = extractIssueExcerpt(txt);

      if kwargs.stop_on_first_issue && (has_dt_min_warning || has_maxsubstep)
         rows = rows(1:icase);
         break
      end
   end

   % Return one compact report table plus the selector metadata that produced
   % it so the run can be reproduced later.
   report = struct();
   report.summary = struct2table(rows);
   report.tier = kwargs.tier;
   report.smbmodel = kwargs.smbmodel;
   report.solver = kwargs.solver;
   report.simyear = kwargs.simyear;
   report.smoke_sites = reshape(kwargs.smoke_sites, [], 1);
   report.full_sites = reshape(kwargs.full_sites, [], 1);
   report.timestamp_utc = datetime('now', 'TimeZone', 'UTC');
   report.matlab_version = string(version);
   report.host = string(computer);

   if ~isblanktext(kwargs.output_file)
      outdir = fileparts(char(kwargs.output_file));
      if ~isempty(outdir) && exist(outdir, 'dir') ~= 7
         mkdir(outdir);
      end
      save(char(kwargs.output_file), 'report');
   end

   disp(report.summary(:, {'case_id', 'elapsed_s', ...
      'has_dt_min_warning', 'has_maxsubstep'}))
end

function excerpt = extractIssueExcerpt(txt)
   %EXTRACTISSUEEXCERPT Keep only the checksubstep-related lines.

   lines = string(splitlines(string(txt)));
   keep = contains(lines, 'timestep = ') | ...
      contains(lines, 'n_subfail == maxsubstep');
   excerpt = strjoin(lines(keep), newline);
end
