function result = releaseMetadata(mode, kwargs)
   %RELEASEMETADATA Prepare, observe, or finalize release metadata.
   %
   %  result = icemodel.internal.releaseMetadata("prepare", ...
   %     version="1.2", date_released="2026-08-01")
   %  result = icemodel.internal.releaseMetadata("observe")
   %  result = icemodel.internal.releaseMetadata("finalize", ...
   %     version_doi="10.5281/zenodo.12345678")
   %
   % This maintainer-only command edits CITATION.cff and reads public GitHub
   % and Zenodo state. It never merges, tags, pushes, publishes, or creates a
   % pull request. Tests may replace the read-only services through the
   % function-handle options below.

   arguments
      mode (1, 1) string {mustBeMember(mode, ...
         ["prepare", "observe", "finalize"])}
      kwargs.version (1, 1) string = ""
      kwargs.date_released (1, 1) string = ...
         string(datetime("today", "Format", "yyyy-MM-dd"))
      kwargs.version_doi (1, 1) string = ""
      kwargs.cff_file (1, 1) string = ...
         string(icemodel.internal.fullpath("CITATION.cff"))
      kwargs.repository (1, 1) string = "mgcooper/icemodel"
      kwargs.timeout_seconds (1, 1) double ...
         {mustBeNonnegative, mustBeLessThanOrEqual( ...
         kwargs.timeout_seconds, 600)} = 600
      kwargs.poll_interval_seconds (1, 1) double ...
         {mustBePositive} = 15
      kwargs.fetcher (1, 1) function_handle = @fetchJson
      kwargs.waiter (1, 1) function_handle = @pause
      kwargs.clock (1, 1) function_handle = @monotonicSeconds
      kwargs.validator (1, 1) function_handle = @validateCff
      kwargs.preflight (1, 1) function_handle = @checkReleasePreconditions
   end

   % All modes require one readable CFF and a syntactically valid concept DOI.
   cff_file = kwargs.cff_file;
   if ~isfile(cff_file)
      error('icemodel:internal:releaseMetadata:cffMissing', ...
         'CFF file not found: %s', cff_file)
   end
   concept_doi = "10.5281/zenodo.11539329";

   % Prepare requires an explicit target. Read-only observation and DOI
   % finalization are always bound to the version persisted in CITATION.cff.
   persisted_version = string(icemodel.internal.readCffVersion(cff_file));
   if mode == "prepare"
      release_version = kwargs.version;
      if release_version == ""
         error('icemodel:internal:releaseMetadata:versionRequired', ...
            'prepare requires an explicit version')
      end
   else
      release_version = persisted_version;
      if kwargs.version ~= "" && kwargs.version ~= persisted_version
         error('icemodel:internal:releaseMetadata:cffVersionMismatch', ...
            'Requested version %s does not match CITATION.cff version %s', ...
            kwargs.version, persisted_version)
      end
   end
   validateReleaseVersion(release_version)

   % Dispatch only metadata-local behavior. Publication authority remains
   % outside this command.
   switch mode
      case "prepare"
         result = prepareRelease(cff_file, release_version, ...
            kwargs.date_released, kwargs.preflight, kwargs.validator);
      case "observe"
         result = observeRelease(release_version, kwargs.repository, ...
            concept_doi, kwargs.timeout_seconds, ...
            kwargs.poll_interval_seconds, kwargs.fetcher, kwargs.waiter, ...
            kwargs.clock);
      case "finalize"
         result = finalizeRelease(cff_file, release_version, ...
            kwargs.version_doi, concept_doi, kwargs.fetcher, ...
            kwargs.validator);
   end
end

function result = prepareRelease(cff_file, version, date_released, ...
      preflight, validator)
   %PREPARERELEASE Check local state and stage validated CFF metadata.

   validateReleaseDate(date_released)
   project_dir = string(fileparts(cff_file));
   if project_dir == ""
      project_dir = string(pwd);
   end
   ok = preflight(project_dir, version);
   if ~isequal(ok, true)
      error('icemodel:internal:releaseMetadata:preconditionFailed', ...
         'Release preflight did not report success')
   end

   % Update exactly the two top-level release scalars before validating the
   % staged file. Invalid candidate metadata never replaces the live CFF.
   text = fileread(cff_file);
   text = replaceTopLevelScalar(text, "version", version);
   text = replaceTopLevelScalar(text, "date-released", ...
      "'" + date_released + "'");
   installValidatedCff(cff_file, text, validator)
   resetVersionIfCanonical(cff_file)

   result = struct("mode", "prepare", "status", "ready", ...
      "version", version, "date_released", date_released, ...
      "cff_file", cff_file);
end

function result = observeRelease(version, repository, concept_doi, ...
      timeout_seconds, interval_seconds, fetcher, waiter, clock)
   %OBSERVERELEASE Poll public GitHub and Zenodo state within one hard bound.

   validateRepository(repository)
   github_url = "https://api.github.com/repos/" + repository ...
      + "/releases/tags/v" + version;
   zenodo_url = zenodoConceptUrl(concept_doi);
   started = clock();
   attempts = 0;
   github_seen = false;
   zenodo_seen = false;
   version_doi = "";

   % A valid but incomplete public state is resumable. Malformed or failed
   % public reads return an error result instead of masquerading as pending.
   while true
      elapsed = max(0, clock() - started);
      if elapsed >= timeout_seconds
         result = observationResult("pending", version, attempts, elapsed, ...
            github_seen, zenodo_seen, version_doi, "", "");
         return
      end
      attempts = attempts + 1;
      try
         remaining = timeout_seconds - elapsed;
         github_seen = githubReleaseSeen( ...
            fetcher, github_url, version, remaining);
         elapsed = max(0, clock() - started);
         if elapsed >= timeout_seconds
            result = observationResult("pending", version, attempts, ...
               elapsed, github_seen, zenodo_seen, version_doi, "", "");
            return
         end
         remaining = timeout_seconds - elapsed;
         [zenodo_seen, version_doi] = zenodoVersionSeen( ...
            fetcher(zenodo_url, remaining), version, concept_doi);
      catch exception
         result = observationResult("error", version, attempts, ...
            max(0, clock() - started), github_seen, zenodo_seen, ...
            version_doi, ...
            string(exception.identifier), string(exception.message));
         return
      end

      elapsed = max(0, clock() - started);
      if elapsed >= timeout_seconds
         result = observationResult("pending", version, attempts, elapsed, ...
            github_seen, zenodo_seen, version_doi, "", "");
         return
      end
      if github_seen && zenodo_seen
         result = observationResult("ready", version, attempts, elapsed, ...
            true, true, version_doi, "", "");
         return
      end

      % Never request a wait that crosses the caller's bounded deadline.
      waiter(min(interval_seconds, timeout_seconds - elapsed))
   end
end

function result = finalizeRelease(cff_file, version, version_doi, ...
      concept_doi, fetcher, validator)
   %FINALIZERELEASE Verify Zenodo lineage and add both CFF DOI identifiers.

   validateZenodoDoi(version_doi, "version_doi")
   if version_doi == concept_doi
      error('icemodel:internal:releaseMetadata:versionDoiIsConceptDoi', ...
         'version_doi must identify a version record, not the concept record')
   end

   % The public concept query is the authority for DOI, version, and concept
   % lineage. Do not write the CFF until all three agree.
   response = fetcher(zenodoConceptUrl(concept_doi), 30);
   verifyZenodoVersionDoi(response, version, version_doi, concept_doi)

   original = fileread(cff_file);
   [updated, added_concept, added_version] = ensureDoiIdentifiers( ...
      original, concept_doi, version_doi, version);
   if added_concept || added_version
      installValidatedCff(cff_file, updated, validator)
      resetVersionIfCanonical(cff_file)
   else
      % Finalization always validates, including a repeated idempotent call.
      validator(cff_file);
   end

   result = struct("mode", "finalize", "status", "finalized", ...
      "version", version, "concept_doi", concept_doi, ...
      "version_doi", version_doi, "added_concept", added_concept, ...
      "added_version", added_version, "cff_file", cff_file);
end

function seen = githubReleaseSeen(fetcher, url, version, timeout_seconds)
   %GITHUBRELEASESEEN Return true only for the requested published tag.

   try
      response = fetcher(url, timeout_seconds);
   catch exception
      if contains(string(exception.identifier), "404") ...
            || contains(string(exception.message), "404")
         seen = false;
         return
      end
      rethrow(exception)
   end

   if isstruct(response) && isfield(response, "not_found") ...
         && isequal(response.not_found, true)
      seen = false;
      return
   end
   if ~isstruct(response) || ~isscalar(response) ...
         || ~isfield(response, "tag_name")
      error('icemodel:internal:releaseMetadata:githubMalformed', ...
         'GitHub release response is malformed')
   end
   if string(response.tag_name) ~= "v" + version
      error('icemodel:internal:releaseMetadata:githubTagMismatch', ...
         'GitHub returned tag %s instead of v%s', ...
         string(response.tag_name), version)
   end

   % Drafts and prereleases are valid API responses but are not stable public
   % release success.
   seen = true;
   for field = ["draft", "prerelease"]
      if ~isfield(response, field)
         continue
      end
      value = response.(field);
      if ~isscalar(value) || ~(islogical(value) || isnumeric(value))
         error('icemodel:internal:releaseMetadata:githubMalformed', ...
            'GitHub %s state is malformed', field)
      end
      seen = seen && ~logical(value);
   end
end

function [seen, version_doi] = zenodoVersionSeen(response, version, concept_doi)
   %ZENODOVERSIONSEEN Find the requested version under the expected concept.

   records = zenodoRecords(response);
   seen = false;
   version_doi = "";
   for n = 1:numel(records)
      record_version = zenodoRecordVersion(records(n));
      if normalizeVersion(record_version) ~= version
         continue
      end
      verifyZenodoConcept(records(n), concept_doi)
      if ~isfield(records(n), "doi") || strlength(string(records(n).doi)) == 0
         error('icemodel:internal:releaseMetadata:zenodoMalformed', ...
            'Zenodo version record has no DOI')
      end
      validateZenodoDoi(string(records(n).doi), "record DOI")
      seen = true;
      version_doi = string(records(n).doi);
      return
   end
end

function verifyZenodoVersionDoi(response, version, version_doi, concept_doi)
   %VERIFYZENODOVERSIONDOI Prove one DOI's version and concept lineage.

   records = zenodoRecords(response);
   for n = 1:numel(records)
      if isfield(records(n), "doi") ...
            && string(records(n).doi) == version_doi
         if normalizeVersion(zenodoRecordVersion(records(n))) ~= version
            error('icemodel:internal:releaseMetadata:zenodoVersionMismatch', ...
               'Zenodo DOI %s does not report version %s', ...
               version_doi, version)
         end
         verifyZenodoConcept(records(n), concept_doi)
         return
      end
   end
   error('icemodel:internal:releaseMetadata:versionDoiNotFound', ...
      'Zenodo concept query did not expose version DOI %s', version_doi)
end

function records = zenodoRecords(response)
   %ZENODORECORDS Validate and return the public concept-query records.

   if ~isstruct(response) || ~isscalar(response) ...
         || ~isfield(response, "hits") || ~isstruct(response.hits) ...
         || ~isscalar(response.hits) || ~isfield(response.hits, "hits")
      error('icemodel:internal:releaseMetadata:zenodoMalformed', ...
         'Zenodo concept response is malformed')
   end
   records = response.hits.hits;
   if isempty(records)
      records = struct.empty();
   elseif ~isstruct(records)
      error('icemodel:internal:releaseMetadata:zenodoMalformed', ...
         'Zenodo hit list is malformed')
   end
end

function version = zenodoRecordVersion(record)
   %ZENODORECORDVERSION Read one required Zenodo metadata.version value.

   if ~isfield(record, "metadata") || ~isstruct(record.metadata) ...
         || ~isscalar(record.metadata) ...
         || ~isfield(record.metadata, "version")
      error('icemodel:internal:releaseMetadata:zenodoMalformed', ...
         'Zenodo record has no metadata.version')
   end
   version = string(record.metadata.version);
end

function verifyZenodoConcept(record, concept_doi)
   %VERIFYZENODOCONCEPT Require both DOI and record-id concept identities.

   concept_id = extractAfter(concept_doi, "10.5281/zenodo.");
   if ~isfield(record, "conceptdoi") ...
         || string(record.conceptdoi) ~= concept_doi ...
         || ~isfield(record, "conceptrecid") ...
         || string(record.conceptrecid) ~= concept_id
      error('icemodel:internal:releaseMetadata:zenodoConceptMismatch', ...
         'Zenodo record does not belong to concept DOI %s', concept_doi)
   end
end

function [text, added_concept, added_version] = ensureDoiIdentifiers( ...
      text, concept_doi, version_doi, version)
   %ENSUREDOIIDENTIFIERS Add missing concept and version CFF DOI entries.

   doi_tokens = regexp(text, ...
      '(?m)^[ \t]+value:[ \t]*[''"]?([^''"\s#]+)', 'tokens');
   doi_values = strings(numel(doi_tokens), 1);
   for n = 1:numel(doi_tokens)
      doi_values(n) = string(doi_tokens{n}{1});
   end
   added_concept = ~any(strcmpi(doi_values, concept_doi));
   added_version = ~any(strcmpi(doi_values, version_doi));
   if ~added_concept && ~added_version
      return
   end

   % Match the file's existing newline style and insert before the next
   % top-level key so the identifiers remain valid CFF sequence entries.
   eol = newline;
   if contains(text, sprintf('\r\n'))
      eol = sprintf('\r\n');
   end
   entries = "";
   if added_concept
      entries = entries + doiEntry( ...
         "The concept DOI representing all versions of IceModel", ...
         concept_doi, eol);
   end
   if added_version
      entries = entries + doiEntry( ...
         "The version DOI for IceModel v" + version, version_doi, eol);
   end

   identifier_starts = regexp(text, ...
      '(?m)^identifiers:[ \t]*(?:#.*)?\r?$', 'start');
   if numel(identifier_starts) > 1
      error('icemodel:internal:releaseMetadata:identifiersMalformed', ...
         'CFF contains more than one top-level identifiers key')
   end
   if isempty(identifier_starts)
      insertion = regexp(text, '(?m)^repository-code:', 'start', 'once');
      if isempty(insertion)
         error('icemodel:internal:releaseMetadata:identifiersMalformed', ...
            'CFF has neither identifiers nor repository-code')
      end
      entries = "identifiers:" + eol + entries;
   else
      top_level = regexp(text, ...
         '(?m)^[A-Za-z][A-Za-z0-9-]*:[^\r\n]*(?:\r?\n|$)', 'start');
      insertion = top_level(find(top_level > identifier_starts, 1));
      if isempty(insertion)
         insertion = strlength(string(text)) + 1;
         if ~endsWith(text, {newline, char(13)})
            entries = eol + entries;
         end
      end
   end

   % Character indexing keeps the untouched CFF bytes stable.
   entries = char(entries);
   text = [text(1:insertion - 1), entries, text(insertion:end)];
end

function entry = doiEntry(description, doi, eol)
   %DOIENTRY Render one DOI using the CFF identifiers sequence form.

   entry = "  - description: """ + description + """" + eol ...
      + "    type: doi" + eol ...
      + "    value: " + doi + eol;
end

function result = observationResult(status, version, attempts, elapsed, ...
      github_seen, zenodo_seen, version_doi, error_id, error_message)
   %OBSERVATIONRESULT Build the stable result contract for observe mode.

   result = struct("mode", "observe", "status", status, ...
      "version", version, "attempts", attempts, ...
      "elapsed_seconds", elapsed, ...
      "github_release_seen", github_seen, ...
      "zenodo_record_seen", zenodo_seen, "version_doi", version_doi, ...
      "error_identifier", error_id, "error_message", error_message);
end

function text = replaceTopLevelScalar(text, key, value)
   %REPLACETOPLEVELSCALAR Replace exactly one top-level CFF scalar line.

   pattern = "(?m)^" + key + "[ \t]*:[^\r\n]*";
   matches = regexp(text, pattern, 'match');
   if numel(matches) ~= 1
      error('icemodel:internal:releaseMetadata:cffScalarMalformed', ...
         'CFF must contain exactly one top-level %s field', key)
   end
   text = regexprep(text, pattern, key + ": " + value, 'once');
end

function installValidatedCff(cff_file, text, validator)
   %INSTALLVALIDATEDCFF Validate a same-directory stage before replacement.

   stage_folder = string(fileparts(cff_file));
   if stage_folder == ""
      stage_folder = string(pwd);
   end
   stage = string(tempname(stage_folder)) + ".cff";
   cleanup = onCleanup(@() deleteIfFile(stage));
   writeText(stage, text)
   validator(stage);
   [ok, message] = movefile(stage, cff_file, 'f');
   if ~ok
      error('icemodel:internal:releaseMetadata:cffInstallFailed', ...
         'Could not install validated CFF: %s', message)
   end
end

function writeText(filename, text)
   %WRITETEXT Write one UTF-8-compatible CFF text payload.

   fid = fopen(filename, 'w');
   if fid < 0
      error('icemodel:internal:releaseMetadata:cffWriteFailed', ...
         'Could not open staged CFF: %s', filename)
   end
   cleanup = onCleanup(@() fclose(fid));
   count = fwrite(fid, char(text), 'char');
   if count ~= strlength(string(text))
      error('icemodel:internal:releaseMetadata:cffWriteFailed', ...
         'Could not write the complete staged CFF: %s', filename)
   end
end

function deleteIfFile(filename)
   %DELETEIFFILE Remove only the command's unresolved temporary stage.

   if isfile(filename)
      delete(filename)
   end
end

function resetVersionIfCanonical(cff_file)
   %RESETVERSIONIFCANONICAL Refresh the process cache after a live CFF edit.

   canonical = string(icemodel.internal.fullpath("CITATION.cff"));
   if canonicalPath(cff_file) == canonicalPath(canonical)
      icemodel.internal.version('reset');
   end
end

function path = canonicalPath(filename)
   %CANONICALPATH Resolve relative components and symbolic links.

   file = javaObject('java.io.File', char(string(filename)));
   path = string(file.getCanonicalPath());
end

function ok = checkReleasePreconditions(project_dir, version)
   %CHECKRELEASEPRECONDITIONS Require clean dev state and an unused tag.

   branch = runGit(project_dir, "branch --show-current");
   if branch ~= "dev"
      error('icemodel:internal:releaseMetadata:notDevBranch', ...
         'Release preparation must start on dev, not %s', branch)
   end
   dirty = runGit(project_dir, "status --porcelain");
   if dirty ~= ""
      error('icemodel:internal:releaseMetadata:trackedChanges', ...
         'Release preparation requires a clean tracked worktree')
   end
   existing_tag = runGit(project_dir, "tag --list v" + version);
   remote_tag = runGit(project_dir, ...
      "ls-remote --tags origin refs/tags/v" + version);
   if existing_tag ~= "" || remote_tag ~= ""
      error('icemodel:internal:releaseMetadata:tagExists', ...
         'Release tag v%s already exists', version)
   end
   ok = true;
end

function output = runGit(project_dir, arguments)
   %RUNGIT Run one read-only Git query in the candidate project.

   command = "git -C " + shellQuote(project_dir) + " " + arguments;
   [status, text] = system(command);
   if status ~= 0
      error('icemodel:internal:releaseMetadata:gitFailed', ...
         'Git preflight failed: %s', strtrim(string(text)))
   end
   output = strtrim(string(text));
end

function validateCff(filename)
   %VALIDATECFF Run the repository's CFF schema validator.

   command = "uvx cffconvert --validate --infile " + shellQuote(filename);
   [status, text] = system(command);
   if status ~= 0
      error('icemodel:internal:releaseMetadata:cffInvalid', ...
         'CFF validation failed: %s', strtrim(string(text)))
   end
end

function response = fetchJson(url, timeout_seconds)
   %FETCHJSON Read one public JSON endpoint with a bounded request timeout.

   request_timeout = min(30, max(0.001, timeout_seconds));
   options = weboptions("ContentType", "json", "Timeout", request_timeout, ...
      "HeaderFields", ["User-Agent", "icemodel-release-maintainer"]);
   response = webread(char(url), options);
end

function seconds_now = monotonicSeconds()
   %MONOTONICSECONDS Return elapsed seconds from a process-local timer.

   persistent timer
   if isempty(timer)
      timer = tic;
   end
   seconds_now = toc(timer);
end

function url = zenodoConceptUrl(concept_doi)
   %ZENODOCONCEPTURL Build the fixed public concept-record query.

   concept_id = extractAfter(concept_doi, "10.5281/zenodo.");
   url = "https://zenodo.org/api/records?q=conceptrecid%3A" + concept_id ...
      + "&sort=mostrecent&size=10";
end

function version = normalizeVersion(version)
   %NORMALIZEVERSION Remove Zenodo's optional leading v from a version.

   version = string(version);
   if startsWith(lower(version), "v")
      version = extractAfter(version, 1);
   end
end

function validateRepository(repository)
   %VALIDATEREPOSITORY Require one URL-safe GitHub owner/name pair.

   if isempty(regexp(repository, ...
         '^[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+$', 'once'))
      error('icemodel:internal:releaseMetadata:repositoryInvalid', ...
         'repository must have GitHub owner/name form')
   end
end

function validateReleaseVersion(version)
   %VALIDATERELEASEVERSION Require a compact semantic release identifier.

   if isempty(regexp(version, ...
         '^\d+\.\d+(?:\.\d+)?(?:[-+][0-9A-Za-z.-]+)?$', 'once'))
      error('icemodel:internal:releaseMetadata:versionInvalid', ...
         'Release version is invalid: %s', version)
   end
end

function validateReleaseDate(date_released)
   %VALIDATERELEASEDATE Require one real ISO calendar date.

   if isempty(regexp(date_released, '^\d{4}-\d{2}-\d{2}$', 'once'))
      error('icemodel:internal:releaseMetadata:dateInvalid', ...
         'date_released must have YYYY-MM-DD form')
   end
   try
      parsed = datetime(date_released, "InputFormat", "yyyy-MM-dd", ...
         "Format", "yyyy-MM-dd");
   catch
      error('icemodel:internal:releaseMetadata:dateInvalid', ...
         'date_released is not a valid calendar date')
   end
   if string(parsed) ~= date_released
      error('icemodel:internal:releaseMetadata:dateInvalid', ...
         'date_released is not a valid calendar date')
   end
end

function validateZenodoDoi(doi, field_name)
   %VALIDATEZENODODOI Require the DOI shape minted by this Zenodo concept.

   if isempty(regexp(doi, '^10\.5281/zenodo\.\d+$', 'once'))
      error('icemodel:internal:releaseMetadata:doiInvalid', ...
         '%s is not a Zenodo DOI: %s', field_name, doi)
   end
end

function quoted = shellQuote(value)
   %SHELLQUOTE Protect one path from the host command shell.

   value = char(string(value));
   if ispc
      if contains(value, '%')
         error('icemodel:internal:releaseMetadata:unsafeWindowsPath', ...
            'Release command paths cannot contain %% on Windows')
      end
      quoted = string([char(34), value, char(34)]);
      return
   end
   embedded_quote = char([39, 34, 39, 34, 39]);
   value = strrep(value, char(39), embedded_quote);
   quoted = string([char(39), value, char(39)]);
end
