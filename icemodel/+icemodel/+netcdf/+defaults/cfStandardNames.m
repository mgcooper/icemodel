function names = cfStandardNames(kwargs)
   %CFSTANDARDNAMES Load the official CF Standard Name Table as a set.
   %
   %  names = icemodel.netcdf.defaults.cfStandardNames()
   %  names = ... cfStandardNames(refresh=true)
   %  names = ... cfStandardNames(file="...cf-standard-name-table.xml")
   %
   % Returns a string array of every recognized CF standard name (the
   % <entry> ids plus the deprecated-but-recognized <alias> ids) from the
   % official CF Standard Name Table. This is the authoritative set the
   % canonical variable map (icemodel.netcdf.defaults.variables) validates
   % its standard_name fields against, so CF names are pulled from the table
   % PROGRAMMATICALLY rather than hand-maintained.
   %
   % The table XML is fetched once from cfconventions.org and CACHED in a
   % gitignored cache directory (.cache/cf-standard-names/). The fetch is
   % network-dependent; if the cache is absent AND the network is blocked,
   % the function falls back to the small committed fixture under
   % test/fixtures/cf-standard-names/ and warns. Refresh the cache with
   % refresh=true when a new CF table version is published.
   %
   % Name-value
   %  refresh : (false) force a fresh webread, overwriting the cache
   %  file    : explicit path to a table XML (bypasses fetch+cache)
   %
   % See also: icemodel.netcdf.defaults.variables,
   %  icemodel.netcdf.defaults.variable

   arguments
      kwargs.refresh (1, 1) logical = false
      kwargs.file (1, 1) string = ""
   end

   persistent CACHED
   if ~kwargs.refresh && strlength(kwargs.file) == 0 && ~isempty(CACHED)
      names = CACHED;
      return
   end

   if strlength(kwargs.file) > 0
      xmlfile = kwargs.file;
   else
      xmlfile = resolveTable(kwargs.refresh);
   end

   names = parseTable(xmlfile);

   if strlength(kwargs.file) == 0
      CACHED = names;
   end
end

%% Local functions
function xmlfile = resolveTable(refresh)
   %RESOLVETABLE Path to a usable CF table XML, fetching+caching if needed.

   url = ['https://cfconventions.org/Data/cf-standard-names/current/' ...
      'src/cf-standard-name-table.xml'];

   cachefile = fullfile(cacheDir(), 'cf-standard-name-table.xml');

   if ~refresh && isfile(cachefile)
      xmlfile = cachefile;
      return
   end

   % Fetch and cache. webread is used (not curl/wget, which the sandbox
   % blocks). On failure, fall back to the committed fixture.
   try
      opts = weboptions('Timeout', 30, 'ContentType', 'text');
      xmltext = webread(url, opts);
      if ~isfolder(cacheDir())
         mkdir(cacheDir())
      end
      fid = fopen(cachefile, 'w');
      assert(fid > 0, 'icemodel:netcdf:cfStandardNames:cacheWriteFailed', ...
         'cannot open CF table cache file for writing: %s', cachefile)
      closer = onCleanup(@() fclose(fid));
      fwrite(fid, xmltext);
      clear closer
      xmlfile = cachefile;
   catch err
      fixture = fullfile(fixtureDir(), 'cf-standard-name-table.xml');
      if isfile(fixture)
         warning('icemodel:netcdf:cfStandardNames:usingFixture', ...
            ['CF table fetch failed (%s); using committed fixture %s. ' ...
            'Refresh the cache with refresh=true when network is available.'], ...
            err.message, fixture)
         xmlfile = fixture;
      else
         rethrow(err)
      end
   end
end

function names = parseTable(xmlfile)
   %PARSETABLE Read <entry> and <alias> ids from a CF table XML.

   S = readstruct(xmlfile);

   entryids = string({S.entry.idAttribute});
   if isfield(S, 'alias') || isprop(S, 'alias')
      aliasids = string({S.alias.idAttribute});
   else
      aliasids = strings(1, 0);
   end

   names = unique([entryids(:); aliasids(:)]);
end

function d = cacheDir()
   %CACHEDIR Gitignored cache directory under the repo root.
   d = fullfile(repoRoot(), '.cache', 'cf-standard-names');
end

function d = fixtureDir()
   %FIXTUREDIR Committed CF-table fixture directory.
   d = fullfile(repoRoot(), 'test', 'fixtures', 'cf-standard-names');
end

function r = repoRoot()
   %REPOROOT Repository root inferred from this file's location.
   here = fileparts(mfilename('fullpath'));
   % .../icemodel/+icemodel/+netcdf/+defaults -> up 4 to repo root
   r = fileparts(fileparts(fileparts(fileparts(fileparts(here)))));
end
