function pathlist = getpath(pathtype, sitename, smbmodel, userdata, simyears, varargin)
   %GETPATH Return canonical icemodel data and run paths.
   %
   %  pathlist = icemodel.getpath('data')
   %  pathlist = icemodel.getpath('demo')
   %  pathlist = icemodel.getpath('test')
   %  pathlist = icemodel.getpath('input')
   %  pathlist = icemodel.getpath('eval')
   %  pathlist = icemodel.getpath('userdata')
   %  pathlist = icemodel.getpath('output', sitename, smbmodel)
   %  pathlist = icemodel.getpath('restart', sitename, smbmodel)
   %  pathlist = icemodel.getpath('output', sitename, smbmodel, userdata)
   %  pathlist = icemodel.getpath('output', sitename, smbmodel, userdata, ...
   %     simyears, testname)
   %  pathlist = icemodel.getpath('output', sitename, smbmodel, userdata, ...
   %     simyears, testname, startdate, enddate)
   %
   % This is the canonical path getter used by configureRun and by downstream
   % workflows that need stable path derivation without re-creating the
   % corresponding OPTS struct. DEMO returns the repo-local demo root and TEST
   % returns the repo-local test root; runtime data paths remain config-backed
   % through icemodel.config(...).

   narginchk(1, Inf)

   if nargin < 2, sitename = ''; end
   if nargin < 3, smbmodel = ''; end
   if nargin < 4, userdata = ''; end
   if nargin < 5 || isempty(simyears)
      simyears = '';
   else
      simyears = arrayfun(@num2str, simyears(:), 'un', false);
   end

   omitpart = @(value) isempty(value) || isblanktext(value);

   switch char(pathtype)
      case 'data'
         pathlist = appendParts(resolveDataPath(), varargin);

      case 'input'
         pathlist = appendParts(resolveInputPath(), varargin);

      case 'eval'
         pathlist = appendParts(resolveEvalPath(), varargin);

      case 'demo'
         pathlist = appendParts(resolveDemoPath(), varargin);

      case 'test'
         pathlist = appendParts(resolveTestPath(), varargin);

      case 'userdata'
         pathlist = appendParts(resolveUserdataPath(), varargin);

      case 'output'
         [window, varargin] = extractWindowTag(varargin);
         parts = [{resolveOutputPath(), sitename, smbmodel, userdata}, ...
            varargin, {window}];
         parts = parts(~cellfun(omitpart, parts));

         if isempty(simyears)
            pathlist = fullfile(parts{:});
         else
            pathlist = fullfile(parts{:}, simyears);
         end

      case 'restart'
         [window, varargin] = extractWindowTag(varargin);
         parts = [{resolveOutputPath(), sitename, smbmodel}, varargin, ...
            {window, 'restart'}];
         parts = parts(~cellfun(omitpart, parts));
         pathlist = fullfile(parts{:});

      otherwise
         error('unrecognized pathtype: %s', pathtype)
   end
end

function [tag, rest] = extractWindowTag(args)
   %EXTRACTWINDOWTAG Peel trailing datetime pair into a YYYYMMDD-YYYYMMDD tag.

   tag = '';
   rest = args;
   if numel(args) < 2
      return
   end
   tail2 = args(end-1:end);
   if ~all(cellfun(@(v) isa(v, 'datetime') && isscalar(v) && ~isnat(v), tail2))
      return
   end
   tag = char(sprintf('%s-%s', ...
      string(tail2{1}, 'yyyyMMdd'), string(tail2{2}, 'yyyyMMdd')));
   rest = args(1:end-2);
end

function pathlist = appendParts(pathname, parts)
   %APPENDPARTS Build one path from a base folder and optional trailing parts.

   parts = [{pathname}, parts];
   parts = parts(~cellfun(@(value) isempty(value) || isblanktext(value), parts));
   pathlist = fullfile(parts{:});
end

function pathname = resolveDataPath()
   %RESOLVEDATAPATH Return the configured or default data root.

   pathname = getenv('ICEMODEL_DATA_PATH');
   if isblanktext(pathname)
      pathname = icemodel.internal.fullpath('data');
   end
end

function pathname = resolveDemoPath()
   %RESOLVEDEMOPATH Return the canonical repo-local demo root.

   pathname = icemodel.internal.fullpath('demo');
end

function pathname = resolveInputPath()
   %RESOLVEINPUTPATH Return the configured or default input root.

   pathname = getenv('ICEMODEL_INPUT_PATH');
   if isblanktext(pathname)
      pathname = fullfile(resolveDataPath(), 'input');
   end
end

function pathname = resolveEvalPath()
   %RESOLVEEVALPATH Return the configured or default evaluation root.

   pathname = getenv('ICEMODEL_EVAL_PATH');
   if isblanktext(pathname)
      pathname = fullfile(resolveDataPath(), 'eval');
   end
end

function pathname = resolveTestPath()
   %RESOLVETESTPATH Return the repo-local formal test tree root.

   pathname = icemodel.internal.fullpath('test');
end

function pathname = resolveUserdataPath()
   %RESOLVEUSERDATAPATH Return the configured or default userdata root.

   pathname = getenv('ICEMODEL_USERDATA_PATH');
   if isblanktext(pathname)
      pathname = fullfile(resolveInputPath(), 'userdata');
   end
end

function pathname = resolveOutputPath()
   %RESOLVEOUTPUTPATH Return the configured or default output root.

   pathname = getenv('ICEMODEL_OUTPUT_PATH');
   if isblanktext(pathname)
      pathname = fullfile(resolveDataPath(), 'output');
   end
end
