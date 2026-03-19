function pathlist = setpath(pathtype, sitename, smbmodel, userdata, simyears, varargin)
   %SETPATH Return canonical icemodel data and run paths.
   %
   %  pathlist = icemodel.setpath('data')
   %  pathlist = icemodel.setpath('input')
   %  pathlist = icemodel.setpath('eval')
   %  pathlist = icemodel.setpath('test')
   %  pathlist = icemodel.setpath('userdata')
   %  pathlist = icemodel.setpath('output', sitename, smbmodel)
   %  pathlist = icemodel.setpath('restart', sitename, smbmodel)
   %  pathlist = icemodel.setpath('output', sitename, smbmodel, userdata)
   %  pathlist = icemodel.setpath('output', sitename, smbmodel, userdata, ...
   %     simyears, testname)
   %
   % This is the canonical path builder used by configureRun and by
   % downstream workflows that need stable path derivation without re-creating
   % the corresponding OPTS struct.

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

      case 'test'
         pathlist = appendParts(icemodel.internal.fullpath('test'), varargin);

      case 'userdata'
         pathlist = appendParts(resolveUserdataPath(), varargin);

      case 'output'
         parts = [{resolveOutputPath(), sitename, smbmodel, userdata}, varargin];
         parts = parts(~cellfun(omitpart, parts));

         if isempty(simyears)
            pathlist = fullfile(parts{:});
         else
            pathlist = fullfile(parts{:}, simyears);
         end

      case 'restart'
         parts = [{resolveOutputPath(), sitename, smbmodel}, varargin, {'restart'}];
         parts = parts(~cellfun(omitpart, parts));
         pathlist = fullfile(parts{:});

      otherwise
         error('unrecognized pathtype: %s', pathtype)
   end
end

function pathlist = appendParts(pathname, parts)
   parts = [{pathname}, parts];
   parts = parts(~cellfun(@(value) isempty(value) || isblanktext(value), parts));
   pathlist = fullfile(parts{:});
end

function pathname = resolveDataPath()
   pathname = getenv('ICEMODEL_DATA_PATH');
   if isblanktext(pathname)
      pathname = icemodel.internal.fullpath('data');
   end
end

function pathname = resolveInputPath()
   pathname = getenv('ICEMODEL_INPUT_PATH');
   if isblanktext(pathname)
      pathname = fullfile(resolveDataPath(), 'input');
   end
end

function pathname = resolveEvalPath()
   pathname = getenv('ICEMODEL_EVAL_PATH');
   if isblanktext(pathname)
      pathname = fullfile(resolveDataPath(), 'eval');
   end
end

function pathname = resolveUserdataPath()
   pathname = getenv('ICEMODEL_USERDATA_PATH');
   if isblanktext(pathname)
      pathname = fullfile(resolveInputPath(), 'userdata');
   end
end

function pathname = resolveOutputPath()
   pathname = getenv('ICEMODEL_OUTPUT_PATH');
   if isblanktext(pathname)
      pathname = fullfile(resolveDataPath(), 'output');
   end
end
