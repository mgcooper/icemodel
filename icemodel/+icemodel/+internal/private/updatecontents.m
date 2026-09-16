function updatecontents(folder)
   %UPDATECONTENTS Write the Contents.m file of a folder and its subfolders.
   %
   %  updatecontents(FOLDER) writes FOLDER/Contents.m for the namespace
   %  folder FOLDER. The file lists every file of FOLDER and of each
   %  namespace ("+") or class ("@") folder below it, one row per file with
   %  the H1 line of that file. A group of rows starts with the folder name.
   %  `help` prints the file.
   %
   %  The listing skips hidden files and folders, files whose names contain
   %  "~" or ".mex", the Contents.m file itself, and every other subfolder,
   %  such as private/. A private function is not a namespace member, so it
   %  gets no row.
   %
   %  Do not edit Contents.m by hand or with the MATLAB Contents Report.
   %  Run icemodel.internal.makecontents instead.
   %
   %  This function adapts updateContents from the IOSR Matlab Toolbox
   %  (Copyright 2016 University of Surrey), through the matfunclib template
   %  tool (+tbx/+internal/makecontents.m and private/updatecontents.m) and
   %  the fixes of its groupstats copy.
   %
   % See also: makecontents

   folder = char(folder);
   [~, name] = fileparts(folder);

   % Row groups: the folder itself, then its namespace and class folders.
   dirs = [{folder}; listingFolders(folder)];
   rows = cell(numel(dirs), 1);
   h1 = cell(numel(dirs), 1);
   for d = 1:numel(dirs)
      [rows{d}, h1{d}] = folderRows(dirs{d}, folder);
   end
   rows = vertcat(rows{:});
   h1 = vertcat(h1{:});

   % Pad each file name to the longest name that carries an H1 line.
   width = max([0; cellfun(@length, rows(~cellfun(@isempty, h1)))]);

   % Compose every line through deblank, so the file carries no trailing
   % whitespace. An editor that strips trailing whitespace then leaves the
   % file unchanged.
   body = strings(numel(rows), 1);
   for row = 1:numel(rows)
      if isempty(h1{row})
         body(row) = deblank(['%   ' rows{row}]);
      else
         body(row) = deblank(['%   ' rows{row} ...
            repmat(' ', 1, width - length(rows{row})) ' - ' h1{row}]);
      end
   end
   header = ["% " + upper(name); "%"; ...
      "%   Contents file for " + upper(name) + " and its subfolders."];
   footer = ["%"; "%   updatecontents.m generated this file on " ...
      + string(datetime('now', 'Format', 'dd MMM yyyy')) + " at " ...
      + string(datetime('now', 'Format', 'HH:mm:ss')) + "."];
   writelines([header; body; footer], fullfile(folder, 'Contents.m'));
end

function folders = listingFolders(folder)
   %LISTINGFOLDERS Return the namespace and class folders below FOLDER.
   %
   % Every name between FOLDER and the subfolder must start with "+" or "@",
   % so the listing does not name the files of a private or data folder as
   % members. The paths are full and in case-insensitive order.

   listing = dir(fullfile(folder, '**'));
   listing = listing([listing.isdir] & ~startsWith({listing.name}, '.'));
   folders = fullfile({listing.folder}, {listing.name});
   relative = extractAfter(folders, length(folder) + 1);
   parts = cellfun(@(p) strsplit(p, filesep), relative, ...
      'UniformOutput', false);
   member = cellfun(@(p) all(startsWith(p, ["+", "@"])), parts);
   folders = unique(folders(member));
   [~, order] = sort(lower(folders));
   folders = reshape(folders(order), [], 1);
end

function [rows, h1] = folderRows(dirname, top)
   %FOLDERROWS Return the listing rows and H1 lines of one folder.
   %
   % A folder with files gives a blank row, a folder-name row, and one row
   % per file. A folder with no listed file gives no rows. An m-file row
   % carries the dotted namespace name, for example icemodel.column.solve.

   listing = dir(dirname);
   names = {listing(~[listing.isdir]).name}';
   names = names(~startsWith(names, '.') & ~contains(names, '~') ...
      & ~contains(names, '.mex') & ~strcmp(names, 'Contents.m'));
   rows = cell(0, 1);
   h1 = cell(0, 1);
   if isempty(names)
      return
   end
   [~, order] = sort(lower(names));
   names = names(order);

   % The namespace name is the trailing run of "+" and "@" folder names,
   % because MATLAB forms namespace names from those folders only.
   parts = strsplit(dirname, filesep);
   plain = find(~startsWith(parts, ["+", "@"]), 1, 'last');
   prefix = [strjoin(cellfun(@(p) p(2:end), parts(plain + 1:end), ...
      'UniformOutput', false), '.') '.'];

   lines = cell(numel(names), 1);
   for f = 1:numel(names)
      lines{f} = h1Line(fullfile(dirname, names{f}));
      [~, base, ext] = fileparts(names{f});
      if strcmpi(ext, '.m')
         names{f} = [prefix base];
      end
   end
   [~, topname] = fileparts(top);
   rows = [{''}; {upper([topname extractAfter(dirname, length(top))])}; names];
   h1 = [{''}; {''}; lines];
end

function line = h1Line(filename)
   %H1LINE Return the H1 line of an m-file, or '' for any other file.
   %
   % The H1 line is the first line of the help block, without its leading
   % percent signs, the leading function name, and a trailing period. Its
   % first letter is upper case. The help block starts on the first
   % non-blank line, or on the line after the function or classdef line and
   % its continuation lines. A file with code there has no help block, so a
   % later comment that documents the code gives no H1 line.

   line = '';
   [~, name, ext] = fileparts(filename);
   if ~strcmp(ext, '.m')
      return
   end
   text = strtrim(readlines(filename));
   text = text(text ~= "");

   % Step over the declaration, which can continue on more lines with "...".
   head = 1;
   if ~isempty(text) ...
         && ~isempty(regexp(text(1), '^(function|classdef)\>', 'once'))
      while head < numel(text) && endsWith(text(head), '...')
         head = head + 1;
      end
      head = head + 1;
   end
   if head > numel(text) || ~startsWith(text(head), '%')
      return
   end

   % Drop the percent signs and spaces, then the name the line starts with.
   line = char(regexprep(text(head), '^[%\s]+', ''));
   if startsWith(lower(line), lower(name))
      line = strtrim(line(length(name) + 1:end));
   end

   % Drop the period before the empty check, because a comment that holds
   % only the name and a period, such as %FOO., leaves no H1 text.
   if endsWith(line, '.')
      line = line(1:end - 1);
   end
   if isempty(line)
      return
   end
   line(1) = upper(line(1));
end
