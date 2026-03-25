function printFilePath(pathname, action)
   %PRINTFILEPATH Print a file path truncated to the test/ directory.
   %
   %  icemodel.test.helpers.printFilePath(pathname)

   if nargin < 2
      action = "load";
   end

   % Truncate the pathname to the test/ directory for compact viewing.
   testdir = icemodel.getpath('test');
   parentdir = fileparts(testdir);
   truncated = extractAfter(string(pathname), parentdir);
   if ismissing(truncated) || strlength(truncated) == 0
      truncated = string(pathname);
   elseif startsWith(truncated, filesep)
      truncated = extractAfter(truncated, 1);
   end

   % Print the file path to the console.
   switch action
      case "load"
         fprintf('  Loaded: %s\n', truncated);
      case "save"
         fprintf('  Saved: %s\n', truncated);
   end
end
