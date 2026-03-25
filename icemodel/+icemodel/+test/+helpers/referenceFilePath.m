function pathname = referenceFilePath(name, kwargs)
   %REFERENCEFILEPATH Return the canonical reference file path.
   %
   %  pathname = icemodel.test.helpers.referenceFilePath("runoff")
   %  pathname = icemodel.test.helpers.referenceFilePath("runoff", simyear=2016)

   arguments
      name (1, :) string {mustBeMember(name, ["runoff"])} ...
         = "runoff"

      kwargs.simyear double ...
         = []
   end

   simyear = kwargs.simyear;
   testdir = icemodel.getpath('test');
   refdir = fullfile(testdir, 'references');

   switch name
      case "runoff"
         pathname = resolveRunoffPath(refdir, simyear);
   end
end

function pathname = resolveRunoffPath(refdir, simyear)
   %RESOLVERUNOFFPATH Resolve the default runoff reference file(s) to load.

   % Try the combined file first.
   pathname = fullfile(refdir, 'runoff_reference.mat');
   if exist(pathname, 'file') == 2
      return
   end

   % Try a year-specific file.
   if ~isempty(simyear)
      pathname = fullfile(refdir, ...
         sprintf('runoff_reference_%d.mat', simyear));
      if exist(pathname, 'file') == 2
         return
      end
   end

   % Fall back to any year-specific files.
   files = dir(fullfile(refdir, 'runoff_reference_*.mat'));
   if isempty(files)
      error('icemodel:test:noReference', ...
         'No runoff reference file found in %s', refdir)
   end
   pathname = fullfile(files(1).folder, files(1).name);
end
