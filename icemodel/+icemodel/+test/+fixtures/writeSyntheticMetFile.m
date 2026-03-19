function filepath = writeSyntheticMetFile(metdir, simyear, varargin)
   %WRITESYNTHETICMETFILE Write a synthetic met file for tests.
   %
   %  filepath = icemodel.test.fixtures.writeSyntheticMetFile(metdir, simyear)

   % Package the caller overrides before forwarding to the canonical builder,
   % so file and in-memory synthetic met data stay defined in one place.
   args = [{'metdir', metdir}, varargin];
   [~, filepath] = icemodel.test.fixtures.makeSyntheticMetFile(simyear, ...
      args{:});
end
