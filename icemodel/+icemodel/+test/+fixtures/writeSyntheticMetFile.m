function filepath = writeSyntheticMetFile(metdir, simyear, varargin)
%WRITESYNTHETICMETFILE Write a synthetic met file for tests.
%
%  filepath = icemodel.test.fixtures.writeSyntheticMetFile(metdir, simyear)

   args = [{'metdir', metdir}, varargin];
   [~, filepath] = icemodel.test.fixtures.makeSyntheticMetFile(simyear, ...
      args{:});
end
