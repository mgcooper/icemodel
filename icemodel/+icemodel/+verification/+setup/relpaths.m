function rel = relpaths(filenames, base)
   %RELPATHS Reduce absolute staged paths to base-relative names for JSON.
   %
   %  rel = icemodel.verification.setup.relpaths(filenames, base)
   %
   %  Shared by the firn importers (importPromiceSites, importSumup) to record
   %  staged met/userdata files in the manifest colocation record as
   %  base-relative names. Storing relative names keeps the manifest portable
   %  inside the staged data tree; listcases/loadColocatedData resolve them back
   %  to absolute paths at read time.
   %
   %  Inputs
   %    filenames : string/char/cellstr  absolute staged file paths.
   %    base      : string  the input subdir (met/ or userdata/) to strip.
   %
   %  Returns
   %    rel : 1xN string  the base-relative file names.
   %
   % See also: icemodel.verification.setup.importPromiceSites,
   %  icemodel.verification.setup.importSumup
   filenames = string(filenames);
   base = string(base);
   rel = erase(filenames, base + filesep);
   rel = reshape(rel, 1, []);
end
