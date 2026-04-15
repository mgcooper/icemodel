function writeoutput(ice1, ice2, opts, thisyear, time, swd, lwd, albedo)
   %WRITEOUTPUT Post-process and save ice1/ice2 output to disk.
   %
   %  icemodel.writeoutput(ice1, ice2, opts, thisyear, time, swd, lwd, albedo)
   %
   % Post-processes the raw annual output structs via icemodel.postprocess and
   % saves the result as two .mat files under opts.pathoutput/simyear/:
   %   ice1_<casename>.mat  — 1-D (per-timestep scalar) outputs
   %   ice2_<casename>.mat  — 2-D (column profile × timestep) outputs
   %
   % A backup is made of any existing file when opts.backupflag is true.
   %
   %#codegen

   if ~opts.saveflag
      return
   end

   % Post-process raw annual output.
   [ice1, ice2] = icemodel.postprocess(ice1, ice2, opts, swd, lwd, albedo, time);

   % Resolve output directory for this simulation year.
   filepath = fullfile(opts.pathoutput, int2str(opts.simyears(thisyear)));

   % Write the 1-D scalar output file.
   filename = fullfile(filepath, ['ice1_' opts.casename '.mat']);
   backupfile(filename, opts.backupflag);
   save(filename, 'ice1', '-v7.3');

   % Write the 2-D column profile output file.
   filename = fullfile(filepath, ['ice2_' opts.casename '.mat']);
   backupfile(filename, opts.backupflag);
   save(filename, 'ice2', '-v7.3');
end
