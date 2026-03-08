function saveRunOpts(opts)
   %SAVERUNOPTS Save the resolved OPTS struct for a run.
   %
   %  icemodel.saveRunOpts(opts)

   optsfile = fullfile(opts.pathoutput, 'opts', ['opts_' opts.casename '.mat']);
   backupfile(optsfile, opts.backupflag);
   save(optsfile, 'opts');
end
