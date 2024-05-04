function WRITEOUTPUT(ice1, ice2, opts, thisyear, time, swd, lwd, albedo)
   %WRITEOUTPUT Save output to disk.
   %
   %#codegen

   % Save the data
   if opts.saveflag

      % Post process.
      [ice1, ice2] = POSTPROC(ice1, ice2, opts, swd, lwd, albedo, time);

      % Set the output path for this year.
      filepath = fullfile(opts.pathoutput, int2str(opts.simyears(thisyear)));

      % Write the 1-d file.
      filename = ['ice1_' opts.casename '.mat'];
      filename = fullfile(filepath, filename);
      backupfile(filename, opts.backupflag);
      save(filename, 'ice1', '-v7.3');

      % Write the 2-d file.
      filename = ['ice2_' opts.casename '.mat'];
      filename = fullfile(filepath, filename);
      backupfile(filename, opts.backupflag);
      save(filename, 'ice2', '-v7.3');
   end

   % [ice1,ice2] = POSTPROC2(T_sfc,T_ice,frac_ice,frac_liq,df_liq, ...
   %    Time(isbetween(Time,t1,t2)),opts);
end
