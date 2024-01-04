function WRITEOUTPUT(ice1, ice2, opts, thisyear, time, swd, lwd, albedo)
   %WRITEOUTPUT Save output to disk.
   %
   %#codegen

   % Save the data
   if opts.savedata

      % Post process
      [ice1, ice2] = POSTPROC(ice1, ice2, opts, swd, lwd, albedo, time);

      % Write the 1-d file
      fname = ['ice1_' opts.casename];
      save( ...
         fullfile(opts.pathoutput, int2str(opts.simyears(thisyear)), fname), ...
         'ice1', '-v7.3');

      % Write the 2-d file
      fname = ['ice2_' opts.casename];
      save( ...
         fullfile(opts.pathoutput, int2str(opts.simyears(thisyear)), fname), ...
         'ice2', '-v7.3');
   end

   % [ice1,ice2] = POSTPROC2(T_sfc,T_ice,frac_ice,frac_liq,df_liq, ...
   %    Time(isbetween(Time,t1,t2)),opts);
end
