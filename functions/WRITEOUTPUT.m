function WRITEOUTPUT(ice1, ice2, opts, thisyear, time, swd, lwd, albedo)

% save the data
if opts.savedata

   % post process
   [ice1, ice2] = POSTPROC(ice1, ice2, opts, swd, lwd, albedo, time);
                  
   % write the files
   save(fullfile(opts.pathoutput, int2str(opts.simyears(thisyear)), ...
      ['ice1_' opts.casename]), 'ice1', '-v7.3');
   
   save(fullfile(opts.pathoutput, int2str(opts.simyears(thisyear)), ...
      ['ice2_' opts.casename]), 'ice2', '-v7.3');
end

% [ice1,ice2] = POSTPROC2(T_sfc,T_ice,frac_ice,frac_liq,df_liq, ...
%    Time(isbetween(Time,t1,t2)),opts);