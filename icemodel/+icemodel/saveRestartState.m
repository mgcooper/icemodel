function saveRestartState(opts, simyear, T, f_ice, f_liq, Ts)
%SAVERESTARTSTATE Save the year-boundary state needed for a restart.
%
%  icemodel.saveRestartState(opts, simyear, T, f_ice, f_liq, Ts)

   restart = struct();
   restart.simyear = simyear;
   restart.T = T;
   restart.f_ice = f_ice;
   restart.f_liq = f_liq;
   restart.Ts = Ts;
   restart.casename = string(opts.casename);
   restart.smbmodel = string(opts.smbmodel);
   restart.opts = opts;
   restart.timestamp_utc = datetime('now', 'TimeZone', 'UTC');

   filepath = icemodel.restartfile(opts, simyear);
   backupfile(filepath, opts.backupflag);
   save(filepath, 'restart');
end
