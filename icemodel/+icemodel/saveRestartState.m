function saveRestartState(opts, simyear, T, f_ice, f_liq, Ts, varargin)
%SAVERESTARTSTATE Save the year-boundary state needed for a restart.
%
%  icemodel.saveRestartState(opts, simyear, T, f_ice, f_liq, Ts)
%  icemodel.saveRestartState(opts, simyear, T, f_ice, f_liq, Ts, r_eff)

   restart = struct();
   restart.simyear = simyear;
   restart.T = T;
   restart.f_ice = f_ice;
   restart.f_liq = f_liq;
   restart.Ts = Ts;
   if nargin > 6
      restart.r_eff = varargin{1};
   end
   restart.casename = string(opts.casename);
   restart.smbmodel = string(opts.smbmodel);
   restart.opts = opts;
   restart.timestamp_utc = datetime('now', 'TimeZone', 'UTC');

   filepath = icemodel.restartfile(opts, simyear);
   backupfile(filepath, opts.backupflag);
   save(filepath, 'restart');
end
