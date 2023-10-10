function [Tair,rh,wspd,Qsi,Qli,Pa,albedo,imonth,snowd] = SUBMET( ...
      met,itime,iter,dt_new,dt,maxiter)

   % I THINK THIS IS ONLY NEEDED IF WE WANT TO PASS IN MET DATA WITH AN HOURLY
   % OR LONGER TIMESTEP AND THEN INTERPOLATE TO SUB-HOURLY. IF WE USE 15
   % MINUTE DATA WE CAN ASSUME IT IS CONSTANT ACROSS THE TIMESTEP

   %   itime = a running calendar of model time                    [days]
   %   iter = met iteration
   %   dt = timestep of met forcing (1 hour, 15 min, etc.)      [seconds]
   %   dt_iter = time elapsed since the start of this met timestep   [seconds]
   %   dt_new = duration of the dynamic sub-step                    [seconds]

   %  Example:

   %   |---------o----------------|
   %  t1       itime              t2

   %  where t1 and t2 are met forcing timestamps and itime is model time

   %  dTair = Tair(t2) - Tair(t1);
   %  dt_iter = itime - t1;
   %  dt_met = t2-t1;
   %  Tair = Tair(t1) + dTair*dt_iter/dt;

   %  The met time step is one hour. The model time (itime) is usually at a
   %  timestep within the met step. dt_iter is that time minus the start. the
   %  function interpolates between the start of the met step and the start of
   %  the next met step

   Tair = met.tair(iter);
   rh = met.rh(iter);
   wspd = met.wspd(iter);
   Qsi = met.swd(iter);
   Qli = met.lwd(iter);
   Pa = met.psfc(iter);
   albedo = met.albedo(iter);
   snowd = met.snowd(iter);
   % xTsfc = T_old(1);
   % xTsfc = min(Tsfc,Tf);
   % xTsfc = bc_type*[min(Tsfc,Tf), T_old(1)];

   if itime == met.Time(iter) % || iter == maxiter
      %    itime = itime + seconds(dt_new);
      imonth = month(met.Time(iter));
      return
   else
      % interpolate the met data to the model time

      % met values at i and i+1
      imet = [ Tair rh wspd Qsi Qli albedo snowd Pa;                 ...
         met.tair(iter+1) met.rh(iter+1) met.wspd(iter+1)      ...
         met.swd(iter+1) met.lwd(iter+1) met.albedo(iter+1)    ...
         met.snowd(iter+1) met.psfc(iter+1)                    ];

      imet = imet(1,:) + diff(imet)*seconds(itime-met.Time(iter-1))/dt;
      % itime = itime + seconds(dt_new);
      imonth = month(met.Time(iter));

      % pull out the forcings
      Tair = imet(1);
      rh = max(imet(2),0.0);
      wspd = max(imet(3),0.0);
      Qsi = max(imet(4),0.0);    % wdir is 4 and precip is 5
      Qli = max(imet(5),0.0);
      albedo = max(imet(6),0.0);    % elevation is 9 and ? is 10
      snowd = max(imet(7),0.0);
      Pa = max(imet(8),0.0);


      % % previous method that did not use datetime/duration for itime
      %    dt_iter = itime-met.date(iter);
      %    imet = imet(1,:) + diff(imet)*dt_iter/dt;
      %    itime = itime + dt_new/86400.00;           % running model time in days
      %    imonth = month(met.Time(iter));
      %  imonth = month(datetime(itime,'convertfrom','datenum'));


      % seems like I could use met.Time(iter) + seconds(dt_sum) to find the
      % model time rather than itime, but maybe itime is needed to keep the
      % timing exact
   end
end

% due to numerical round-off, when model time hits a met timestep, the
% model time can be a nanosecond or so off from model time, so dt_iter will
% be on the order 1e-10, which has zero effect, and it would be expensive
% to round it every time, and using a higher level function like datetime
% will create trouble later when converting to julia

% % this is the version before switching to timetble:
% % interpolate the met data between timesteps
%     dt_iter = itime-met.date(iter);
%     imet = met(iter:iter+1,:);     % met values at i and i+1
%     imet = imet(1,:) + diff(imet)*dt_iter/dt_met;
%   % tmet = met(iter:iter+1,end);
%   % imet = interp1(tmet,imet,itime);
%     itime = itime + dt_sub/86400.00;           % running model time in days
%
%     % pull out the forcings
%     Tair = imet(1);
%     rh = max(imet(2),0.0);
%     wspd = max(imet(3),0.0);
%     Qsi = max(imet(6),0.0);    % wdir is 4 and precip is 5
%     Qli = max(imet(7),0.0);
%     albedo = max(imet(8),0.0);    % elevation is 9 and ? is 10
%   % xTsfc = T_old(1);
%   % xTsfc = min(Tsfc,Tf);
%   % xTsfc = bc_type*[min(Tsfc,Tf); T_old(1)];
