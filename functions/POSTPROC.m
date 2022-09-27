function [ice1,ice2,met] = POSTPROC(ice1,ice2,met,opts)

% minimum required inputs:

   % ice1.Tsfc
   % ice2.Tice
   % ice2.f_liq
   % ice2.f_ice
   % ice2.Sc
   
% load physical constants   
   load('PHYSCONS','Tf','emissSB','cv_ice','cv_liq','ro_ice','ro_liq',  ...
                  'ro_air','k_liq','Ls','Lf','Rv','emiss');

% anything we want at the 15 min timestep can be added to diags
%    diags.Time     =  met.Time;
%    diags.Qc       =  enbal.Qc;
   
% Compute a mesh for plotting
   dz             =  opts.dz_thermal;
   Z              =  opts.z0_thermal;
   
% compute bulk density (kg/m3), heat capacity (J/kg/K), thermal K (W/m/K)
   T              =  ice2.Tice;
   f_liq          =  ice2.f_liq;
   f_ice          =  ice2.f_ice;
   [k_eff,k_vap]  =  GETGAMMA(T,f_liq,f_ice,ro_ice,k_liq,Ls,Rv,Tf);
   ro_sno         =  f_ice.*ro_ice+f_liq.*ro_liq+(1.0-f_liq-f_ice).*ro_air;
   cp_sno         =  (cv_ice.*f_ice+cv_liq.*f_liq)./ro_sno;
   
   ice2.k_eff     =  k_eff;                  % eff. thermal k
   ice2.k_vap     =  k_vap;                  % water vapor diffusivity
   ice2.cp_sno    =  cp_sno;                 % sp. heat cap.
   ice2.ro_sno    =  ro_sno;                 % ice density
   ice2.Qsub      =  ice2.Sc.*dz;            % subsurf absorbed shortwave
   ice2.Z         =  (dz/2:dz:Z-dz/2)';
   
% Compute the radiative heat fluxes
   ice1.Tsfc      = min(ice1.Tsfc,Tf);             % surface temp,cap at Tf
   ice1.albedo    = met.albedo;                    % albedo
   ice1.Qsi       = met.swd;                       % shortwave down
   ice1.Qsr       = ice1.Qsi.*ice1.albedo;         % shortwave up
   ice1.Qli       = emiss.*met.lwd;                % longwave down
   ice1.Qle       = emissSB.*ice1.Tsfc.^4;         % longwave up
   ice1.Qsn       = (1-ice1.albedo).*ice1.Qsi;     % net shortwave 
   ice1.Qln       = ice1.Qli-ice1.Qle;             % net longwave
   ice1.Qn        = ice1.Qsn+ice1.Qln;             % net radiation
   ice1.Qsip      = (1-ice1.chi).*ice1.Qsi;        % penetrated shortwave
   ice1.Qabs      = (1-ice1.albedo).*ice1.Qsi;     % absorbed shortwave
   ice1.Qsrf      = ice1.chi.*ice1.Qabs;           % skin shortwave
   ice1.Qsub      = (1-ice1.chi).*ice1.Qabs;       % subsurf shortwave
   ice1.Qbal      = ice1.Qsub+ice1.Qsrf-ice1.Qabs; % balance

% for reference, Qsub can also be computed this way
 % enbal.Qsub     = sum(ice2.Sc.*dz)';             % subsurf shortwave
 
% repeat for met
   met.swu        =  met.swd.*met.albedo;
   met.swn        =  met.swd.*(1-met.albedo);
   met.lwu        =  emissSB.*met.tsfc.^4;
   met.lwd        =  emiss.*met.lwd;
   met.lwn        =  met.lwd - met.lwu;
   met.netr       =  met.swn + met.lwn;

   % do this after computing lwu etc
   met.tsfc       =  met.tsfc-Tf;
   met.tair       =  met.tair-Tf;
   
% Convert temperature from Kelvin to Celsius
   ice1.Tsfc      =  ice1.Tsfc-Tf;             % surface temp,cap at Tf
   ice2.Tice      =  ice2.Tice-Tf;              % ice temperature
   
% calculate surface and subsurface runoff
   if opts.skinmodel == true
      ice1  =  SRFRUNOFF(ice1,ro_liq,Ls,Lf,opts.dt);
   else
      ice1  =  ICERUNOFF(ice1,ice2,opts);
   end
   
% % temporary - take Tflag out of enbal
%    Tflag          = enbal.Tflag;
%    enbal          = rmfield(enbal,'Tflag');
%    diags.Tflag    = Tflag;
   
% convert enbal and ice1 to timetables
   dates       = met.Time; dates.TimeZone = 'UTC';
   ice1        = struct2table(ice1);
   ice1        = table2timetable(ice1,'RowTimes',dates);
   ice1        = retime(ice1,'hourly','mean');
   met         = retime(met,'hourly','mean');   

% Retime ice2 to hourly and round the data
   fields      = fieldnames(ice2);
   oldtime     = dates;
   newtime     = ice1.Time;

   for mm = 1:numel(fields)

      thisfield = fields{mm};

      if strcmp(thisfield,'Z') || strcmp(thisfield,'LCflag')
         continue
      end
      dat               = transpose(ice2.(thisfield));
      newdat            = interp1(oldtime,dat,newtime);
      ice2.(thisfield)  = transpose(newdat);

      if ismember(thisfield,{'f_ice','f_liq','k_vap','k_eff'})
         ice2.(thisfield) = round(ice2.(thisfield),5);
      elseif ismember(thisfield,{'h_melt','h_freeze','Tice'})
         ice2.(thisfield) = round(ice2.(thisfield),3);
      elseif ismember(thisfield,{'cp_sno','ro_sno'})
         ice2.(thisfield) = round(ice2.(thisfield),1);
      elseif ismember(thisfield,{'Qsub','Sc','errT','errH','df_liq','df_drn'})
         ice2.(thisfield) = round(ice2.(thisfield),8);
      end
   end
   
   % repeat rounding for ice1. in this case, we round all fields to 5
   % digits, so i use 'fields' in the logical checks below, but if i had
   % vars i wanted rounded to different precision, would need to replace 
   fields = ice1.Properties.VariableNames;
   for mm = 1:numel(fields)
      thisfield = fields{mm};      
      if ismember(thisfield,fields)
         ice1.(thisfield) = round(ice1.(thisfield),5);
      end
   end
   
   % deal with retiming LCFlag to hourly
   if isfield(ice2,'LCflag')
      LCflag   = false(size(ice2.f_ice));
      idx      = 0;
      for mm = 1:4:size(ice2.LCflag,2)-3
         idx = idx+1;
         LCflag(:,idx) = sum(ice2.LCflag(:,mm:mm+3),2)>0;
      end
   
      ice1.LCflag = transpose(LCflag(1,:));
      ice2        = rmfield(ice2,'LCflag');
   end
   
   ice2.Time   = ice1.Time;

% rename enbal to match the naming conventions i use everywhere else
   oldvars  =  {'Qsi','Qsr','Qsn','Qli','Qle','Qln','Qh','Qe','Qc','Qn','Tsfc'};
   newvars  =  {'swd','swu','swn','lwd','lwu','lwn','shf','lhf','chf','netr','tsfc'};
   ice1    =  renamevars(ice1,oldvars,newvars);      
   
   
% The surface energy balance: Qm = chi*Qsi*(1-albedo) + Qln + Qh + Qe + Qc
% The subsurface shortwave balance: Qsub = (1-chi)*Qsi*(1-albedo) 
%
%                      Qsi  Qsr
%                       |    ^               
%  skin (wall)          v    |     Qabs = Qsi*(1-albedo) = Qsrf + Qsub
% ----------------------|---/---   Qsrf = chi*Qsi*(1-albedo) = chi*Qabs
%  surface layer        v  /    \
% -------------------- Qsip----- > Qsub = Qsi*(1-chi)*(1-albedo)
%  subsurface layers    |       /       = Qsip*(1-albedo) 
%                       v      /        = Qabs*(1-chi)  
% ------------------------------        
% 
% The total absorbed shortwave radiation equals the sum of the 'skin'
% absorbed radiation (Qsrf) and the subsurface absorbed radiation (Qsub)
% The subsurface absorbed radiation is the integral from z=0 to z=infty of
% dQ/dz*dz i.e. sum(Sc.*dz) with Sc the source term in units W/m3. chi is
% just the ratio of Qsrf to Qabs i.e. how much of the total absorbed energy
% is absorbed by the wall. Qsip is the penetrating radiation, some of which
% is absorbed and some of which is reflected (contributes to Qsr)
% in short, if we call the net solar downflux in the ice Q, then:
% Qsub = dQ is the absorbed solar downflux in the ice in each layer
% Sc = dQ/dz is the source term (divergence of the net solar downflux)
    
   
   
   
   
   
   
% Save the ice surface data
%    ice1.ro_sno1   = tocol(ro_ice.*ice2.f_ice(1,:));   % top cell density
%    ice1.Tice1     = tocol(ice2.Tice(1,:));            % top cell temp
%    ice1.dH_sum    = tocol(sum(ice2.dH));              % delta enthalpy
%    ice1.dH_avg    = tocol(mean(ice2.dH));
   
% Convert the layer combination and melt-freeze flags to logicals
%  ice1.ablation  = (ice2.h_melt(1,:)'+ice1.surf_subl).*ro_liq./ro_ice
%  ice1.LCflag    = logical(min(sum(ice1.LCflag),1));
%  ice1.MFflag    = logical(min(sum(ice1.MFflag),1));
%  ice1.LCflag    = sum(ice2.LCflag);    % layer combo flag

   % vars2   = fieldnames(ice2);

   % % remove error following layer combinations
   % iLC         = find(LCflag) + 1;
   % Terr(iLC,:) = nan;
   % Terr_max    = max(Terr);
   % figure; plot(Terr_max);



   
% % delete this is FLIQ2ROF is working correctly   
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % with d_liq set zero every timestep, this is correct:
%    d_rof       = cumsum( sum( dz.*ice2.df_liq ) );
%    d_drn       = cumsum( sum( dz.*ice2.df_drn ) );
%    
%  % if d_liq = 0.0 every step, rather than column of zeros
%  % d_drn       = cumsum( dz.*ice2.df_drn ) ;
%  
%    d_mlt       = nan(size(d_rof));
%    d_frz       = nan(size(d_rof));
%    for n = 1:numel(d_rof)
%        imlt        = ice2.df_liq(:,n)>0;
%        ifrz        = ice2.df_liq(:,n)<0;
%        d_mlt(n)    = sum(dz.*ice2.df_liq(imlt,n));
%        d_frz(n)    = sum(-dz.*ice2.df_liq(ifrz,n));
% 
%        % this mimics the h_freeze and h_melt that used to be in ice2
%        ice2.h_melt(imlt,n)     = dz.*ice2.df_liq(imlt,n);
%        ice2.h_freeze(ifrz,n)   = -dz.*ice2.df_liq(ifrz,n);
%    end
%    d_mlt       = cumsum(d_mlt);
%    d_frz       = cumsum(d_frz);
% 
%    ice1.depth_melt      =  d_mlt(:);  % cumulative melt
%    ice1.depth_freeze    =  d_frz(:);  % cumulative freeze
%    ice1.depth_drain     =  d_drn(:);  % cumulative drain
%    ice1.runoff          =  d_rof(:);
%    
%    % if we remove the cumsum here, we don't need the diff() in ICERUNOFF
%    ice2.h_melt          =  cumsum(ice2.h_melt,2);
%    ice2.h_freeze        =  cumsum(ice2.h_freeze,2);
%    
%    % note: cumsum(sum(ice2.h_melt)) = ice1.depth_melt
%    % figure;plot(sum(ice2.h_melt)); hold on; plot(ice1.depth_melt,':')
%    % but we want ice2.h_melt to know the layers that melt/freeze
   
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
   
   

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % if d_liq is set zero every timestep, we need cumsum, otherwise we dont
% % if d_liq is NOT set zero every timestep, we need diff(d_liq) to get melt/freeze
%    
%    d_rof       = sum( dz.*ice2.df_liq );
%    d_drn       = sum( dz.*ice2.df_drn );
%    df_liq      = [zeros(size(ice2.Z)) diff(ice2.df_liq,[],2)];
%    d_mlt       = nan(size(d_rof));
%    d_frz       = nan(size(d_rof));
%    for n = 1:numel(d_rof)
%        imlt        = df_liq(:,n)>0;
%        ifrz        = df_liq(:,n)<0;
%        d_mlt(n)    = sum(dz.*df_liq(imlt,n));
%        d_frz(n)    = sum(-dz.*df_liq(ifrz,n));
% 
%        % this mimics the h_freeze and h_melt that used to be in ice2
%        ice2.h_melt(imlt,n)     = dz.*df_liq(imlt,n);
%        ice2.h_freeze(ifrz,n)   = -dz.*df_liq(ifrz,n);
%    end
%    d_mlt       = cumsum(d_mlt);
%    d_frz       = cumsum(d_frz);
% 
%    ice1.depth_melt     =   d_mlt(:);  % cumulative melt
%    ice1.depth_freeze   =   d_frz(:);  % cumulative freeze
%    ice1.depth_drain    =   d_drn(:);  % cumulative drain
%    ice1.runoff         =   d_rof(:);
%    ice2.h_melt         =   cumsum(ice2.h_melt,2);
%    ice2.h_freeze       =   cumsum(ice2.h_freeze,2);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~