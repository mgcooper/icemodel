%--------------------------------------------------------------------------
function h = plotRunoff(Runoff,Discharge,Catchment,varargin)
         
%--------------------------------------------------------------------------
%  need to enforce racmo, mar, merra ordering

   p               = inputParser;
   p.FunctionName  = 'plotRunoff';
   p.CaseSensitive = false;
   p.KeepUnmatched = true;

   addRequired(   p, 'Runoff',                     @(x)istimetable(x)   );
   addRequired(   p, 'Discharge',                  @(x)istimetable(x)   );
   addRequired(   p, 'Catchment',                  @(x)isstruct(x)      );
   addParameter(  p, 'plotcolumn',     true,       @(x)islogical(x)     );
   addParameter(  p, 'plotpatch',      false,      @(x)islogical(x)     );
   addParameter(  p, 'plotsurf',       false,      @(x)islogical(x)     );
   addParameter(  p, 'plotmelt',       false,      @(x)islogical(x)     );
   addParameter(  p, 'plotcomp',       false,      @(x)islogical(x)     );
   addParameter(  p, 'plotmean',       false,      @(x)islogical(x)     );
   addParameter(  p, 'plotmembers',    false,      @(x)islogical(x)     );
   addParameter(  p, 'plotensemble',   false,      @(x)islogical(x)     );
   addParameter(  p, 'refstart',       false,      @(x)islogical(x)     );
   addParameter(  p, 'legendtext',     {''},       @(x)iscell(x)        );
   addParameter(  p, 'sitename',       '',         @(x)ischar(x)        );
   addParameter(  p, 't1',             NaT,        @(x)isdatetime(x)    );
   addParameter(  p, 't2',             NaT,        @(x)isdatetime(x)    );
   
   
   parse(p,Runoff,Discharge,Catchment,varargin{:});
   
   opts.plot_patch      =  p.Results.plotpatch;
   opts.plot_column     =  p.Results.plotcolumn;
   opts.plot_surf       =  p.Results.plotsurf;
   opts.plot_melt       =  p.Results.plotmelt;
   opts.plot_comp       =  p.Results.plotcomp;
   opts.plot_mean       =  p.Results.plotmean;
   opts.plot_members    =  p.Results.plotmembers;
   opts.plot_ensemble   =  p.Results.plotensemble;
   opts.refstart        =  p.Results.refstart;
   opts.ltext           =  p.Results.legendtext;
   opts.sitename        =  p.Results.sitename;
   opts.t1              =  p.Results.t1;
   opts.t2              =  p.Results.t2;
   
   
   % plot_comp means compare two version of icemodel 
   % plot_
   
%--------------------------------------------------------------------------
   warning off
%    narginchk(3,6)
% 
%    numargin    = nargin;
%    argsin      = varargin;
%    opts        = set_opts(numargin,argsin);

   if opts.plot_surf == true
      opts.plot_column = false;
   end

   if opts.plot_surf == true
       Runoff.icemodel =   Runoff.skinmodelRunoff;
   elseif opts.plot_melt == true
       Runoff.icemodel =   Runoff.icemodelMelt;
   elseif opts.plot_column == true
       Runoff.icemodel =   Runoff.icemodelRunoff;
   end

%--------------------------------------------------------------------------

sitename    =  Catchment.sitename;

% % % % %  prep the runoff for plotting
[Q,R,txt]   =  prep_runoff(Runoff,Discharge,Catchment,opts);
                                                        % % % % % %

% % % % %  make the plot
h       =   make_plot(Q,R,txt,opts);
                                                        % % % % % % 

% this sends back the adcp data and the model data in the same format
data  = [Q.obs R.med'];
data  = array2table(data,'VariableNames',txt);
data  = table2timetable(data,'RowTimes',datetime(Q.t,'ConvertFrom','datenum'));
    
% return the data
h.data      = data;
% h.errbar    = herr;
% h.lines     = l;

% custom processing for slv plots
   if any(strcmp(sitename,{'slv1','slv2'}))
      h.h(1).Marker           = 'o';
      h.h(1).MarkerSize       = 10;
      h.h(1).Color            = [0.635,0.078,0.184];
      h.h(1).MarkerFaceColor  = [0.635,0.078,0.184];
      h.h(1).MarkerEdgeColor  = 'none';
      h.h(1).LineWidth        = 2;
   end

axis tight   
figformat(  'linelinewidth',2,'axeslinewidth',1,'labelinterpreter',     ...
            'tex','axesinterpreter','tex','textinterpreter','tex',   ...
            'legendinterpreter','tex','legendlocation','northwest',    ...
            'suppliedaxis',h.ax,'suppliedfigure',h.f);
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function [Q,R,txt] = prep_runoff(Runoff,Discharge,Catchment,opts)
    
% check timezone consistency
   if ~isequal(Runoff.Time.TimeZone,Discharge.Time.TimeZone)
      if isempty(Discharge.Time.TimeZone)
         Discharge.Time.TimeZone = Runoff.Time.TimeZone;
      end
   end
   
   sitename =  Catchment.sitename;
   
   % put the runoff into an array
   racmo    =  Runoff.racmoRunoff;
   mar      =  Runoff.marRunoff;
   merra    =  Runoff.merraRunoff;
   icemod   =  Runoff.icemodel;
   
   % check if icemodel is an ensemble
   if width(icemod)>1
      % we need the colors 
   end
   
   % the lakes require unique processing
   if any(strcmp(sitename,{'slv1','slv2'}))
      if opts.refstart == true
         t1    =  Discharge.Time(find(~isnan(Discharge.lake_area_km2),1,'first'));
      else
         t1    =  datetime(2015,6,25,'TimeZone','UTC');
      end
      t2       =  datetime(2015,8,5,'TimeZone','UTC');
      [si,ei]  =  dateInds(t1,t2,Runoff.Time);
      Y        =  [racmo(si:ei) mar(si:ei) merra(si:ei) icemod(si:ei,:)];
      RR       =  retime(Discharge,Runoff.Time,'fillwithmissing');
      RL       =  RR.lake_volume_min_km3(si:ei).*1e9;
      RH       =  RR.lake_volume_max_km3(si:ei).*1e9;
      Q.obs    =  (RL+RH)./2;
      Q.err    =  (RH-RL)./2;
      Q.t      =  datenum(RR.Time(si:ei));
      C1       =  Catchment.min.ease.area;
      C2       =  Catchment.max.ease.area;
      subQ0    =  false;
      ltxt     = 'SLV';
      
      Catchment.med.ease.area = (C1+C2)/2;
      
      if opts.refstart == true
         Q.obs = Q.obs-Q.obs(1,:);
      end
      
   else
      
      % get the dates to extract the rcm data during the observation period
      if strcmp(sitename,'ak4')
         
         yyyy     =  year(Runoff.Time(1));
         
         if isnat(opts.t1)
            t1    =  datetime(yyyy,6,1,'TimeZone','UTC');
         else
            t1    =  opts.t1;
         end
         
         if isnat(opts.t2)
            t2    =  datetime(yyyy,9,1,'TimeZone','UTC');
         else
            t2    =  opts.t2;
         end
         
         si       =  find(Discharge.Time == t1);
         ei       =  find(Discharge.Time == t2);
         
         % convert the observed discharge to cumulative runoff in m3
         Q.obs    =   cumsum(Discharge.Qm3(si:ei)); Q.obs = Q.obs-Q.obs(1);
         Q.err    =   cumsum((Discharge.Qm3H(si:ei)-Discharge.Qm3L(si:ei))./2);
         Q.t      =   datenum(Discharge.Time(si:ei));
         
         % reset si,ei to get the rcm data
         si       =  find(Runoff.Time == t1);
         ei       =  find(Runoff.Time == t2);
      else
         
         % convert the observed discharge to cumulative runoff in m3
         Q.obs    =   cumsum(Discharge.Qm3); Q.obs = Q.obs-Q.obs(1);
         Q.err    =   cumsum((Discharge.Qm3H-Discharge.Qm3L)./2);
         Q.t      =   datenum(Discharge.Time);
         
         si       =  find(Runoff.Time == Discharge.Time(1));
         ei       =  find(Runoff.Time == Discharge.Time(end));
         
      end
      
      Y        =  [racmo(si:ei) mar(si:ei) merra(si:ei) icemod(si:ei,:)];
      subQ0    =  true;
      ltxt     = 'ADCP';
      

   end

% convert runoff to km3 using the catchment areas (km2)
% put the data in rows so patch is easier to deal with
   for n = 1:size(Y,2)
      if subQ0 == true
         R.min(n,:) = Catchment.min.ease.area.*(cumsum(Y(:,n))-Y(1,n));
         R.med(n,:) = Catchment.med.ease.area.*(cumsum(Y(:,n))-Y(1,n));
         R.max(n,:) = Catchment.max.ease.area.*(cumsum(Y(:,n))-Y(1,n));
      else
         R.min(n,:) = Catchment.min.ease.area.*cumsum(Y(:,n));
         R.med(n,:) = Catchment.med.ease.area.*cumsum(Y(:,n));
         R.max(n,:) = Catchment.max.ease.area.*cumsum(Y(:,n));
      end
   end
      
% next makes legend text
   if ~isempty(opts.ltext)

      txt   = opts.ltext;

   else

      if opts.plot_ensemble == true
      %Y   =   [racmo.runoff(si:ei) mar.runoff(si:ei) icemod];
         txt  =   {ltxt,'RCMs','IceModel'};  % set the legend text

         if opts.plot_surf == true
            txt   =   {ltxt,'RCMs','SkinModel'};
         end 

      else

         if opts.plot_comp == true
            if opts.plot_surf == true
               txt   =   {ltxt,'RACMO','MAR','MERRA','SkinModel v2','SkinModel'};
            else
               txt   =   {ltxt,'RACMO','MAR','MERRA','IceModel v2','IceModel'};
            end
         else
            if opts.plot_surf == true

               txt   =   {ltxt,'RACMO','MAR','MERRA','SkinModel'};

            else
               txt   =   {ltxt,'RACMO','MAR','MERRA','IceModel'};
            end
         end   
      end
   end

end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


function h = make_plot(Q,R,ltext,opts)

sitename =  opts.sitename;

% Note: R goes: racmo, mar, merra, icemodel, skinmodel
% load('distinguishablecolors.mat','dc'); c = dc;   
N        =  size(R.med,1); 
c        =  distinguishable_colors(N+2);
tpatch   =  [Q.t' fliplr(Q.t')];

% ADCP (pf = plot frequency, lw = line width, ms = marker size)
pf = 1; if ismember(sitename,{'upperBasin','ak4','kanl'}); pf = 24; end % daily
lw = 1; if ismember(sitename,{'upperBasin'}); lw = 1; end
ms = 10;if ismember(sitename,{'upperBasin'}); lw = 10; end


c1       =   [0.635,0.078,0.184];
h.f      =   macfig('monitor','mac','size','full');
herr     =   myerrorbar(c1,Q.t(1:pf:end),Q.obs(1:pf:end),Q.err(1:pf:end),    ...
               'LineWidth',lw,'MarkerSize',ms); hold on; 
h.ax     =  gca;

% plot_mean    =  patch each RCM min/max, lines on top
% plot_members =  patch average of all RCM's, lines on top, + separate 
%                 patch for icemodel
% plot_

% should have an option to plot as patches or lines


% plot_patch means the data passed in must have 

if opts.plot_patch == true

   % patch each rcm, with the average
   if opts.plot_mean == true

       for n = 1:N
           ypatch  =   [R.max(n,:) fliplr(R.min(n,:))];
           p(n)    =   patch('XData',tpatch,'YData',ypatch,'FaceColor',    ...
                       c(n,:),'FaceAlpha',0.15,'EdgeColor','none'); 
           l(n)    =   plot(Q.t,R.med(n,:),'Color',c(n,:));
       end
       lh          =   [herr p];

   % plot one patch for mean min/max then the members on top    
   elseif opts.plot_members == true

       ypatch      =   [max(R.max(1:N-1,:)) fliplr(min(R.min(1:N-1,:)))];
       p1          =   patch('XData',tpatch,'YData',ypatch,'FaceColor',    ...
                       c(1,:),'FaceAlpha',0.15,'EdgeColor','none'); 
       for n = 1:N-1
           l(n)    =   plot(Q.t,R.med(n,:),'Color',[c(n,:) 0.5]);
       end

   % separate patch for icemodel
       ypatch  =   [R.max(N,:) fliplr(R.min(N,:))];
       p2      =   patch('XData',tpatch,'YData',ypatch,'FaceColor',  ...
                   c(N,:),'FaceAlpha',0.15,'EdgeColor','none');
       l(N)    =   plot(Q.t,R.med(N,:),'Color',c(N,:));
       lh      =   [herr l];

   elseif opts.plot_ensemble == true

      % patch icemodel
      yL       =   mean(R.min(N-1:N,:));
      yH       =   mean(R.max(N-1:N,:));
      ypatch   =   [yH fliplr(yL)];
      p2       =   patch('XData',tpatch,'YData',ypatch,'FaceColor',c(1,:), ...
                    'FaceAlpha',0.15,'EdgeColor','none');

      % lines for each model
      for n = 1:(N-2)/2
        l(n)   =   plot(Q.t,R.med(n,:),'Color',c(n,:));
      end

      % lines for each emulator
      for n = (N-2)/2+1:(N-2)
        l(n)   =   plot(Q.t,R.med(n,:),'Color',c(n-(N-2)/2,:),'LineStyle',':');
      end

      lh    = [herr l p2]; % ltext = ''; % not sure why this was ''
   end
   
% plot lines, no patches   
else

   if opts.plot_ensemble == true
      
      % note: these assume the last two entries are AWS and MODIS
      
      % lines for each model
      for n = 1:(N-2)/2
        l(n)   =   plot(Q.t,R.med(n,:),'Color',c(n,:));
      end

      % lines for each emulator
      for n = (N-2)/2+1:(N-2)
        l(n)   =   plot(Q.t,R.med(n,:),'Color',c(n-(N-2)/2,:),'LineStyle',':');
      end

      % lines for the icemodel AWS/MODIS forcing
      caws     =  c(N+1,:); % rgb('dark grey'); % c(N-1,:);
      cmod     =  c(N+2,:); % rgb('dark grey'); % c(N-1,:);
      l(N-1)   =  plot(Q.t,R.med(N-1,:),'Color',caws,'LineStyle','--');
      l(N)     =  plot(Q.t,R.med(N,:),'Color',cmod,'LineStyle','--');
      
      lh    = [herr l]; % ltext = ''; % not sure why this was ''
      
   else
      
      for n = 1:N
         l(n)  =   plot(Q.t,R.med(n,:),'Color',c(n,:));
      end
         lh    =   [herr l];
   end
end
   
legend(lh,ltext,'location','northwest')
ylabel('Cumulative Runoff [m^3]');
% ylabel('Cumulative Runoff [m$^3$]');


% if any(ismember(datevec(Q.t),2015))
%     h.ax.XLim   =   [Q.t(1)-6/24 Q.t(end)+12/24];
%     h.ax.YLim   =   [-0.2e6 6e6];
% elseif any(ismember(datevec(Q.t),2016))
%     h.ax.XLim   =   [Q.t(1)-13/24 Q.t(end)+13/24];
%     h.ax.YLim   =   [-0.05e7 2e7];
% end

h.h      = lh;       % lh has all handles
h.ltext  = ltext;

% h.ax.XTickLabel       =   datestr(h.ax.XTick,'mmm-DD');

datetick;


end

% % this clips SEB to the obsevation time limits
% if plot_comp == false
%     seb         =   seb.adcp;
%     [si,ei]     =   dateInds(Q.Time(1),Q.Time(end),seb.time);
%     seb_runoff  =   seb.runoff(si:ei);
% end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


% function opts = set_opts(numargin,argsin)
%     
%     opts.plot_alpha         =   0.25;
%     opts.plot_mean          =   true;
%     opts.plot_members       =   false;
%     opts.plot_raw           =   false;
%     opts.plot_sensitivity   =   false;
%     opts.plot_surf          =   false;
%     opts.plot_column        =   false;
%     opts.plot_comp          =   false;
%     opts.plot_melt          =   false;
%     
%     if numargin >= 4 && strcmp(argsin{1},'members')
%         opts.plot_mean      =   false;
%         opts.plot_members   =   true;
%     elseif numargin >= 4 && strcmp(argsin{1},'raw')
%         opts.plot_mean      =   false;
%         opts.plot_raw       =   true;
%     elseif numargin >= 4 && strcmp(argsin{1},'comp')
%         opts.plot_mean      =   false;
%         opts.plot_raw       =   true;
%         opts.plot_comp      =   true;
%     elseif numargin >= 4 && strcmp(argsin{1},'sensitivity')
%         opts.plot_mean      =   false;
%         opts.plot_raw       =   true;    
%         opts.plot_sensitivity =  true;
%     end
% 
%     if numargin == 5 && strcmp(argsin{2},'surf')
%         opts.plot_surf      =   true;
%     elseif numargin == 5 && strcmp(argsin{2},'melt')
%         opts.plot_melt      =   true;
%     elseif numargin == 5 && strcmp(argsin{2},'column')
%         opts.plot_column    =   true;
%     end
%     
% end