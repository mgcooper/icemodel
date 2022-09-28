function h = plotPromice(Ablation,varargin)
   
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
   p               = inputParser;
   p.FunctionName  = 'plotPromice';
   p.CaseSensitive = false;
   p.KeepUnmatched = true;

   addRequired(   p, 'Ablation',                   @(x)istimetable(x)   );
   addParameter(  p, 'refstart',       NaT,        @(x)isdatetime(x)    );
   addParameter(  p, 'sitename',       '',         @(x)ischar(x)        );
   
   
   parse(p,Ablation,varargin{:});
   
   refstart    =  p.Results.refstart;
   sitename    =  p.Results.sitename;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
   
    
    load('defaultcolors.mat','dc');
    
    if isempty(Ablation.Time.TimeZone)
        Ablation.Time.TimeZone = 'UTC';
    end

    Time    = Ablation.Time;
    yyyy    = year(Ablation.Time(1));
    
    if isnat(refstart)
      t1    = datetime(yyyy,5,1,0,0,0,'TimeZone','UTC');
    else
      t1    = refstart;
    end
       
    t2      = datetime(yyyy,10,1,0,0,0,'TimeZone','UTC');  
    tplot   = t1:caldays(1):t2;
    
    % find the first valid promice observation to override t1
    ifirst  = find(~isnan(Ablation.promice(t1:t2)),1,'first');
    t1      = tplot(ifirst);

%     if tplot(ifirst)<t1
%       t1    = tplot(ifirst);
%     end

    iplot   = isbetween(Ablation.Time,t1,t2);    
    x_aws   = Time(iplot);
    x_rcm   = Time(iplot);
    
    y_aws   = Ablation.promice(iplot);

    melt_rcm(:,1)   = cumsum(Ablation.icemodelMelt(iplot));
    melt_rcm(:,2)   = cumsum(Ablation.marMelt(iplot));
    melt_rcm(:,3)   = cumsum(Ablation.merraMelt(iplot));
    melt_rcm(:,4)   = cumsum(Ablation.racmoMelt(iplot));

    roff_rcm(:,1)   = cumsum(Ablation.icemodelRunoff(iplot));
    roff_rcm(:,2)   = cumsum(Ablation.marRunoff(iplot));
    roff_rcm(:,3)   = cumsum(Ablation.merraRunoff(iplot));
    roff_rcm(:,4)   = cumsum(Ablation.racmoRunoff(iplot));

% % don't think this is needed anymore, i_ref is just the first index    
%     t_ref       = datetime(yyyy,6,1,'TimeZone','UTC');
%     i_ref       = x_aws == t_ref;
%     
%     if isnan(y_aws(i_ref))
%         i_ref   = find(i_ref)+find(~isnan(y_aws(find(i_ref):end)),1,'first');
%         t_ref   = x_aws(i_ref);
%     end

    i_ref       = 1;
    y_aws       = y_aws-y_aws(i_ref);
    
    % repeat for the model data
%   i_ref       = x_rcm == t_ref;
    i_ref       = 1;
    melt_rcm    = melt_rcm-melt_rcm(i_ref,:);
    roff_rcm    = roff_rcm-roff_rcm(i_ref,:);

    
    % convert aws ablation to water equivalent
    y_aws       = y_aws.*800/1000;
        
    h.f         = macfig('monitor','mac','size','full');
    h.aws       = plot(x_aws,y_aws); hold on; 

    for n = 1:size(roff_rcm,2)
        h.roff(n)   = plot(x_rcm,roff_rcm(:,n),'Color',dc(n+1,:)); 
        h.melt(n)   = plot(x_rcm,melt_rcm(:,n),':','Color',dc(n+1,:)); 
    end

%     need to add this
%     hline(promice_ablation)
    
    legend( 'AWS Ablation',                         ...
            'IceModel (Runoff)','IceModel (Melt)',  ...
            'MAR (Runoff)', 'MAR (Melt)',           ...
            'MERRA (Runoff)', 'MERRA (Melt)',       ...
            'RACMO (Runoff)', 'RACMO (Melt)',       ...
            'Location','nw');
         
         
   ylabel('ablation [m]');
   
   h.ax  = gca;
   
   figformat(  'linelinewidth',2,'axeslinewidth',1,'labelinterpreter',  ...
               'tex','axesinterpreter','tex','textinterpreter','tex',   ...
               'legendinterpreter','tex','legendlocation','northwest',  ...
               'suppliedaxis',h.ax,'suppliedfigure',h.f);
        
end