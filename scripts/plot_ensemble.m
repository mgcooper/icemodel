clean

pathdata    =  '/Users/coop558/mydata/icemodel/output/v11/run2/';
pathsave    =  '/Users/coop558/mydata/icemodel/figs/v11/run2/';
patheval    =  '/Users/coop558/mydata/icemodel/eval/v11/run2/';

savefigs    =   false;
savedata    =   false;
siteNames   =  {'ak4'};
forcingData =  {'kanl'};
userVars    =  {'albedo'};
meltModel   =  {'icemodel'};
simYears    =  2009:2016;


switch siteNames{1}
   case {'behar','slv1','slv2'}
      
      if strcmp(forcingData,'kanm')
         % this is the case where AWS is the forcing
         userData =  {'none','merra','racmo','mar','modis'};
         order    =  [5 1 2 4 3];    % racmo, mar, merra, aws, modis
         
      elseif strcmp(forcingData,'mar')
         % this is the case where MAR is the forcing
         userData =  {'none','merra','racmo','kanm','modis'};
         order    =  [5 4 2 1 3];    % racmo, mar, merra, aws, modis
      end
      
   case {'ak4','upperBasin'}
      
      if strcmp(forcingData,'kanl')
         % this is the case where AWS is the forcing
         userData =  {'none','merra','racmo','mar','modis'};
         order    =  [5 1 2 4 3];    % racmo, mar, merra, aws, modis
         
      elseif strcmp(forcingData,'mar')
         % this is the case where MAR is the forcing
         userData =  {'none','merra','racmo','kanl','modis'};
         order    =  [5 4 2 1 3];    % racmo, mar, merra, aws, modis
      end
end

ensemble    =  ensembleList(  forcingData,userData,userVars,meltModel,  ...
                              simYears,siteNames);

% need to send icemodel emulators in the following order: racmo, mar,
% merra, aws, modis, to get the legend right. To see the 'order', check
% 'mGroup.userData' in the loop below. note that the case where 'userData'
% is 'none' corresponds to the forcingData albedo (MAR or AWS)

% init the percent dif. should be 8 simulations (racmo, mar, merra, +3
% icemodel emulators, +icemodel aws, +icemodel modis)
pdif  = nan(numel(simYears),8);

for m = 1:numel(simYears)
   
   simYear  =  simYears(m);
   mYear    =  ismember(ensemble.allcombos.simYears,string(simYear));
   mGroup   =  ensemble.allcombos(mYear,:);
   numSims  =  sum(mYear);
   
   % define the start/end time to evaluate the models
   t1       =  datetime(str2double(simYear),6,1,'TimeZone','UTC');
   t2       =  datetime(str2double(simYear),10,1,'TimeZone','UTC');
   
   for n = 1:numSims
      
      siteName    = char(mGroup.siteNames(n));
      simYear     = char(mGroup.simYears(n));
      meltModel   = char(mGroup.meltModel(n));
      forcingData = char(mGroup.forcingData(n));
      userData    = char(mGroup.userData(n));
      userVars    = char(mGroup.userVars(n));
      
      % prep the file name prefix used to save 
      fprfx       = [meltModel '_' siteName '_' simYear '_' forcingData];
      fdata       = [fprfx '_swap_' upper(userData) '_' userVars '.mat'];
      
      % load the data
      if exist([pathdata fdata],'file') == 2
         load([pathdata fdata]);
      else
         continue
      end
      
      if n == 1 || ~exist('Rcol','var') 
         Rcol  =  nan(height(ice1),numSims);
         Rsrf  =  nan(height(ice1),numSims);
         M     =  nan(height(ice1),numSims);
      end
      
      Rcol(:,n)   = ice1.column_runoff;
      Rsrf(:,n)   = ice1.surf_runoff;
      M(:,n)      = ice1.depth_melt;
      
   end
   
%    % for the case where we don't have racmo/merra
%    if size(Rcol,2) == 3
%       Rcol(:,4:5) = nan.*Rcol(:,1:2);
%       Rsrf(:,4:5) = nan.*Rsrf(:,1:2);
%       M(:,4:5)    = nan.*M(:,1:2);
%    end
   
   Rcol  = Rcol(:,order);
   Rsrf  = Rsrf(:,order);
   M     = M(:,order);
   
   
   % put R and M into 'ice1' for compatibility with prepRunoff
   ice1.column_runoff   = Rcol;
   ice1.surf_runoff     = Rsrf;
   ice1.depth_melt      = M;
   
   % load mar, merra, racmo
   [Runoff,Discharge,Catchment]  = prepRunoff(opts,ice1);
   % Ablation                      = prepAblation(opts,ice1);
   
   ltext    = {'ADCP','RACMO','MAR','MERRA','ICE (RACMO)','ICE (MAR)',  ...
               'ICE (MERRA)','ICE (AWS)','ICE (MODIS)'};
   
   h = plotRunoff( Runoff,Discharge,Catchment,  ...
                                          'legendtext',     ltext,      ...
                                          'plotensemble',   true,       ...
                                          'sitename',       siteName,   ...
                                          't1',             t1,         ...
                                          't2',             t2          );
   
   set(h.h(end-1:end),'LineWidth',4)
   title([siteName ', ' simYear ', ' forcingData ' forcing']);
   
   % compute the % difference 
   Rmod        = table2array(h.data(:,2:end));
   Robs        = table2array(h.data(:,1));
   iend        = find(~isnan(Robs),1,'last');
   pdif(m,:)   = 100.*(Rmod(iend,:) - Robs(iend,:))./Robs(iend,:);
   
   if savefigs == true
      if ~exist(pathsave,'dir'); mkdir(pathsave); end
      fsave = [pathsave fprfx '_ensemble.png'];
      exportgraphics(h.f,fsave,'Resolution','300');
   end
   
   if savedata == true
      if ~exist(patheval,'dir'); mkdir(patheval); end
      data  = h.data;
      save([patheval fprfx '_ensemble.mat'],'data'); clear data
   end
   
   
%    pause;
   clear Rcol Rsrf M;
   
   if m<numel(simYears) && n<numel(numSims); close all; end
   
end

varnames = makevalidvarnames(ltext(2:end));
pdif     = array2table(pdif,'VariableNames',varnames);

if savedata == true
   fprfx = [meltModel '_' siteName '_' num2str(simYears(1)) '_'   ...
               num2str(simYears(end)) '_' forcingData];
   save([patheval fprfx '_ensemble_eval.mat'],'pdif');
end

test = nanmean(table2array(pdif),1)

mean(test(1:3))
mean(test(4:end))


% if isequal(userData,{'none','merra','racmo','kanm','modis'}) || ...
%    isequal(userData,{'none','merra','racmo','kanl','modis'})
%    % this is the case where MAR is the forcing
%    
%    order = [5 4 2 1 3];    % racmo, mar, merra, aws, modis
%    
% elseif isequal(userData,{'none','merra','racmo','mar','modis'})
%    % this is the case where the AWS is the forcing
%    
%    order = [5 1 2 4 3];    % racmo, mar, merra, aws, modis
%    
% elseif isequal(userData,{'none','mar','modis'})
%    % this is the case where the AWS is the forcing
%    
%    order = [5 1 4 3 2];    % racmo, mar, merra, aws, modis
% end







