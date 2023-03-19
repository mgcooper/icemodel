function [met,opts] = METINIT(opts)
    
% re-activate this when/if i put all the metfiles in the input directory
%load([opts.path.input 'met/' opts.metfname],'met')

if opts.sitename == "sector"

   % load the first met file
   load(opts.metfname,'met');

   if opts.calendar_type == "noleap"
      met = rmleapinds(met);
   end
   
   % subset the met file to the requested simyears
   met = met(ismember(year(met.Time),opts.simyears),:);
   opts.maxiter = height(met)/opts.numyears;
   opts.dt = seconds(met.Time(2)-met.Time(1));

   % check for bad albedo data
   bi = find(met.modis<=0 | met.modis>=1);
   if ~isempty(bi)
      met.modis(bi) = met.albedo(bi);
      warning('bad albedo')
   end

   if string(opts.userdata) == "modis"
      met.albedo = met.modis;
   end
   
   return
   
end
   
   
% these load met-station forcing data for nearby catchments

% slv1
if strcmp(opts.metfname,      'met_slv1_kanm_2015_15m.mat')
   load([opts.path.metdata    'met_KANM_KANM_2015_15m.mat'],'met')
elseif strcmp(opts.metfname,  'met_slv1_kanm_2015_1hr.mat')
   load([opts.path.metdata    'met_KANM_KANM_2015_1hr.mat'],'met')

% slv2
elseif strcmp(opts.metfname,  'met_slv2_kanm_2015_15m.mat')
   load([opts.path.metdata    'met_KANM_KANM_2015_15m.mat'],'met')
elseif strcmp(opts.metfname,  'met_slv2_kanm_2015_1hr.mat')
   load([opts.path.metdata    'met_KANM_KANM_2015_1hr.mat'],'met')

% behar   
elseif strcmp(opts.metfname,  'met_behar_kanm_2015_1hr.mat')
   load([opts.path.metdata    'met_KANM_KANM_2015_1hr.mat'],'met')
elseif strcmp(opts.metfname,  'met_behar_kanm_2016_1hr.mat')
   load([opts.path.metdata    'met_KANM_KANM_2016_1hr.mat'],'met')   
elseif strcmp(opts.metfname,  'met_behar_kanm_2015_15m.mat')
   load([opts.path.metdata    'met_KANM_KANM_2015_15m.mat'],'met')
elseif strcmp(opts.metfname,  'met_behar_kanm_2016_15m.mat')
   load([opts.path.metdata    'met_KANM_KANM_2016_15m.mat'],'met')

% ak4
elseif strcmp(opts.metfname,  'met_ak4_kanl_2009_15m.mat')
   load([opts.path.metdata    'met_KANL_KANL_2009_15m.mat'],'met')
elseif strcmp(opts.metfname,  'met_ak4_kanl_2010_15m.mat')
   load([opts.path.metdata    'met_KANL_KANL_2010_15m.mat'],'met')
elseif strcmp(opts.metfname,  'met_ak4_kanl_2011_15m.mat')
   load([opts.path.metdata    'met_KANL_KANL_2011_15m.mat'],'met')
elseif strcmp(opts.metfname,  'met_ak4_kanl_2012_15m.mat')
   load([opts.path.metdata    'met_KANL_KANL_2012_15m.mat'],'met')
elseif strcmp(opts.metfname,  'met_ak4_kanl_2013_15m.mat')
   load([opts.path.metdata    'met_KANL_KANL_2013_15m.mat'],'met')   
elseif strcmp(opts.metfname,  'met_ak4_kanl_2014_15m.mat')
   load([opts.path.metdata    'met_KANL_KANL_2014_15m.mat'],'met')
elseif strcmp(opts.metfname,  'met_ak4_kanl_2015_15m.mat')
   load([opts.path.metdata    'met_KANL_KANL_2015_15m.mat'],'met')
elseif strcmp(opts.metfname,  'met_ak4_kanl_2016_15m.mat')
   load([opts.path.metdata    'met_KANL_KANL_2016_15m.mat'],'met')   
elseif strcmp(opts.metfname,  'met_ak4_kanl_2009_1hr.mat')
   load([opts.path.metdata    'met_KANL_KANL_2009_1hr.mat'],'met')
elseif strcmp(opts.metfname,  'met_ak4_kanl_2010_1hr.mat')
   load([opts.path.metdata    'met_KANL_KANL_2010_1hr.mat'],'met')
elseif strcmp(opts.metfname,  'met_ak4_kanl_2011_1hr.mat')
   load([opts.path.metdata    'met_KANL_KANL_2011_1hr.mat'],'met')
elseif strcmp(opts.metfname,  'met_ak4_kanl_2012_1hr.mat')
   load([opts.path.metdata    'met_KANL_KANL_2012_1hr.mat'],'met')
elseif strcmp(opts.metfname,  'met_ak4_kanl_2013_1hr.mat')
   load([opts.path.metdata    'met_KANL_KANL_2013_1hr.mat'],'met')
elseif strcmp(opts.metfname,  'met_ak4_kanl_2014_1hr.mat')
   load([opts.path.metdata    'met_KANL_KANL_2014_1hr.mat'],'met')
elseif strcmp(opts.metfname,  'met_ak4_kanl_2015_1hr.mat')
   load([opts.path.metdata    'met_KANL_KANL_2015_1hr.mat'],'met')
elseif strcmp(opts.metfname,  'met_ak4_kanl_2016_1hr.mat')
   load([opts.path.metdata    'met_KANL_KANL_2016_1hr.mat'],'met')

% upperBasin
elseif strcmp(opts.metfname,  'met_upperBasin_kanl_2016_1hr.mat')
   load([opts.path.metdata    'met_KANL_KANL_2016_1hr.mat'],'met')
elseif strcmp(opts.metfname,  'met_upperBasin_kanl_2016_15m.mat')
   load([opts.path.metdata    'met_KANL_KANL_2016_15m.mat'],'met')

% final check - if the metfile does not exist, return w/ error   
elseif exist(opts.metfname,'file') ~= 2
   met         =  [];
   opts.error  =  'metfile does not exist';
   return;

% catchment-scale RCM forcing file
else
   load([opts.path.metdata opts.metfname],'met')
end

% next swaps out a variable from userData in place of the default met data
if opts.useUserData == true
    
   simYear     =   opts.yyyy;
   sitename    =   opts.sitename;
   userData    =   opts.forcingUserData;
   userVars    =   opts.forcingUserVars;
   userPath    =   opts.path.userdata;
   numUserVars =   numel(userVars);
    
% load the forcingUserData file and retime to match the metfile
   userFile    =   [userPath userData '_' sitename '_' simYear '.mat'];
   
   if exist(userFile,'file') ~= 2
      opts.error  =  'metfile does not exist';
      return;
   end
   
   load(userFile,'Data');
   if Data.Time(2)-Data.Time(1) ~= met.Time(2)-met.Time(1)
      Data  = retime(Data,met.Time,'linear');
   end

   % swap out the data in the metfile for the user data
   for n = 1:numUserVars

      thisVar  = userVars{n};
      swapVar  = thisVar;

      % MODIS requires changing the variable name to albedo
      if lower(userData) == "modis" && lower(thisVar) == "albedo"
         thisVar = 'albedo';
         swapVar = 'modis';
      end

      met.(thisVar) = Data.(swapVar);
   end
end

if isempty(met.Time.TimeZone)
   met.Time.TimeZone = 'UTC';
end

%%%% test
% met = retime(met,'hourly','mean');
%%%% test

% met            = rmleapinds(met);
opts.maxiter   = size(met,1);
opts.dt        = seconds(met.Time(2)-met.Time(1));
   
% % check:
%    figure; plot(Data.Time,Data.albedo); hold on; 
%    plot(Data.Time,Data.modis); legend('MAR','MODIS');
%    oldalbedo  = met.albedo;
%    figure; plot(met.Time,oldalbedo); hold on;
%    plot(met.Time,met.albedo); legend('old albedo','new albedo');

% turned this off b/c it makes the model slower, better to use 15 m metfiles
% interpolate the forcings to the chosen timestep
%    if opts.dt ~= 3600
%       switch opts.dt
%          case 900
%             met = interpMet(met,'15m');
%       end
%       met.date =  datenum(met.Time);
%    end


% older checks:
% % convert air temperature from celsius to kelvin
%     met(:,5) = met(:,5) + Tfp;
%     
% % convert srf temperature from celsius to kelvin
%     met(:,14) = met(:,14) + Tfp;
    
% % TEMP - needed for new timetable format
%     if istimetable(met)
%         met     = timetable2table(met);
%         vars    = met.Properties.VariableNames;
%         met     = table2array(met(:,2:end));
%         met     = met(:,[1:14 19]); % 
%     end
    
    
% % this means the timestep is < 1hr
%     if nhrs<1
%         t0      =   datenum([met(1,1:4) 0 0]);
%         tF      =   datenum([met(end,1:4) 0 0]);
%         dt      =   opts.dt/3600/24;
%         told    =   t0:1/24:tF;
%         tnew    =   t0:dt:tF;
%         metold  =   met(:,5:end);
%         metnew  =   interp1(told,metold,tnew);
%         tnew    =   datevec(tnew);
%         met     =   [tnew(:,1:4) metnew];
%     end
    

% % this is no longer needed if met_15m or met_1hr is used    
%     one_hr  =   1/24;
%     t1      =   datenum([met(1  , 1:4) 0 0]);
%     t2      =   datenum([met(end, 1:4) 0 0]);
%     tmet    =   t1:one_hr:t2; tmet = tmet(:);
%     met     =   [met tmet];
    
% % was gonna use this to deal with precision of the running model time
% calendar, but that problem is fixed by the met interpolation method
% one_hr  =   round(1/24,10);
   
   
% % the original albedo emulator swap   
%     switch opts.albedoData
%         case 'merra'
%             load([opts.path.input 'external/merra_' opts.yyyy],'merra');
%             met.albedo =   merra.albedo;
%         case 'mar'
%             load([opts.path.input 'external/mar_' opts.yyyy],'mar');
%             met.albedo =   mar.albedo;
%         case 'racmo'
%             load([opts.path.input 'external/racmo_' opts.yyyy],'racmo');
%             met.albedo =   racmo.albedo;
%         otherwise
%     end 
       
    