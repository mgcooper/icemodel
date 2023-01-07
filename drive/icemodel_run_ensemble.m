clean

% this demonstrates how to set up an ensemble of runs at one site

%------------------------------------------------------------------------------
%  set the run-specific model configuration
%------------------------------------------------------------------------------
savedata    =  true;
dryrun      =  false;
sitenames   =  {'ak4','upperbasin','behar','slv1','slv2'};
forcingdata =  {'mar','kanl','kanm'};
userdata    =  {'racmo','mar','merra','kanm','kanl','modis','none'};
uservars    =  {'albedo'};
meltmodels  =  {'skinmodel'};
startyear   =  2009;
endyear     =  2016;
simyears    =  startyear:endyear;

% this builds an array of all combinations of the above settings.
ensemble    =  ensembleList(  forcingdata,userdata,uservars,meltmodels, ...
                              simyears,sitenames);

% setenv('ICEMODELOUTPUTPATH',[getenv('HOME') '/work/data/icemodel/output/v10b/']);
setenv('ICEMODELOUTPUTPATH','/Volumes/Samsung_T5b/icemodel/output/ensemble/v10b/');
%------------------------------------------------------------------------------
% notes on ensemble options
%------------------------------------------------------------------------------

% valid userdata:
% {'racmo','mar','merra','kanm','kanl','modis','none'}

% valid sitenames:
% {'upperbasin','behar','slv1','slv2','ak4'}

% valid forcingdata:
% {'mar','kanm','kanl'}

%------------------------------------------------------------------------------
% run all combos
%------------------------------------------------------------------------------
nruns = 0;
for n = 1:ensemble.numcombos
   
   sitename    = char(ensemble.allcombos.sitenames(n));
   simyear     = char(ensemble.allcombos.simyears(n));
   meltmodel   = char(ensemble.allcombos.meltmodels(n));
   forcingdata = char(ensemble.allcombos.forcingdata(n));
   userdata    = char(ensemble.allcombos.userdata(n));
   uservars    = char(ensemble.allcombos.uservars(n));
 % uservars    = cellstr(ensemble.allcombos.uservars(n));
   
   % this skips cases like 'icemodel_upperbasin_2016_kanl_swap_kanl_albedo'
   % this does not skip the cases where the sitename and forcing are the
   % same, which are the cases where kanm/kanl are used to force basin runs
   if strcmpi(userdata,forcingdata)                                     ||    ...
      (strcmpi(sitename,'ak4') && strcmpi(userdata,'kanm'))             ||    ...
      (strcmpi(sitename,'upperbasin') && strcmpi(userdata,'kanm'))      ||    ...
      (strcmpi(sitename,'behar') && strcmpi(userdata,'kanl'))           ||    ...
      (strcmpi(sitename,'slv1') && strcmpi(userdata,'kanl'))            ||    ...
      (strcmpi(sitename,'slv2') && strcmpi(userdata,'kanl'))            ||    ...
      (strcmpi(sitename,'slv1') && str2double(simyear) ~= 2015)         ||    ...
      (strcmpi(sitename,'slv2') && str2double(simyear) ~= 2015)         ||    ...
      (strcmpi(sitename,'behar') && str2double(simyear) < 2015)         ||    ...
      (strcmpi(sitename,'upperbasin') && str2double(simyear) < 2016)    ||    ...
      (strcmpi(userdata,'racmo') && str2double(simyear) < 2012)         ||    ...
      (strcmpi(userdata,'merra') && str2double(simyear) < 2012)
      continue;
   end

   % set the model options
   opts = icemodel_opts(sitename,meltmodel,forcingdata,userdata,uservars,   ...
                        str2double(simyear),str2double(simyear));

   if strcmp(opts.msg,'met file not found')
      continue;
   end

   % display the model configuration
   disp([meltmodel ', ' sitename ', '  forcingdata ', ' userdata ', ' simyear])
   nruns = nruns + 1;
%    disp(opts.fsave)
%------------------------------------------------------------------------------
%     RUN THE MODEL
%------------------------------------------------------------------------------
   if dryrun == false

      if opts.meltmodel == "skinmodel"
         [ice1,ice2,met,opts] = skinmodel(opts);
      else
         [ice1,ice2,met,opts] = icemodel(opts);
      end

      % save the data
      if savedata == true
         if ~exist(opts.pathoutput,'dir'); mkdir(opts.pathoutput); end 
         save(opts.fsave,'ice1','ice2','opts');
      end

   end

end
nruns
%------------------------------------------------------------------------------
% quick eval
%------------------------------------------------------------------------------
if dryrun == false
   [Runoff,Discharge,Catchment]  = prepRunoff(opts,ice1);
   AblationHourly                = prepAblation(opts,ice1,'hourly');
   AblationDaily                 = prepAblation(opts,ice1,'daily');
   
   % plot the runoff
   if opts.meltmodel == "skinmodel"
	   h1 = plotRunoff(Runoff,Discharge,Catchment,'plotsurf',true,'sitename',  ...
                        sitename,'userdata',userdata,'forcingdata',forcingdata);
   else
      h1 = plotRunoff(Runoff,Discharge,Catchment,'sitename',sitename,         ...
                        'userdata',userdata,'forcingdata',forcingdata);
   end
   
   % PLOT ABLATION
   t1 = datetime(str2double(simyear),6,1,0,0,0,'TimeZone','UTC');
   plotPromice(AblationHourly,'refstart',t1);
   plotPromice(AblationDaily,'refstart',t1);
   
   % plot enbal - need option to plot met station data when forcing is rcm 
   plotEnbal(ice1,met);
end