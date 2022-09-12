clean

% _new as of july 2022 to resave the v10/run1

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% icemodel
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
savedata    = true;
pathdata    =  '/Users/coop558/mydata/icemodel/output/v10/run1/';
siteNames   =  {'behar','slv1','slv2'};
forcingData =  {'kanm','mar'};
userVars    =  {'albedo'};
meltModel   =  {'icemodel'};
simYears    =  2015:2016;

cd(pathdata)


switch siteNames{1}
   case {'behar','slv1','slv2'}
         userData =  {'none','merra','racmo','kanm','modis','mar'};
   case {'ak4','upperBasin'}
      userData =  {'none','merra','racmo','mar','modis','kanl'};
end

ensemble    =  ensembleList(  forcingData,userData,userVars,meltModel,  ...
                              simYears,siteNames);


for m = 1:numel(simYears)

   simYear  =  simYears(m);
   mYear    =  ismember(ensemble.allcombos.simYears,string(simYear));
   mGroup   =  ensemble.allcombos(mYear,:);
   numSims  =  sum(mYear);

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
      
      if strcmpi(userData,forcingData) || exist(fdata,'file') ~= 2
         continue
      end

      % load the data
      load([pathdata fdata]);

      fields   = fieldnames(ice2);

      t1       = datetime(year(enbal.Time(1)),1,1,0,0,0);
      t2       = datetime(year(enbal.Time(1)),12,31,23,45,0);
      oldtime  = tocol(t1:minutes(15):t2);
      newtime  = enbal.Time;

      oldtime.TimeZone = 'UTC';

      for mm = 1:numel(fields)

         thisfield = fields{mm};

         if strcmp(thisfield,'Z') || strcmp(thisfield,'LCflag')
            continue
         end
         dat               = transpose(ice2.(thisfield));
         newdat            = interp1(oldtime,dat,newtime);
         ice2.(thisfield)  = transpose(newdat);

         if ismember(thisfield,{'f_ice','f_liq','k_vap','k_eff','Qp'})
            ice2.(thisfield) = round(ice2.(thisfield),5);
         elseif ismember(thisfield,{'cp_sno','ro_sno','Tice'})
            ice2.(thisfield) = round(ice2.(thisfield),1);
         elseif ismember(thisfield,{'Qp','Terr','dH','df_liq','df_drn'})
            ice2.(thisfield) = round(ice2.(thisfield),8);
         end
      end

      if savedata == true
         save([pathdata fdata],'enbal','ice1','ice2','opts');
      end

   end  
end

