clean

savedata = true;

%-----------------------------------------------------
% upperbasin, 2016, mar, kanl, icemodel
%-----------------------------------------------------
% sitename = {'upperbasin'};
% forcings = {'mar', 'kanl'};
% userdata = {'none','merra','racmo','mar','modis', 'kanl'};
% uservars = {'albedo'};
% simmodel = {'icemodel', 'skinmodel'};
% simyears = 2016:2016;

%-----------------------------------------------------
% DONE behar, 2016, mar, kanm
sitename = {'behar'};
forcings = {'kanm', 'mar'};
userdata = {'none'};
uservars = {'none'};
simmodel = {'skinmodel'};
simyears = 2015:2016;

%% Set the parameters

seb_solver = [1, 2, 3];
maxiter = [1, 10, 100];

params.seb_solver = seb_solver;
params.maxiter = maxiter;
parameters = fieldnames(params);

%% Create the ensemble
ensemble = ensembleList( ...
   forcings, userdata, uservars, simmodel, simyears, sitename, ...
   seb_solver, maxiter);

% Initialize the results
results = ensemble.allcombos;

% Add any custom run-specific settings to the results table (if not in opts)
results.seb_substep = false(height(results), 1);

%% run all combos
for n = 1:ensemble.numcombos

   sitename = char(ensemble.allcombos.sitename(n));
   thisyear = char(ensemble.allcombos.simyears(n));
   simmodel = char(ensemble.allcombos.simmodel(n));
   forcings = char(ensemble.allcombos.forcings(n));
   userdata = char(ensemble.allcombos.userdata(n));
   uservars = char(ensemble.allcombos.uservars(n));

   % this skips cases like 'icemodel_upperBasin_2016_kanl_swap_KANL_albedo'
   % this does not skip the cases where the sitename and forcing are the
   % same, which are the cases where kanm/kanl are used to force basin runs
   if strcmpi(userdata, forcings)
      continue;
   end

   % display the run info
   disp([simmodel ', ' sitename ', '  forcings ', ' userdata ', ' thisyear])

   thisyear = str2double(thisyear);
   
   % set the model options
   opts = icemodel.setopts(simmodel, sitename, thisyear, forcings, ...
      userdata, uservars, savedata);

   % Override the options with the parameters
   for m = 1:numel(parameters)
      value = ensemble.allcombos.(parameters{m})(n);
      if isnumeric(params.(parameters{m}))
         opts.(parameters{m}) = str2double(value);
      else
         opts.(parameters{m}) = value;
      end

      % display the run info
      disp([parameters{m} ' = ' char(value)])
   end

   % RUN THE MODEL
   switch simmodel
      case 'icemodel'
         tic; [ice1, ice2] = icemodel(opts); results.time(n) = toc;

      case 'skinmodel'
         tic; [ice1, ice2] = skinmodel(opts); results.time(n) = toc;
   end

   % post process
   [ice1, ice2, met] = POSTPROC(ice1, ice2, opts, thisyear);

   % quick eval
   [Runoff, Discharge, Catchment] = prepRunoff(opts, ice1);
   AblationDaily = prepAblation(opts, ice1, 'daily');
   
   % this controls the time period over which ablation and runoff are plotted
   t1 = datetime(thisyear,6,1,0,0,0,'TimeZone','UTC');
   t2 = datetime(thisyear,9,1,0,0,0,'TimeZone','UTC');

   % evaluate the runoff
   [h1, data1] = plotRunoff(Runoff, Discharge, Catchment,  ...
      'sitename', sitename, 'userdata', userdata, 'forcingdata', forcings, ...
      't1', t1, 't2', t2, 'plotsurf', opts.simmodel == "skinmodel", ...
      'visible', 'off');
   
   % evaluate the melt
   [h2, data2] = plotPromice(AblationDaily, 'refstart', t1, 'visible', 'off');
   
   % Agreement with field obs
   results.runoff(n) = (data1{end, end} - data1.ADCP(end)) ./ data1.ADCP(end);
   results.ablation(n) = (data2.IcemodelRunoff(end) - data2.AWS(end)) / data2.AWS(end);
   
   % Cumulative melt and runoff
   results.totalmelt(n) = ice1.melt(end);
   results.totalrunoff(n) = ice1.runoff(end);
   
   % Agreement with weather station data
   tmp = icemodel.evaluate(ice1, met);
   fields = fieldnames(tmp);
   for m = 1:numel(fields)
      results.(fields{m})(n) = tmp.(fields{m});
   end
end

results.time = results.time / opts.spinup_loops;

%% save the results
if savedata == true
   pathname = fullfile(getenv('ICEMODELOUTPUTPATH'), 'sensitivity', sitename);
   if ~isfolder(pathname)
      mkdir(pathname)
   end
   filename = fullfile(pathname, [simmodel '_' mkfiledate]);
   save(filename, 'results', 'opts')
   writetable(results, filename, "FileType", "spreadsheet");
end
