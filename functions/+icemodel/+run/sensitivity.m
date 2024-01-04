clean

savedata = true;

% NOTE that for these tests, I turned off the -(opts.sebsolver) thing in
% setopts, which forces skinmodel simulations to use the -solver, to force
% SEBSOLVE to use 1 outer iter. This was intentional so I can actually test
% skinmodel with outer seb iters and w/o. But the actual tests I ran all used
% +solvers, and the maxiter 1,10,100 applied to SKINSOLVE. Note that -solvers
% also force iterflag == false. THIS EXPLAINS WHY I WAS SURPRISED THAT ALL THE
% TESTS RAN WITHOUT ERROR - because it's the SEBSOLVE maxiter that causes
% problems w/SFCTEMP. BUT ACTUALLY MAYBE IT DOESN'T EXPLAIN IT - b/c the problem
% is when outer iterations ARE used w/SFCTEMP, I thought it would fail, but by
% turning off the auto-reset to -solver in setopts, I allowed SFCTEMP to run
% with outer iterations and got no errors. So maybe the errors were related to
% the xTs reset I removed? 

% what should be tested?
% - fmin
% - fcp
% - nspins
% - grain size
% - z0
% - dz_thermal/spectral
% - z_thermal/spectral
% - sebsolver
% - bctype
% - maxiter
% - tol
% - 

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
% DONE behar, 2015, mar, kanm
% DONE behar, 2016, mar, kanm
% DONE upper, 2016, mar, kanm
sitename = {'upperbasin', 'ak4'};
forcings = {'kanl', 'mar'};
userdata = {'merra', 'racmo', 'mar', 'modis', 'kanl', 'none'};
uservars = {'albedo'};
simmodel = {'icemodel'};
simyears = 2016:2016;

%% Set the parameters

% names = {'seb_solver', 'maxiter'};
% values = {[1, 2, 3], [1, 10, 100]};

names = {'z_wind', 'bc_type'};
values = {[2, 10], [1, 2]};

for n = 1:numel(names)
   params.(names{n}) = values{n};
end

%% Create the ensemble

ensemble = ensembleList( ...
   forcings, userdata, uservars, simmodel, simyears, sitename, params);

% Remove cases where forcings == userdata
ensemble.allcombos = ensemble.allcombos( ...
   ensemble.allcombos.forcings ~= ensemble.allcombos.userdata, :);
ensemble.numcombos = height(ensemble.allcombos);

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

   if strcmpi(userdata, forcings)
      continue
   end

   % display the run info
   disp([simmodel ', ' sitename ', '  forcings ', ' userdata ', ' thisyear])

   thisyear = str2double(thisyear);
   
   % set the model options
   opts = icemodel.setopts(simmodel, sitename, thisyear, forcings, ...
      userdata, uservars, savedata);

   % Override the options with the parameters
   for m = 1:numel(names)
      value = ensemble.allcombos.(names{m})(n);
      if isnumeric(params.(names{m}))
         opts.(names{m}) = str2double(value);
      else
         opts.(names{m}) = value;
      end

      % display the run info
      disp([names{m} ' = ' char(value)])
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
   
   % Update the experiment registry?
   
end
