function grid(kwargs)
   % Run grid points. See second section to run grid points for each catchment.

   arguments
      kwargs.savedata = true
      kwargs.sitename = 'sector'
      kwargs.forcings = 'mar'
      kwargs.userdata = 'mar' % don't use "none" - set equal to forcings
      kwargs.uservars = 'albedo'
      kwargs.simmodel = 'skinmodel'
      kwargs.simyears = 2008:2018
      kwargs.gridnums = loadIceMask(varnames="ID")
      kwargs.testname = 'zobs2'
   end
   [savedata, sitename, forcings, userdata, uservars, simmodel, ...
      simyears, gridnums, testname] = dealout(kwargs);

   % Set the path where the output will be saved.
   pathdata = icemodel.setpath('output', kwargs.sitename, kwargs.simmodel, ...
      kwargs.userdata, [], kwargs.testname);

   % Divide the gridpoints into jobs for running in chunks.
   [S, E] = icemodel.chunkgridcell(1, numel(gridnums), 1, 1);

   %% config

   % Initialize the model options.
   opts = icemodel.setopts(simmodel, sitename, simyears, forcings, ...
      userdata, uservars, savedata, [], testname);

   % Display the run information.
   disp([simmodel ' ' userdata ' ' int2str(simyears(1)) ':'   ...
      int2str(simyears(end)) ' ' int2str(S) ':' int2str(E)])

   % Run the model for each point.
   for n = S:E

      thispoint = gridnums(n);

      opts.casename = int2str(thispoint);
      opts.metfname = {fullfile(opts.pathinput, 'met', 'sector', ...
         ['met_' int2str(thispoint) '.mat'])};

      % run the model
      switch simmodel
         case 'icemodel'
            icemodel(opts);
         case 'skinmodel'
            skinmodel(opts);
      end
   end

   % load the met data and run the post processing
   % [ice1, ice2, met] = icemodel.loadresults(opts);
   %
   % % for testing
   % opts.sitename = 'behar';

   %% Override the full grid points with catchment points

   % allsites = {'behar', 'slv1', 'slv2', 'ak4', 'upperbasin'};
   % simyears = {2015:2016, 2015, 2015, 2009:2016, 2016};

   % allsites = {'slv1', 'slv2', 'ak4', 'upperbasin'};
   % simyears = {2015, 2015, 2009:2016, 2016};
   % gridnums = getSiteMetIndex(allsites);

   % run the model for each point (thispoint = 1415, thispoint = 1503 rio behar)
   for n = 1:numel(allsites)

      thissite = allsites{n};
      testname = thissite;
      thispoint = gridnums.(thissite);

      % initialize the model options
      opts = icemodel.setopts(simmodel, sitename, simyears{n}, forcings, ...
         userdata, uservars, savedata, [], testname);

      opts.metfname = {fullfile(opts.pathinput, 'met', 'sector', ...
         ['met_' int2str(thispoint) '.mat'])};

      % run the model
      tic; icemodel(opts); toc
   end
end
