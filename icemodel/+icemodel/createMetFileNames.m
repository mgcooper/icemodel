function metfname = createMetFileNames(opts)
   %CREATEMETFILENAMES Create icemodel met file names for model OPTS.
   %
   %  metfname = icemodel.createMetFileNames(opts)
   %
   % See also: icemodel.setopts icemodel.configureRun

   sitename = opts.sitename;
   forcings = opts.forcings;
   simyears = opts.simyears;

   % Deal with the case where met-station forcing data (as opposed to gridded
   % climate model forcing data) is requested for a nearby catchment by
   % replacing the catchment name in the metfile with the met station name.
   % For example, if sitename=="behar" and forcingdata=="kanm", this sets the
   % metfile name to met_kanm_kanm_YYYY rather than met_behar_kanm_YYYY, to
   % negate the need to create a second (identical) met_behar_kanm_YYYY file.
   if strcmpi(forcings, 'kanl') ...
         && ismember(sitename, {'ak4', 'upperbasin'})
      metname = 'kanl';
   elseif strcmpi(forcings, 'kanm') ...
         && ismember(sitename, {'slv1', 'slv2', 'behar'})
      metname = 'kanm';
   else
      metname = sitename;
   end

   switch opts.dt
      case 900
         dtstr = '15m.mat';
      case 3600
         dtstr = '1hr.mat';
      otherwise
         error('unsupported dt for met file naming: %g', opts.dt)
   end

   metfname = cell(1, numel(simyears));
   for n = 1:numel(simyears)
      simyear = num2str(simyears(n));
      metfname{n} = ['met_' metname '_' forcings '_' simyear '_' dtstr];
   end
end
