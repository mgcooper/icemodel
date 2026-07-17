function filename = metfilename(site, forcings, t1, t2, dt)
   %METFILENAME Build a standard icemodel met-file name.
   %
   %  filename = icemodel.forcing.helpers.metfilename(site, forcings, ...
   %     t1, t2, dt)
   %
   % Two naming forms are supported, matching what
   % icemodel.createMetFileNames parses on the read side:
   %
   %  Window form (t1 and t2 are datetimes):
   %     met_<site>_<forcings>_<YYYYMMDD>_<YYYYMMDD>_<dt>.mat
   %
   %  Legacy per-year form (t1 is a year number, t2 is []):
   %     met_<site>_<forcings>_<YYYY>_<dt>.mat
   %
   % DT is the forcing timestep in seconds (900, 1800, or 3600) or the literal
   % filename suffix ("15m", "30m", or "1hr"). The 30-minute form supports the
   % proven native Samimi Dye-2 cadence; repository writers still default to 15m.
   %
   % See also: icemodel.createMetFileNames,
   %  icemodel.forcing.helpers.writemet

   arguments
      site (1, 1) string
      forcings (1, 1) string
      t1
      t2
      dt
   end

   try
      dtstr = char(icemodel.forcing.helpers.metTimestepSuffix(dt));
   catch err
      if string(err.identifier) ~= ...
            "icemodel:forcing:metTimestepSuffix:unsupportedTimestep"
         rethrow(err)
      end
      % Preserve the established public writer-side error identifier.
      error('icemodel:forcing:metfilename:unsupportedTimestep', ...
         'unsupported dt for met file naming')
   end

   if isa(t1, 'datetime')
      if ~isa(t2, 'datetime')
         error('icemodel:forcing:metfilename:badWindow', ...
            'window form requires t1 and t2 both datetimes');
      end
      filename = sprintf('met_%s_%s_%s_%s_%s.mat', site, forcings, ...
         char(t1, 'yyyyMMdd'), char(t2, 'yyyyMMdd'), dtstr);

   elseif isnumeric(t1) && isscalar(t1) && isempty(t2)
      filename = sprintf('met_%s_%s_%d_%s.mat', site, forcings, t1, dtstr);

   else
      error('icemodel:forcing:metfilename:badTimeArguments', ...
         't1/t2 must be datetimes (window form) or a year and [] (per-year form)');
   end
   filename = string(filename);
end
