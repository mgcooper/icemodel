function sites = promicesite(source_dir)
   %PROMICESITE Auto-discovered PROMICE station-id namelist (single source of truth).
   %
   %  sites = icemodel.verification.namelists.promicesite()
   %  sites = icemodel.verification.namelists.promicesite(source_dir)
   %
   % Outputs
   %  sites   String ROW of PROMICE station ids discovered from the on-disk
   %          hourly L3 product (<STATION>_hour.nc), e.g. "KAN_L".
   %
   % Role
   %  THE single source of truth for the full PROMICE station list. Unlike the
   %  static snowmipsite list, the PROMICE set is DISCOVERED from the staged L3
   %  product, so it tracks whatever stations are present without a hand-
   %  maintained list (a new station file is picked up automatically). Any
   %  consumer needing the full set - importPromiceSites' default-sites path,
   %  analysis scripts, tests - calls this rather than re-globbing or hardcoding.
   %  For the per-site catalog (zone, eval_target, coords) use
   %  icemodel.verification.helpers.promicesiteinfo.
   %
   % Input
   %  source_dir  Optional PROMICE product directory. Default is
   %              data/verification/promice; the /hour subfolder is appended when
   %              present (the hourly product carries the <STATION>_hour.nc files).
   %
   % See also: icemodel.verification.helpers.promicesiteinfo,
   %  icemodel.verification.setup.importPromiceSites,
   %  icemodel.verification.namelists.snowmipsite

   arguments
      source_dir (1, 1) string = ""
   end

   if source_dir == ""
      source_dir = string(fullfile(icemodel.internal.fullpath('data'), ...
         'verification', 'promice'));
   end
   if isfolder(fullfile(source_dir, 'hour'))
      source_dir = fullfile(source_dir, 'hour');
   end
   files = dir(fullfile(source_dir, '*_hour.nc'));
   if isempty(files)
      error('icemodel:verification:namelists:promicesite:noStations', ...
         'no <STATION>_hour.nc files found under %s', source_dir)
   end
   sites = reshape(string(erase({files.name}, "_hour.nc")), 1, []);
end
