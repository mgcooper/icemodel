function proplist = completions(funcname)
   %COMPLETIONS function completions
   %
   %#codegen
   switch lower(funcname)

      case 'completions'
         tmp = dir( ...
            fullfile(icemodel.internal.fullpath, 'icemodel', '+icemodel', '*.m'));
         proplist = strrep({tmp.name}, '.m', '');

      case 'config'
         proplist = {'demo'}.';
         % To re-enable completions:
         % {"name":"casename", "kind":"namevalue", "type":["choices=icemodel.completions('config')"]},

      case 'cvconvert'
         proplist = {'volumefraction', 'bulkdensity', 'mass', 'volume', ...
            'massfraction', 'totaldensity'}.';

      case 'physicalconstant'
         proplist = {'Tf', 'Lv', 'Lf', 'Ls', 'ro_air', 'ro_ice', 'ro_liq', ...
            'cp_air', 'cp_liq', 'cp_ice', 'k_liq', 'k_ice', 'Rd', 'Rv', 'SB', ...
            'emiss', 'gravity', 'kappa', 'kappa_p', 'epsilon', 'P0', 'es0', ...
            'S0', 'N0', 'psychro', 'dalr', 'malr', 'fcp', 'scale_ht', 'hrsperday', ...
            'secperhr', 'roLv', 'roLs', 'roLf', 'cv_air', 'cv_ice', 'cv_liq', ...
            'ro_iwe', 'ro_wie', 'emissSB', 'fcpsq', 'secperday', 'c0'}.';

      case 'solver'
         proplist = cellstr(string(icemodel.namelists.solver()));

      case 'testtier'
         proplist = {'smoke', 'full', 'all'}.';

      case 'testsmbmodel'
         proplist = [{'all'}; cellstr(icemodel.namelists.smbmodel())];

      case 'rollingbaseline'
         proplist = {'rolling'}.';

      case 'smbmodel'
         proplist = cellstr(icemodel.namelists.smbmodel());

      case 'forcings'
         proplist = cellstr(icemodel.namelists.forcings());

      case 'sitename'
         proplist = cellstr(icemodel.namelists.sitename());

      case 'userdata'
         proplist = cellstr(icemodel.namelists.userdata());

      case 'uservars'
         proplist = cellstr(icemodel.namelists.uservars());
   end
end
