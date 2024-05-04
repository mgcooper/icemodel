function proplist = completions(funcname)
   %COMPLETIONS function completions
   %
   %#codegen
   switch lower(funcname)

      case 'completions'
         tmp = dir( ...
            fullfile(icemodel.projectpath, 'functions', '+icemodel', '*.m'));
         proplist = strrep({tmp.name}, '.m', '');

      case 'cvconvert'
         proplist = {'volumefraction', 'bulkdensity', 'mass', 'volume', ...
            'massfraction', 'totaldensity'}.';

      case 'physicalconstant'
         proplist = {'Tf', 'Lv', 'Lf', 'Ls', 'ro_air', 'ro_ice', 'ro_liq', ...
            'cp_air', 'cp_liq', 'cp_ice', 'k_liq', 'k_ice', 'Rd', 'Rv', 'SB', ...
            'emiss', 'gravity', 'kappa', 'kappa_p', 'epsilon', 'P0', 'es0', ...
            'S0', 'N0', 'psychro', 'dalr', 'malr', 'fcp', 'scale_ht', 'hrsperday', ...
            'secperhr', 'roLv', 'roLs', 'roLf', 'cv_air', 'cv_ice', 'cv_liq', ...
            'ro_iwe', 'ro_wie', 'emissSB', 'fcpsq', 'secperday'}.';
   end
end
