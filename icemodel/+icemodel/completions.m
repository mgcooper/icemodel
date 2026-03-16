function proplist = completions(funcname)
   %COMPLETIONS function completions
   %
   %#codegen
   switch lower(funcname)

      case 'completions'
         proplist = cellstr(icemodel.namelists.completions());

      case 'config'
         proplist = cellstr(icemodel.namelists.config());

      case 'cvconvert'
         proplist = cellstr(icemodel.namelists.cvconvert());

      case 'physicalconstant'
         proplist = cellstr(icemodel.namelists.physicalconstant());

      case 'solver'
         proplist = cellstr(string(icemodel.namelists.solver()));

      case 'testtier'
         proplist = cellstr(icemodel.namelists.testtier());

      case 'testsmbmodel'
         proplist = cellstr(icemodel.namelists.testsmbmodel());

      case 'benchmark'
         proplist = cellstr(icemodel.namelists.benchmark());

      case 'unittest'
         proplist = cellstr(icemodel.namelists.unittest());

      case 'rollingbaseline'
         proplist = cellstr(icemodel.namelists.rollingbaseline());

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
