function ro_snow = GRIDFORWARD(ro_snow,grid_thermal,grid_spectral)
   %GRIDFORWARD Transform mass density from thermal grid to spectral grid
   ro_snow=interp1(grid_thermal,ro_snow,grid_spectral,'nearest','extrap');
end
