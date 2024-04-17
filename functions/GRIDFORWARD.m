function ro_sno = GRIDFORWARD(ro_sno, grid_thermal, grid_spectral)
   %GRIDFORWARD Transform mass density from thermal grid to spectral grid
   ro_sno = interp1(grid_thermal, ro_sno, grid_spectral, 'nearest', 'extrap');
end
