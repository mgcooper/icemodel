%--------------------------------------------------------------------------
%   Transform mass density from the thermal grid to the spectral grid
%--------------------------------------------------------------------------
function ro_snow = GRIDFORWARD(ro_snow,grid_thermal,grid_spectral)
    ro_snow=interp1(grid_thermal,ro_snow,grid_spectral,'nearest','extrap');