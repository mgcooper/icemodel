%--------------------------------------------------------------------------
% 	define a z-depth array, the centers of each grid cell
%--------------------------------------------------------------------------
function [z_with_bc,z_without_bc] = GETZ(deltaz,nz)
%--------------------------------------------------------------------------

% Make a z-depth array, in meters.  These are z-levels, the centers
%   of each grid cell.  z_with_bc includes the top and bottom edges of
%   the top and bottom grid cells.

% z_with_bc is identical to y_crds i.e. the control volume centers
    z_with_bc     = nan(nz+2,1);
    z_without_bc  = nan(nz,1);
    z_with_bc(1)  = 0.0;
    for k=2:nz+1
        z_with_bc(k) = deltaz * real(k-1) - deltaz / 2.0;
        z_without_bc(k-1) = z_with_bc(k);
    end
    z_with_bc(nz+2) = deltaz * real(nz);
  