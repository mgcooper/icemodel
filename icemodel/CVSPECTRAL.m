%--------------------------------------------------------------------------
%   build a control volume with quasi-exponential grid spacing
%--------------------------------------------------------------------------

function [  dz_p,                                                       ...
            dz_pbc,                                                     ...
            delz_p,                                                     ...
            f_n,                                                        ...
            z_crds,                                                     ...
            z_walls,                                                    ...
            grid_spectral ]  =  CVSPECTRAL(JJ_spectral,deltaz_spectral)
%--------------------------------------------------------------------------

    deltaz_b        =   0.0;                            % width of c.v. at the boundaries
    dz_p            =   deltaz_spectral.*ones(JJ_spectral,1);
    dz_pbc          =   [deltaz_b;dz_p;deltaz_b];       % array of c.v. widths including boundaries
    delz_neg        =   0.5.*dz_pbc(1:JJ_spectral+1);   % interface-to-previous point
    delz_pos        =   0.5.*dz_pbc(2:JJ_spectral+2);   % interface-to-next point
    delz_p          =   delz_neg + delz_pos;            % distance between grid points
    z_crds          =   [0;cumsum(delz_p)];             % grid point coordinates       
    z_walls         =   [0;cumsum(dz_p)];               % interface coordinates
    f_n             =   delz_pos./delz_p;               % interface conductivity weighting factor

    grid_spectral   =   z_crds(2:end-1);
    
  % round the outputs to the nearest milimeter (except f_n)       
    dz_p            =   roundn(dz_p,-3);
    dz_pbc          =   roundn(dz_pbc,-3);
    delz_p          =   roundn(delz_p,-3);
    f_n             =   roundn(f_n,-1);
    z_crds          =   roundn(z_crds,-3);
    z_walls         =   roundn(z_walls,-3);
    grid_spectral   =   roundn(grid_spectral,-3);
      
