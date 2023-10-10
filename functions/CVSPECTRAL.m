function [dz_p, dz_pbc, delz_p, f_n, z_crds, z_wall, z_node ] = ...
   CVSPECTRAL(N, dz)
%CVSPECTRAL build a control volume mesh for spectral model

dz_b = 0.0;                      % cv width at the boundaries
dz_p = dz.*ones(N,1);
dz_pbc = [dz_b;dz_p;dz_b];       % array of c.v. widths including boundaries
delz_s = 0.5.*dz_pbc(1:N+1);     % interface-to-previous point
delz_n = 0.5.*dz_pbc(2:N+2);     % interface-to-next point
delz_p = delz_s + delz_n;        % distance between grid points
z_crds = [0;cumsum(delz_p)];     % grid point coordinates including boundariess
z_wall = [0;cumsum(dz_p)];       % interface coordinates
z_node = z_crds(2:end-1);        % grid nodes
f_n = delz_n./delz_p;            % interface conductivity weighting factor

% round the outputs to the nearest milimeter (except f_n)
f_n = round(f_n,1);
dz_p = round(dz_p,3);
dz_pbc = round(dz_pbc,3);
delz_p = round(delz_p,3);
z_crds = round(z_crds,3);
z_wall = round(z_wall,3);
z_node = round(z_node,3);

