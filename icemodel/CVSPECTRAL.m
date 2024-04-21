function [dz_cv, delz_p, z_node, z_edge, f, z_node_bc, dz_bc] = CVSPECTRAL(Z, dz)
   %CVSPECTRAL build a control volume mesh for spectral model
   %
   % [dz_cv, delz_p, z_node, z_edge, f, z_node_bc, dz_bc] = CVSPECTRAL(Z, dz)
   %
   % See also: CVMESH

   % Note: this is replaced by CVMESH, but retained for reference.

   N = Z / dz;                         % number of nodes [#]
   dz_cv = dz * ones(N, 1);
   delz_b = 0.0;                       % cv width at the boundaries
   dz_bc = [delz_b; dz_cv; delz_b];    % array of c.v. widths including boundaries
   delz_s = 0.5 * dz_bc(1:N+1);        % interface-to-previous point
   delz_n = 0.5 * dz_bc(2:N+2);        % interface-to-next point
   delz_p = delz_s + delz_n;           % distance between grid points
   z_node_bc = [0; cumsum(delz_p)];    % grid point coordinates including boundariess
   z_edge = [0; cumsum(dz_cv)];        % interface coordinates
   z_node = z_node_bc(2:end-1);        % grid nodes
   f = delz_n ./ delz_p;               % interface conductivity weighting factor

   % round the outputs to the nearest milimeter (except f_n)
   f = round(f, 1);
   dz_cv = round(dz_cv, 3);
   dz_bc = round(dz_bc, 3);
   delz_p = round(delz_p, 3);
   z_edge = round(z_edge, 3);
   z_node = round(z_node, 3);
   z_node_bc = round(z_node_bc, 3);
end
