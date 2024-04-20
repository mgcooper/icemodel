function [dz, delz, z_node, z_edge, f, z_node_bc] = CVMESH(Z, dz, g)
   %CVMESH Compute cell edges and nodes for control volume mesh
   %
   % [dz, delz, z_node, z_edge, f, z_node_bc] = CVMESH(Z, dz) Returns CV
   % thickness (dz), node-to-node distance (delz), node positions (z_node),
   % and edge positions (z_edge) for a domain with upper boundary coordinate
   % z=0, lower boundary coordinate z=Z, and constant CV thickness dz. An
   % interpolation factor f for each CV interface (edge) is also returned using
   % the harmonic mean. Control volume nodes are placed midway between faces,
   % following Patankar's "Practice B".
   %
   % [dz, delz, z_node, z_edge, f, z_node_bc] = CVMESH(Z, dz, g) Uses growth
   % factor 'g' to construct a mesh with exponentially increasing thickness. If
   % g is not provided, a uniform thickness mesh is returned.
   %
   % Inputs:
   %   Z  : Domain thickness [m]
   %   dz : Desired CV thickness [m]
   %   g  : Growth factor by which CV thickness increases each layer [1]
   %
   % Outputs:
   %   dz        : Thickness of CVs [m]
   %   delz      : Distance between CV nodes, node-to-node distance [m]
   %   z_node    : Position of CV nodes [m]
   %   z_edge    : Position of CV edges [m]
   %   f         : Interpolation factor for each CV edge (interface) []
   %   z_node_bc : Position of CV nodes including fictitious boundary nodes [m]
   %
   % The function assumes a Dirichlet boundary condition at the upper boundary
   % and a Neumann boundary condition at the lower boundary. The upper boundary
   % is at z = 0 and the lower boundary is at z = Z. The function computes a
   % constant CV thickness (dz), but adjusts for half CVs at the boundaries.
   % The harmonic mean is used for computing the interpolation factor f for
   % each CV edge.
   %
   % Example usage:
   %   [dz, delz, z_node, z_edge, f, z_node_bc] = CVTHERMAL(10, 1);
   %
   % See also

   % Compute mesh edge positions
   if nargin < 3 || g == 1

      % Linear mesh
      z_edge = (0:dz:Z)';

   else
      % Exponential mesh
      assert(isscalar(g) && isnumeric(g))
      z_edge = EXPMESH(Z, dz, g);

      % z_edge = cumsum([0, dz * g .^ (0:floor(log(Z / dz) / log(g)))])';
      % z_edge = [z_edge(z_edge < Z - eps(z_edge(end))); Z];
   end
   N = numel(z_edge) - 1;

   % Compute node positions and cell thicknesses (dz)
   z_node = z_edge(1:N) + diff(z_edge) / 2;

   % Append fictitious boundary nodes
   z_node_bc = [0; z_node; Z];

   % Compute distance between nodes
   delz = diff(z_node_bc);

   % Compute distance between edges
   dz = diff(z_edge);

   % Compute interpolation factor f for cell edges
   f = zeros(N + 1, 1);
   for n = 1:N+1
      f(n) = (z_node_bc(n+1) - z_edge(n)) / (z_node_bc(n+1) - z_node_bc(n));
   end

   % Plot the mesh
   % icemodel.plotmesh(Z, dz, g)

   % % Round all quantities to the nearest milimeter (except f)
   % [dz, delz, z_node, z_edge, z_node_bc, f] = deal( ...
   %    round(dz, 3), round(delz, 3), round(z_node, 3), round(z_edge, 3), ...
   %    round(z_node_bc, 3), round(f, 3));

   % % These quantities are not computed, but may be useful.
   % dz_bc = [delz_b; dz; delz_b]; % array of c.v. widths including boundaries
   % delz_s = 0.5 * dz_bc(1:N+1); % interface-to-previous point
   % delz_n = 0.5 * dz_bc(2:N+2); % interface-to-next point

   % For reference
   % z_walls = [0; cumsum(dz_spect); sum(dz_spect) + dz_spect(end)];
   % z_walls = round([0; z_spect + dz_spect/2; z_spect(end) + 3*dz_spect/2], 3);
   % z_walls = [z_walls; z_walls(end) + (z_walls(end) - z_walls(end-1))];
   % z_spect = z_walls(1:end-2) + diff(z_walls(1:end-1)) / 2;
end

%% Exponential mesh
function z_edge = EXPMESH(Z, dz, g)
   %EXPMESH Compute mesh cell edges and nodes with exponential thickness.
   %
   % Inputs:
   %   Z    : Domain thickness [m]
   %   dz   : Top CV thickness [m]
   %   g    : Growth factor by which CV thickness increases each layer [1]
   %
   % Outputs:
   %   (Outputs are the same as in the main function but adapted to the
   %    exponential mesh.)
   %
   % See also: CVMESH

   % The "insert a bottom layer with dz thickness" idea was to simplify
   % remeshing to enable an exponential thermal grid, but the easier approach
   % is to use an exponential spectral grid which does not require that, so
   % I commented it out but kept for reference.

   % Preallocate arrays (with an arbitrary large size, will trim later)
   z_edge = zeros(10000, 1);
   z_edge(1) = 0; % Top edge
   n = 0;         % Layer counter

   % Compute edge positions until one exceeds Z
   while true
      n = n + 1;
      if z_edge(n) + dz * g^(n-1) + n * eps(Z) < Z % - dz
         z_edge(n + 1) = z_edge(n) + dz * g^(n-1);
      else
         z_edge(n + 1) = Z; % - dz;
         break
      end
   end
   % Trim to actual number of layers and insert a bottom layer with dz thickness
   % z_edge = [z_edge(1:n+1); z_edge(n+1) + dz];
   z_edge = z_edge(1:n+1);
end

%% Notes
%
% N = Number of internal nodes / c.v.'s (5 in example below)
% N+1 = Number of CV interfaces (6 in example below)
% q = heat flux
% Sc = source term (Sp omitted)
% gamma = conductivity, a function of T defined at T1, T2, ..., TN
% k = the interface conductivity and is calculated from the values at T1,
% T2, ..., TN using eq. 4.9 and the interpolation factor f.
%
%  /////////
%  ----o---- Tsfc ---       (upper boundary T = Tsfc, dz_pbc = 0, delz = 0)
% k1 .....         |        interface 1: k1 = k(T1)
%      |          dz_1
%      o     T1   ---       q1 = k1(Tsfc-T1)/dz_1 - k2(T1-T2)/dz_2 + Sc1*dz
%      |           |
% k2 .....        dz_2      interface 2: k2 = 2*k(T1)*k(T2) / (k(T1)+k(T2))
%      |           |
%      o     T2   ---       q2 = k2(T1-T2)/dz_2 - k3(T2-T3)/dz_3 + Sc2*dz
%      |           |
% k3 .....        dz_3      interface 3: k3 = 2*k(T3)*k(T2) / (k(T3)+k(T2))
%      |           |
%      o     T3   ---       q3 = k3(T2-T3)/dz_3 - k4(T3-T4)/dz_4 + Sc3*dz
%      |           |
% k4 .....        dz_4      interface 4: k4 = 2*k(T4)*k(T3) / (k(T4)+k(T3))
%      |           |
%      o     T4   ---
%      |           |
% k5 .....        dz_5
%      |           |
%      o     T5   ---
%      |           |
% k6 .....        dz_6 (N+1)
%  ----o---- T6 (N+1)      lower boundary dT/dz = 0, dz_pbc = 0, delz = 0
%  /////////
%
% Note that dz for the c.v.'s is constant but dz for the heat flux terms
% includes a 1/2 c.v. at the top and bottom (dz_1 and dz_N+1)
% with the upper and lower boundaries included the actual arrays will have
% an additional level for the upper boundary

%% More notes
% This is the c.v. geometry, labeled as in Patankar, Fig. 4.3. Note that
% unlike Patankar, i is defined at the surface, not at 1/2 c.v. width., and
% I use N for number of internal grid points / control volumes (i.e. the
% upper boundary point B is not included in N)
% ---------------------------------------------
% /|:       :       :       :       :       :|\
% /o:---o---:---o---:---o---:---o---:---o---:o\
% /|:       :       :       :       :       :|\
% B,i   I   w   W   p   P   e   E       N   N+1
% ---------------------------------------------

%  0  .15     .45     .75     1.05    1.35  1.5  z_node_bc (nodes with boundary nodes, N+2)

%  0       .3      .6      .9      1.2      1.5  z_edge (interfaces, N+1)

%|0|--.3---|--.3---|--.3---|--.3---|--.3---|0|  dz (c.v. widths, used for internal temperature, N, or N+2 including fictitious boundary volumes)

%|0|.15-|--.3---|--.3---|--.3---|--.3---|-.15|  delz  (point spacing, used for heat fluxes, N+1, or N+2 including fictitious upper boundary external interface)

%|0|.15-|   |.15|   |.15|   |.15|   |.15|  |-|  delz_pos (interface-to-next point spacing N+1 or N+2 as in delz)

%|0| 1  |  0.5  |  0.5  |  0.5  |  0.5  |  - |  f_n =  (weighting factor, defined at the interfaces, used for heat fluxes, N+1 or N+2 as in delz)
%|0|.15-|.15|   |.15|   |.15|   |.15|   |-|     delz_pos /
%|0|.15-|--.3---|--.3---|--.3---|--.3---|-.15|  delz

% with this arrangement, lets look at the boundary flux:
% we need the boundary conductivity, which we get from the interface
% conductivity equation 4.9:
% ki = 1/[(1-f(i) / kB) + (f(w) / kI)]
% ki = 1/[(1-1 / kB) + ( 0.5/ kI)]
% ki = 1/[0 + ( 0.5/ kI)] = kI, so ki = kI which is good because kB is undefined
% the flux at the top control volume is then:
% qI = kw(TW-TI) - ki(TI-TB) + ScI
% this says the flux through c.v. I is the flux into I at interface w plus
% the source term defined at I minus the flux out of I to the surface
% also note that kw = 2*kI*kW/(kI+kW) which is also what we want - the
% interface conductivity is the harmonic mean of the two adjacent
% conductivity values, which are defined at each c.v. center

% P = grid point
% E = east side grid point, in the positive x-direction from P
% W = west side grid point, in the negative x-direction from P
% e, w = control volume faces (dashed lines in figure 3.2, aka interfaces)
% the source term is integrated from w->e, the heat flux is computed at the
% two faces e.g. (k*dT/dx)e - (k*dT/dx)w + int(S*dx)_w^e = 0
% T is interpreted as a point value at each P and is interpolated between
% points (piecewise-linear profile)

% Patankar on left, Glen on right
% e = n, E = N
% w = s, W = S
% P = j, or k, the point index
% (delx)w or (delx)e = dely_p
