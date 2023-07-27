function [dz, delz, f_n, gridz] = CVTHERMAL(Z,dz)
%CVTHERMAL build control volume used for heat transfer
% 
% JJ = internal temperature nodes / c.v.'s (5 in example below
% JJ+1 = interfaces (6 in example below)
% 
% conductivity (gamma) is a function of T and is defined at T1, T2, ..., TN
% k is the interface conductivity and is calculated from the values at T1,
% T2, ..., TN using eq. 4.9. This puts k on the c.v. interfaces: 
% 
%  /////////
%  ----o---- Tsfc ---       (upper boundary T = Tsfc, dy_pbc = 0, dely = 0) 
% k1 .....         |        interface 1: k1 = k(T1) 
%      |          dy_1
%      o     T1   ---       q1 = k1(Tsfc-T1)/dy_1 - k2(T1-T2)/dy_2 + Sc1*dy
%      |           |
% k2 .....        dy_2      interface 2: k2 = 2*k(T1)*k(T2) / (k(T1)+k(T2))
%      |           |
%      o     T2   ---       q2 = k2(T1-T2)/dy_2 - k3(T2-T3)/dy_3 + Sc2*dy
%      |           |
% k3 .....        dy_3      interface 3: k3 = 2*k(T3)*k(T2) / (k(T3)+k(T2))
%      |           |
%      o     T3   ---       q3 = k3(T2-T3)/dy_3 - k4(T3-T4)/dy_4 + Sc3*dy
%      |           |
% k4 .....        dy_4      interface 4: k4 = 2*k(T4)*k(T3) / (k(T4)+k(T3))
%      |           |
%      o     T4   ---
%      |           |
% k5 .....        dy_5
%      |           |
%      o     T5   ---
%      |           |
% k6 .....        dy_6 (JJ+1)
%  ----o---- T6 (JJ+1)      lower boundary dT/dz = 0, dy_pbc = 0, dely = 0
%  /////////
% 
% Note that dy for the c.v.'s is constant but dy for the heat flux terms
% includes a 1/2 c.v. at the top and bottom (dy_1 and dy_JJ+1)
% with the upper and lower boundaries included the actual arrays will have
% an additional level for the upper boundary:

% activate JJ and deltaz to produce the geometry above:
% JJ          =   5;                    % number of control volumes
% deltaz      =   0.3;                  % width of each control volume
JJ          =   Z/dz;                   % number of nodes           [#]
deltab      =   0.0;                    % width of c.v. at the boundaries
dz          =   dz.*ones(JJ,1);         % array of c.v. widths
dz_pbc      =   [deltab;dz;deltab];     % array of c.v. widths including boundaries
delz_neg    =   0.5.*dz_pbc(1:JJ+1);    % interface-to-previous point
delz_pos    =   0.5.*dz_pbc(2:JJ+2);    % interface-to-next point
delz        =   delz_neg + delz_pos;    % distance between grid points
z_crds      =   [0;cumsum(delz)];       % grid point coordinates       
z_wall      =   [0;cumsum(dz)];         % interface coordinates
f_n         =   delz_pos./delz;         % interface conductivity weighting factor

gridz       =   z_crds(2:end-1);        % just the grid points

% round the outputs to the nearest milimeter (except f_n)       
dz          =   round(dz,3);
dz_pbc      =   round(dz_pbc,3);
delz        =   round(delz,3);
f_n         =   round(f_n,1);
z_crds      =   round(z_crds,3);
z_wall      =   round(z_wall,3);
gridz       =   round(gridz,3);

% See notes at the end for further clarification.

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

%  0  .15     .45     .75     1.05    1.35  1.5  y_crds (grid points, JJ+2)

%  0       .3      .6      .9      1.2      1.5  y_wall (interfaces, JJ+1)

%|0|--.3---|--.3---|--.3---|--.3---|--.3---|0|  dy_pbc  (c.v. widths, used for internal temperature, JJ, or JJ+2 including non-existent boundary volumes)

%|0|.15-|--.3---|--.3---|--.3---|--.3---|-.15|  dely_p  (point spacing, used for heat fluxes, JJ+1, or JJ+2 including non-existent upper boundary external interface)

%|0|.15-|   |.15|   |.15|   |.15|   |.15|  |-|  dely_pos (interface-to-next point spacing JJ+1 or JJ+2 as in dely_p)

%|0| 1  |  0.5  |  0.5  |  0.5  |  0.5  |  - |  f_n =  (weighting factor, defined at the interfaces, used for heat fluxes, JJ+1 or JJ+2 as in dely_p)
%|0|.15-|.15|   |.15|   |.15|   |.15|   |-|     dely_pos /
%|0|.15-|--.3---|--.3---|--.3---|--.3---|-.15|  dely_p 

% with this arrangement, lets look at the boundary flux:
% we need the boundary conductivity, which we get from the interface
% conductivity equation 4.9:
% ki = 1/[(1-f(i) / kB) + (f(w) / kI)]
% ki = 1/[(1-1 / kB) + ( 0.5/ kI)]
% ki = 1/[0 + ( 0.5/ kI)] = kI, so ki = kI which is good because kB is
% undefined
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
