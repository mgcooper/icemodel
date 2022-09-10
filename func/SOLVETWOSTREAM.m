%--------------------------------------------------------------------------
%   Solve the two-stream model following Schlatter, 1972 method
%--------------------------------------------------------------------------

function xynet = SOLVETWOSTREAM(a,r,bulkcoefs,total_solar,albedo,z_walls)
%--------------------------------------------------------------------------

% notation here roughly follows Schlatter
   M        =  length(z_walls)-1;
   xmu      =  bulkcoefs;
   
%  e = A_sub    = sub diagonal
%  f = A_main   = main diagonal
%  g = A_super  = super diagonal
%  b = b_vector = right hand side
%  x = rad      = solution
   
% extend y_wall downward by one c.v.
   dz_bottom   =   z_walls(M+1) - z_walls(M);
   z_walls(M+2)=   z_walls(M+1) + dz_bottom;
   
% since bulk_extcoefs were computed with y_wall, do the same here
   deltaz      =   z_walls(2) - z_walls(1);
   
% initialize the matrix
   e           =   zeros(M+1,1);
   f           =   zeros(M+1,1);
   g           =   zeros(M+1,1);
   b           =   zeros(M+1,1);
   
% Account for the upper boundary condition
   alfa        =   1.0/(a(1)+r(1));
   e(1)        =   0.0;
   f(1)        =   1.0;
   g(1)        =   -alfa/(deltaz+alfa);
   b(1)        =   r(1)*total_solar*deltaz * alfa/(deltaz+alfa);
   
% Fill the vectors between the boundaries
   deltaz      =   z_walls(3:M+2)-z_walls(2:M+1);
%    tmp1        =   deltaz ./ (2.0 .* r(2:M+1));
%    tmp2        =   a(2:M+1) .* (r(3:M+2) - r(1:M));
%    tmp3        =   r(2:M+1) .* (a(3:M+2) - a(1:M));
%    tmp4        =   (2.0 + deltaz.^2 .* xmu(2:M+1).^2);
   e(2:M+1)    =   1.0 + (r(3:M+2) - r(1:M)) ./ (4.0.*r(2:M+1));
%    f(2:M+1)    =   tmp1 .* (tmp2 - tmp3) - tmp4;
   g(2:M+1)    =   1.0 - (r(3:M+2) - r(1:M)) ./ (4.0*r(2:M+1));
   b(2:M+1)    =   0.0;
   
% see if this speeds it up
   f(2:M+1)    =  deltaz./(2.0.*r(2:M+1)).*(a(2:M+1).*(r(3:M+2)-r(1:M))-...
                  r(2:M+1).*(a(3:M+2)-a(1:M)))-(2.0+deltaz.^2.*xmu(2:M+1).^2);
   
% Account for the lower boundary condition
   g(M+1)      =  0.0;
   b(M+1)      =  0.0;
   
% Solve the equation
   x           =  TRISOLVE(e,f,g,b);
   
% Add the boundary conditions to up and reconstruct down.
   [up,down]   =  GETUPDOWN(a,r,x,total_solar,z_walls,M);

% Build an array of source-term values on the c.v. boundaries.
   xynet       =  (up(2:M)+up(3:M+1))./2.0-(down(2:M)+down(3:M+1))./2.0;
   xynet       =  vertcat(up(1)-down(1), xynet);

% Ensure xynet(1) is equal to the total absorbed radiation. This corrects
%   for any (very) small errors in the two-stream model, typically due to a
%   spectral grid that is too shallow to absorb all the radiation
   xynet(1)      =  min(-total_solar*(1-albedo),xynet(1));
