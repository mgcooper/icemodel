%--------------------------------------------------------------------------
% Get the up/down flux at each depth from the two-stream solution vector x
%--------------------------------------------------------------------------
function [up,down] = GETUPDOWN(a,r,x,total_solar,z_walls,M)
%--------------------------------------------------------------------------

% X = up, Y = down
    N       =   M+2;
% Add the boundary conditions to rad.
    up      =   vertcat(x(1:N-1),0);

% Reconstruct y.
    dz      =   z_walls(3:N)-z_walls(2:N-1);
    down    =   (a(2:N-1)+r(2:N-1))./r(2:N-1).*up(2:N-1) -              ...
                        (up(3:N)-up(1:N-2)) ./ (2.0.*dz.*r(2:N-1));
    down    =   vertcat(total_solar,down);
    down    =   vertcat(down,(a(N)+r(N))/r(N)*up(N)-(up(N)-up(N-1))/(dz(1)*r(N)));

% Smooth any small bumps in the up and down curve.  This will assist
%   in any delta up, delta down computations. Do the interior.
    downtmp =   (down(2:N-1) + 0.5 .* (down(1:N-2) + down(3:N)))./ 2.0;
    uptmp   =   (up(2:N-1)   + 0.5 .* (up(1:N-2)   + up(3:N)))  ./ 2.0;

% Do the ends.
    downtmp =   vertcat((down(2) + 0.5 * down(1)) / 1.5, downtmp);
    uptmp   =   vertcat((up(2)   + 0.5 * up(1))   / 1.5, uptmp);
    downtmp =   vertcat(downtmp, (down(N-1) +  0.5 * down(N)) / 1.5);
    uptmp   =   vertcat(uptmp,   (up(N-1)   +  0.5 * up(N))   / 1.5);
    
% Rewrite the arrays. Adjust to get back the Qsi at the surface.
    down    =   vertcat(down(1), downtmp(2:N-1), down(end));
    up      =   vertcat(up(1),   uptmp(2:N-1),   up(end));

% Repeat the smoothing (deal with edge effects)
    downtmp =   (down(2:N-1) + 0.5 * (down(1:N-2) + down(3:N))) / 2.0;
    uptmp   =   (up(2:N-1) + 0.5 * (up(1:N-2) + up(3:N))) / 2.0;
% Do the ends.
    downtmp =   vertcat(down(2)   +  0.5 * down(1) / 1.5, downtmp);
    uptmp   =   vertcat(up(2)     +  0.5 * up(1)   / 1.5, uptmp);
    downtmp =   vertcat(downtmp, down(N-1) +  0.5 * down(N) / 1.5);
    uptmp   =   vertcat(uptmp,   up(N-1)   +  0.5 * up(N)   / 1.5);
    
% Rewrite the arrays.
    down    =   vertcat(down(1), downtmp(2:N-1), down(end));
    up      =   vertcat(up(1),   uptmp(2:N-1),   up(end));

% % Compute the net solar flux (don't call this from the main, just to test)
%     dz      =   z_walls(3:M-1) - z_walls(2:M-2);
%     xynet   =   (up(2:M-2)+up(3:M-1))./dz-(down(2:M-2)+down(3:M-1))./dz;
%     xynet   =   vertcat(up(1)-down(1), xynet, up(M) - down(M));