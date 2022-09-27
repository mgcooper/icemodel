%--------------------------------------------------------------------------
%   Update the bulk extinction coefficient based on albedo and ice density
%--------------------------------------------------------------------------
function [  dQp,                                                        ...
            chi  ]   =  UPDATEEXTCOEFS(grid_spect,JJ_spect,dz_spect,    ...
                        grid_therm,JJ_therm,dz_therm,z_walls,ro_sno,    ...
                        total_solar,spect_lower,spect_upper,albedo,     ...
                        solardwavl,qsfactor,snowd,opts)
%--------------------------------------------------------------------------

% for reference, this is the structure when i used nested functions

if opts.skinmodel == true
    dQp     =  zeros(size(grid_therm));
    chi     =  1.0;
    return;
end
    
% Transform the mass density to the spectral grid resolution
   ro_sno   = interp1(grid_therm,ro_sno,grid_spect,'nearest','extrap');

% Compute the downward bulk extinction coefficient.
%    bulkcoefs = BULKEXTCOEF(dz_spect,ro_sno,spect_lower,spect_upper,solardwavl);

   bulkcoefs = -log((sum(exp(spect_lower.*ro_sno).*solardwavl,2))./     ...
               (sum(exp(spect_upper.*fix(max(ro_sno,300))).*            ...
               solardwavl,2)))./dz_spect;

% Cast in a form that can be used in the two-stream computation (add the
%   boundaries).  Here I have assumed that it is okay to call the
%   boundaries equal to the value at the center of that grid cell (the
%   value prevails thoughout the cell).
   bulkcoefs = vertcat(bulkcoefs,bulkcoefs(end),bulkcoefs(end));
   
% Compute the a and r coefficients from knowledge of the
%   surface albedo and the extinction coefficient.
% 	[a,r]    =  GETAANDR(bulkcoefs,albedo);
   a        =  (1.0 - albedo) ./ (1.0 + albedo) .* bulkcoefs;
   r        =  2.0 .* albedo .* bulkcoefs ./ (1.0 - albedo^2);

% Solve the system of equations for the two-stream model
% 	xynet    =  SOLVETWOSTREAM(a,r,bulkcoefs,total_solar,albedo,z_walls);
   
% notation here roughly follows Schlatter
   M        =  JJ_spect;
   
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
   e(2:M+1)    =   1.0 + (r(3:M+2) - r(1:M)) ./ (4.0.*r(2:M+1));
   f(2:M+1)    =   deltaz./(2.0.*r(2:M+1)).*(a(2:M+1).*(r(3:M+2)-r(1:M))-...
                   r(2:M+1).*(a(3:M+2)-a(1:M)))- ...
                   (2.0+deltaz.^2.*bulkcoefs(2:M+1).^2);
   g(2:M+1)    =   1.0 - (r(3:M+2) - r(1:M)) ./ (4.0*r(2:M+1));
   b(2:M+1)    =   0.0;
   
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


% xynet is the penetrating radiation at each level. In Schlatter it is 
%   just Q. Below, I call the absorbed flux dQp. In Glen's original model,
%   this quantity was not defined, xynet was converted directly to dQ/dz in
%   ICE_HEAT by scaling d(xynet) by Qsip and then to Sc. The only reason it
%   is defined here is because I find it is easier to transform between the
%   spectral grid and the thermal grid in terms of dQ than xynet. 

%   Transform the pentrating radiation back to the thermal grid as dQp
   dQp      =  GRIDINVERSE(xynet,dz_spect,dz_therm,JJ_spect,JJ_therm);

% Compute the amount of solar radiation absorbed in the top grid cell, to
%   be allocated to the SEB. Optionally enhance it by qsfactor to account
%   for impurities on the ice surface or other factors. 
    if albedo > 0.65 % && snowd>0.05
        chi     =   0.9;
    else
        chi     =   dQp(1)/sum(dQp);
        chi     =   min(1.0,max(chi+qsfactor,0.0));
    end
