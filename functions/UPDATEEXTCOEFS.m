%--------------------------------------------------------------------------
%   Update the bulk extinction coefficient based on albedo and ice density
%--------------------------------------------------------------------------
function [  Sc,                                                         ...
            chi  ]   =  UPDATEEXTCOEFS(grid_spect,JJ_spect,dz_spect,    ...
                        grid_therm,JJ_therm,dz_therm,z_walls,dz,ro_sno, ...
                        Qsi,total_solar,spect_lower,spect_upper,albedo, ...
                        solardwavl)
%--------------------------------------------------------------------------

% Transform the mass density to the spectral grid resolution
   ro_sno   = interp1(grid_therm,ro_sno,grid_spect,'nearest','extrap');
   ro_sno   = max(ro_sno,300);
   
% Compute the downward bulk extinction coefficient.
   bulkcoefs = -log((sum(solardwavl.*exp(spect_lower.*ro_sno),2))./     ...
               (sum(solardwavl.*exp(spect_upper.*ro_sno),2)))./dz_spect;
%  NOTE: solardwavl comes multiplied by dwavl, and sum/trapz are identical 

% Cast in a form that can be used in the two-stream computation (add the
%   boundaries).  Here I have assumed that it is okay to call the
%   boundaries equal to the value at the center of that grid cell (the
%   value prevails thoughout the cell).
   bulkcoefs = vertcat(bulkcoefs,bulkcoefs(end),bulkcoefs(end));
   
% Compute the a and r coefficients from knowledge of the
%   surface albedo and the extinction coefficient.
   a  =  (1.0 - albedo) ./ (1.0 + albedo) .* bulkcoefs;
   r  =  2.0 .* albedo .* bulkcoefs ./ (1.0 - albedo^2);

% Solve the system of equations for the two-stream model
   M  =  JJ_spect;   % notation here roughly follows Schlatter
   M1 =  M+1;
   M2 =  M+2;
   
% extend y_wall downward by one c.v.
   z_walls(M2) =  z_walls(M1) + z_walls(M1) - z_walls(M);
   
% since bulk_extcoefs were computed with y_wall, do the same here
   deltaz   =  z_walls(2) - z_walls(1);
   
% initialize the matrix
   e  =  zeros(M1,1);
   f  =  zeros(M1,1);
   g  =  zeros(M1,1);
   b  =  zeros(M1,1);
   
% Account for the upper boundary condition
   alfa  =  1.0/(a(1)+r(1));
   e(1)  =  0.0;
   f(1)  =  1.0;
   g(1)  =  -alfa/(deltaz+alfa);
   b(1)  =  r(1)*total_solar*deltaz * alfa/(deltaz+alfa);
   
% Fill the vectors between the boundaries
   deltaz   =  z_walls(3:M2)-z_walls(2:M1);
   e(2:M1)  =  1.0 + (r(3:M2) - r(1:M)) ./ (4.0.*r(2:M1));
   f(2:M1)  =  deltaz./(2.0.*r(2:M1)).*(a(2:M1).*(r(3:M2)-r(1:M))-   ...
               r(2:M1).*(a(3:M2)-a(1:M)))-(2.0+deltaz.^2.*bulkcoefs(2:M1).^2);
   g(2:M1)  =  1.0 - (r(3:M2) - r(1:M)) ./ (4.0*r(2:M1));
   b(2:M1)  =  0.0;
   
% Account for the lower boundary condition
   g(M1)    =  0.0;
   b(M1)    =  0.0;
   
% Solve the equation
   x = nan(M1,1);
   for k = 2:numel(f)                        % forward elimination
      f(k)  =  f(k)-e(k)/f(k-1)*g(k-1);
      b(k)  =  b(k)-e(k)/f(k-1)*b(k-1);
   end
   x(M1) = b(M1)/f(M1);                      % back substitution
   for k = M:-1:1
      x(k)  =  (b(k)-g(k)*x(k+1))/f(k);
   end 
   
% Add the boundary conditions to up and reconstruct down. X = up, Y = down.
   up       =  vertcat(x(1:M1),0); % Add the boundary conditions to rad.

% Reconstruct y.
   down     =  (a(2:M1)+r(2:M1))./r(2:M1).*up(2:M1)-(up(3:M2)-up(1:M2-2))...
               ./(2.0.*deltaz.*r(2:M1));
   down     =  vertcat(total_solar,down);
   down     =  vertcat(down,(a(M2)+r(M2))/r(M2)*up(M2)-(up(M2)-up(M1))/...
               (deltaz(1)*r(M2)));

% Smooth any small bumps in the up and down curve.  This will assist
%   in any delta up, delta down computations. Do the interior.
   downtmp  =  (down(2:M1) + 0.5 .* (down(1:M2-2) + down(3:M2)))./ 2.0;
   uptmp    =  (up(2:M1)   + 0.5 .* (up(1:M2-2)   + up(3:M2)))  ./ 2.0;

% Do the ends.
   downtmp  =  vertcat((down(2) + 0.5 * down(1)) / 1.5, downtmp);
   uptmp    =  vertcat((up(2)   + 0.5 * up(1))   / 1.5, uptmp);
   downtmp  =  vertcat(downtmp, (down(M1) +  0.5 * down(M2)) / 1.5);
   uptmp    =  vertcat(uptmp,   (up(M1)   +  0.5 * up(M2))   / 1.5);

% Rewrite the arrays. Adjust to get back the Qsi at the surface.
   down     =  vertcat(down(1), downtmp(2:M1), down(end));
   up       =  vertcat(up(1),   uptmp(2:M1),   up(end));

% Repeat the smoothing (deal with edge effects)
   downtmp  =  (down(2:M1) + 0.5 * (down(1:M) + down(3:M2))) / 2.0;
   uptmp    =  (up(2:M1) + 0.5 * (up(1:M) + up(3:M2))) / 2.0;

% Do the ends.
   downtmp  =  vertcat(down(2)   +  0.5 * down(1) / 1.5, downtmp);
   uptmp    =  vertcat(up(2)     +  0.5 * up(1)   / 1.5, uptmp);
   downtmp  =  vertcat(downtmp, down(M1) +  0.5 * down(M2) / 1.5);
   uptmp    =  vertcat(uptmp,   up(M1)   +  0.5 * up(M2)   / 1.5);

% Rewrite the arrays.
   down     =  vertcat(down(1), downtmp(2:M1), down(end));
   up       =  vertcat(up(1),   uptmp(2:M1),   up(end));
   
% Build an array of source-term values on the c.v. boundaries.
   xynet    =  (up(2:M)+up(3:M1))./2.0-(down(2:M)+down(3:M1))./2.0;
   xynet    =  vertcat(up(1)-down(1), xynet);

% Ensure xynet(1) is equal to the total absorbed radiation. This corrects
%   for any (very) small errors in the two-stream model, typically due to a
%   spectral grid that is too shallow to absorb all the radiation
   xynet(1) =  min(-total_solar*(1-albedo),xynet(1));

% xynet is the net solar radiation at each level. In Schlatter it is
%   just Q. The source term (absorbed flux per unit volume) is dQ/dz.
%   Below, I call the absorbed flux dQ, so Sc = dQ/dz = d(xynet)/dz. In
%   Glen's original model, the quantity dQ was not defined, xynet was
%   converted directly to dQ/dz in ICE_HEAT by scaling d(xynet) by Qsip and
%   then to Sc.

%   Transform the pentrating radiation back to the thermal grid as dQ
   % extrapolate xynet to the bottom of the thermal grid
   M3             =  JJ_therm*dz_therm/dz_spect;
   M4             =  round(M3-M,0);
   dxy            =  xynet(M)-xynet(M-1);
   xyextrap       =  zeros(M4,1);
   xyextrap(1)    =  xynet(M) + dxy;
   xyextrap(2:M4) =  xyextrap(1:M4-1) + dxy;
   xynew          =  [xynet;xyextrap];
% upscale it to the thermal grid
   dQ             =  xynew(1:M3-1)-xynew(2:M3);
   dQ(M3)         =  dQ(M3-1);
% get the number of grid cells in the first 3 m segment on each grid
   rz             =  M3/JJ_therm;
% reshape the radiation into equal chunks along the thermal grid
   dQ             =  reshape(dQ,rz,JJ_therm);    
% sum up those equal chunks to get the absorbed radiation in each thermal c.v.
   dQ             =  transpose(sum(dQ,1));

% Compute the amount of solar radiation absorbed in the top grid cell, to
%   be allocated to the SEB. Optionally enhance it by qsfactor to account
%   for impurities on the ice surface or other factors. 
   if albedo>0.65 % && snowd>0.05
      chi   =  0.9;
   else
      chi   =  dQ(1)/sum(dQ);
%       chi   =  min(1.0,max(chi+qsfactor,0.0));
   end
%    chi      =  0;
   Sc       =  -(1.0-chi)*Qsi/total_solar.*dQ./dz;    % [W/m3]

%    if albedo>0.65 % && snowd>0.05
%       chi   =  0.9;
%    else
%       if zD<=dz(1)
%         %Qseb  =  (dQ(1)*Qsi/total_solar).*exp(-zD.*a(1)); 
%          Qseb  =  dQ(1).*exp(-zD.*a(1));
%       else
%          Qseb  =  dQ(1);
%       end
%       chi   =  Qseb/sum(dQ);
%       dQ(1) =  dQ(1)-Qseb;
%    end
%    % compute the solar radiation source-term values in each c.v.
%    Sc       =  -Qsi/total_solar.*dQ./dz;    % [W/m3]

% at this point, sum(dQ) - Qseb = -total_solar*(1-albedo)
   %Qseb    =  -Qsi/total_solar*Qseb;
% these should be equal: [chi,  Qseb/(sum(Sc.*dz)+Qseb)]
   

   
%    % this was me trying to adjust Sc after the fact
%    Qseb     =  dQ(1).*exp(-zD.*a(1));
%    Sc(1)    =  Sc(1) + Qseb/dz(1);
   
% %    % van dalum does it this way, here at the radiation timestep, get tau:
% %    Ftop     =  xynet(1);
% %    tau      =  dz./(log(xynet(1)./xynet));
% %    
% %    % then in the substepping, use tau to determine how much of xynet
% %    % diffuses back to the SEB:
% %    Fseb     =  xynet.*exp(dz./tau);
% 
%    % but I don't need to do it this way, b/c I already know the total
%    % absorbed radiation in each level for the 1-hr timestep, so I just need
%    % to know how much of that radiation diffuses to the surface, which is
%    % just the exponential decay equation evaluated at the diffusion length
%    % scale. The absorbed flux that contributes to the SEB is then:
%    zD       =  sqrt(k_eff.*dt_new./(ro_sno.*cp_sno));
%    Fseb     =  dQ.*exp(-zD.*a);
%    Sc(1)    =  Sc(1) - Fseb/dz(1);
%    % and Fseb get's sent into SFCTEMP (or SFCFLIN)
%    
%    % to implement, first check if zD exceeds the top layer thickness:
%    if sqrt(k_eff(1)*dt/(ro_sno(1)*cp_sno(1))) <= dz(1)
%       Fseb  =  dQ(1)*exp(-sqrt(k_eff(1)*dt/(ro_sno(1)*cp_sno(1)))*a(1));
%    else % assign all of the top layer source term to the SEB
%       Fseb  =  dQ(1);
%    end

% note: using Dalum's method, my first thought was that tau is determined
% using dQ: 
%     tau   =  dz./(log(dQ./(Qsi*(1-albedo))));
% then send this into the substepping to get Fseb at each sub-step:
%     Fseb  =  Qsi*(1-albedo)
% but this is wrong, 

%    % xynet 
%    tau      =  deltaz./(log(xynet./xynet(1)));
%    
%    % it should be the absorbed flux
%    dxy      =  diff(xynet);
%    tau      =  deltaz(1:end-1)./(log(dxy./dxy(1)));
%    
%    figure; plot(grid_spect(1:end-1),-tau); set(gca,'YDir','reverse')
%    
%    tau      =  (z_walls(2:end)-z_walls(1:end-1))./(log(down(2:end)./down(1)));
%    
%    
%    figure; loglog(z_walls(2:end),-tau); set(gca,'YDir','reverse')
%    
%    tau      =  grid_spect./(log(xynet./xynet(1)));
%    
%    figure; plot(grid_spect,tau); set(gca,'YDir','reverse')
%    
%    tau      =  grid_spect./(log(down./down(1)));
%    
%    figure; plot(grid_spect,xynet); set(gca,'YDir','reverse')
%    
%    figure; plot(grid_spect,tau); set(gca,'YDir','reverse')
   
% % these are some tests i was doing to compare sum to trapz    
% wavl = wavl(1,:);
% 
% f1 = solardwavl.*exp(spect_lower.*ro_sno);
% f2 = solardwavl.*exp(spect_upper.*ro_sno);
% 
% f3 = solar.*exp(spect_lower.*ro_sno);
% f4 = solar.*exp(spect_upper.*ro_sno);
% 
% 
% figure; plot(wavl,f1(1,:))
% figure; plot(wavl,f3(1,:))
% 
% test1 = trapz(wavl,f3,2);
% test2 = sum(f3.*dwavl,2);
% figure; myscatter(test1,test2); addOnetoOne;
% 
% test  = -log((trapz(wavl,f3,2))./(trapz(wavl,f4,2)))./dz_spect;
% 
% bulkcoefs = -log((trapz(solardwavl.*exp(spect_lower.*ro_sno),2))./   ...
%       (trapz(solardwavl.*exp(spect_upper.*fix(max(ro_sno,300))),2)))./dz_spect;
% 
% tic; test1 = trapz(solardwavl.*exp(spect_lower.*ro_sno),2); toc;
% tic; test2 = sum(solardwavl.*exp(spect_lower.*ro_sno),2); toc;
% tic; test3 = trapz(wavl,solar.*exp(spect_lower.*ro_sno),2); toc;
