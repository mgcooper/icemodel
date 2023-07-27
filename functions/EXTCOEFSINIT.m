function [  total_solar,                                                ...
            grid_spect,                                                 ...
            z_walls,                                                    ...
            spect_lower,                                                ...
            spect_upper,                                                ...
            solardwavl ]   =  EXTCOEFSINIT(opts,radii,scattercoefs,     ...
                              solar,kabs,kice,dz_spect,JJ_spect,ro_ice);
%EXTCOEFSINIT initialize the extinction coefficients

% Note that the outputs correspond to the center of each level,
%   and the top surface of the top grid cell, and the bottom
%   surface of the bottom grid cell (thus the 502 values, for
%   the current setup).

% There should be no limit to how coarse (like 10 cm) or how fine
%   (like 1 mm) you want to define your levels to be.  There is
%   an upper limit of 500 levels hard coded in the program that
%   could be changed if you want/need to.

% The simulation domain is defined by using a combination of deltaz
%   (the thickness of each level) and nz (the number of levels).

% Note that: down = solar down, up = solar up, up/down = albedo,
%   and up-down = net solar (or - solar_absorbed).

    nwavl = opts.nwavl;

% build a control volume for the spectral model. note that deltaz_spectral
% is an inital value that defines the 1 mm upper node(s) spacing, where
% dz_spectral is the actual c.v. thicknesses, adjusted for
% quasi-exponential spacing
[   ~,                                                                  ...
    ~,                                                                  ...
    ~,                                                                  ...
    ~,                                                                  ...
    ~,                                                                  ...
    z_walls,                                                            ...
    grid_spect ]  =   CVSPECTRAL(JJ_spect,dz_spect);

% ice/snow grain radii
   r_snow         =   radii(opts.i_grainradius) / 1000.0;

% Read in the wavelength-dependent scattering coefficient arrays.
[  g,                                                                   ...
   qext,                                                                ...
   ss_coalb,                                                            ...
   wavelength ]   =   GETSCATTERCOEFS(opts,scattercoefs);

% Generate a delta_wavelength array.
	dwavl          =   GETDWAVL(wavelength,nwavl);

% Produce a downward (spectral) solar spectrum and integrated value. Note
% that the input 'solar' is interpolated here to the 118 bands
[  solar,                                                              ...
   total_solar ]  =   GETSOLAR(solar,nwavl,wavelength,dwavl);

% Compute the spectral extinction coefficients as a function of wavelength.
	spect_coefs    =   SPECTEXTCOEF(opts,qext,g,ss_coalb,r_snow);

% Scale the extinction coefficients by a user-defined absorption
%    coefficient profile to account for impurities in the ice/snow
if opts.kabs_user == true
   spect_coefs    =   SCALESPECTEXTCOEF(spect_coefs,wavelength,kice,kabs);
end

% send these into updateextcoefs to improve speed of exponentiation
   solardwavl     =  solar.*dwavl; 
   spect_walls    =  -z_walls*spect_coefs/ro_ice;
   % spect_walls    =  -z_walls.*repmat(spect_coefs./ro_ice,JJ_spect+1,1);
   spect_lower    =  spect_walls(2:end,:);
   spect_upper    =  spect_walls(1:end-1,:);
   
% % tried using exp here but it's slower
%    spect_walls    =  exp(-z_walls.*repmat(spect_coefs./ro_ice,1,2001));
%    spect_lower    =  spect_walls(:,2:end);
%    spect_upper    =  spect_walls(:,1:end-1);

                         
%--------------------------------------------------------------------------
% This is here for clarification. The xynet above is the column-integrated
% value at each c.v. interface, e.g. xynet(1) is the the column-integrated
% absorbed flux from z(\inf) to z(0), and xynet(2) would be from z(\inf) to
% z(1), where z(1) is the bottom of the top c.v. For me, it is more
% intuitive to see the amount that was absorbed within each layer:

% Calculate the net down, up, and absorbed flux WITHIN each c.v.
%     netdown         =   down(1:nz_spectral+1) - down(2:nz_spectral+2);
%     netup           =   up(1:nz_spectral+1) - up(2:nz_spectral+2);
%     netflux         =   netdown - netup;
%     
%     sum(netflux)
% netdown and netup are what go in and out of each c.v.
% what goes in minus what comes out is what got absorbed
% the sum of what got absorbed should equal (1-albedo)*total_solar
%     sum(netflux)
%     (1-albedo)*total_solar
% 
% % note that netflux can also be expressed as:
%     netflux         =   xynet(1:JJ)-xynet(2:JJ+1);
% 
% % try with the higher resolution spectral grid
%     JJ2 = length(xynet_old)
%     
%     netflux         =   xynet_old(1:JJ2-1)-xynet_old(2:JJ2);
%     
% in HEATSOLVE, this amount, sum(netflux), is scaled by Qsip/total_solar so
% that sum(netflux) equals Qsip:

% consider:
% (1) sum(netflux) = (1-albedo)*total_solar

% but we want sum(netflux) that goes into the numerical calculation to be:
% (2) sum(netflux) = (1-albedo)*Qsi 

% so we multiply (1) by Qsi/total_solar:
% (3) (1-albedo)*total_solar * Qsi/total_solar = (1-albedo)*Qsi

% or as its expressed in the numerical setup:
% (4) Sc = (Qsi/total_solar .* netflux) ./dy_p;

% if total_solar was the actual incoming solar at each timestep, it would
% be simple:
% Sc = dq/dz = netflux./dy_p;


% Now let's say I wanted to let a portion, qsfactor, of the incoming
% radiation be assigned to surface heating, and the rest penetrate. I would
% need to account for that within xynet as well as Sc. 
