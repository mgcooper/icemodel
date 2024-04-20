function [Q0, dz_spect, spect_N, spect_S, solardwavl] = EXTCOEFSINIT(opts, ro_ice)
   %EXTCOEFSINIT initialize the extinction coefficients
   %
   % The outputs correspond to the center of each level, and the top
   % surface  of the top grid cell, and the bottom surface of the
   % bottom grid cell  (thus the 502 values, for the current setup).
   %
   % There should be no limit to how coarse (like 10 cm) or how fine
   % (like 1 mm) you want to define your levels to be. There is an upper
   % limit of 500 levels hard coded in the program that could be changed
   % if you want/need to.
   %
   % The simulation domain is defined by using a combination of deltaz
   % (the thickness of each level) and nz (the number of levels).
   %
   % Note that: down = solar down, up = solar up, up/down = albedo,
   % and up-down = net solar (or - solar_absorbed).
   %
   % See also:

   % The grain radii that can be used for the two-stream spectral model
   radii = [ ...
      0.040, 0.050, 0.065, 0.080, 0.100, ...
      0.120, 0.140, 0.170, 0.200, 0.240, ...
      0.290, 0.350, 0.420, 0.500, 0.570, ...
      0.660, 0.760, 0.870, 1.000, 1.100, ...
      1.250, 1.400, 1.600, 1.800, 2.000, ...
      2.250, 2.500, 2.750, 3.000, 3.500, ...
      4.000, 4.500, 5.000, 5.500, 6.000 ];

   % load the pre-defined Mie scattering values.
   load(fullfile(opts.pathinput, 'spectral', 'mie.mat'), 'mie');

   % load the proto-typical spectral irradiance profile
   load(fullfile(opts.pathinput, 'spectral', 'solar.mat'), 'solar');

   % load the user-defined kabs/kice if provided
   if opts.kabs_user == true
      load(fullfile(opts.pathinput, 'spectral', 'kabs.mat'), 'kabs');
      load(fullfile(opts.pathinput, 'spectral', 'kice.mat'), 'kice');
   else
      kabs = []; kice = [];
   end

   nwavl = opts.nwavl;

   % Build a control volume for the spectral model.
   [dz_spect, ~, ~, z_walls] = CVMESH(opts.z0_spectral, opts.dz_spectral);

   % Retrieve the ice/snow grain radii
   r_snow = radii(opts.i_grainradius) / 1000.0;

   % Read in the wavelength-dependent scattering coefficient arrays.
   [g, qext, ss_coalb, wavelength] = GETSCATTERCOEFS(opts, mie);

   % Generate a delta_wavelength array.
   dwavl = GETDWAVL(wavelength, nwavl);

   % Produce a downward (spectral) solar spectrum and integrated value. Note
   % that the input 'solar' is interpolated here to the 118 bands
   [solar, Q0] = GETSOLAR(solar, nwavl, wavelength, dwavl);

   % Compute the spectral extinction coefficients as a function of wavelength.
   spect_coefs = SPECTEXTCOEF(opts, qext, g, ss_coalb, r_snow);

   % Scale the extinction coefficients by a user-defined absorption coefficient
   % profile to account for impurities in the ice/snow
   if opts.kabs_user == true
      spect_coefs = SCALESPECTEXTCOEF(spect_coefs, wavelength, kice, kabs);
   end

   % Send these into updateextcoefs to improve speed of exponentiation
   solardwavl = solar .* dwavl;
   spect_walls = -z_walls * spect_coefs / ro_ice;
   spect_N = spect_walls(1:end-1, :);
   spect_S = spect_walls(2:end, :);
end

% % tried using exp here but it's slower
% spect_walls = exp(-z_walls.*repmat(spect_coefs./ro_ice,1,2001));
% spect_lower = spect_walls(:,2:end); spect_upper = spect_walls(:,1:end-1);

%--------------------------------------------------------------------------
% This is here for clarification. The xynet above is the column-integrated value
% at each c.v. interface, e.g. xynet(1) is the the column-integrated absorbed
% flux from z(\inf) to z(0), and xynet(2) would be from z(\inf) to z(1), where
% z(1) is the bottom of the top c.v. For me, it is more intuitive to see the
% amount that was absorbed within each layer:

% Calculate the net down, up, and absorbed flux WITHIN each c.v.
%     netdown = down(1:nz_spectral+1) - down(2:nz_spectral+2);
%     netup = up(1:nz_spectral+1) - up(2:nz_spectral+2);
%     netflux = netdown - netup;
%     sum(netflux)
%
% netdown and netup are what go in and out of each c.v. what goes in minus what
% comes out is what got absorbed the sum of what got absorbed should equal
% (1-albedo)*total_solar
%     sum(netflux) (1-albedo)*total_solar
%
% % note that netflux can also be expressed as:
%     netflux = xynet(1:JJ)-xynet(2:JJ+1);
%
% % try with the higher resolution spectral grid
%     JJ2 = length(xynet_old)
%
%     netflux = xynet_old(1:JJ2-1)-xynet_old(2:JJ2);
%
% in HEATSOLVE, this amount, sum(netflux), is scaled by Qsip/total_solar so that
% sum(netflux) equals Qsip:

% consider:
% (1) sum(netflux) = (1-albedo)*total_solar

% but we want sum(netflux) that goes into the numerical calculation to be:
% (2) sum(netflux) = (1-albedo)*Qsi

% so we multiply (1) by Qsi/total_solar:
% (3) (1-albedo)*total_solar * Qsi/total_solar = (1-albedo)*Qsi

% or as its expressed in the numerical setup:
% (4) Sc = (Qsi/total_solar .* netflux) ./dy_p;

% if total_solar was the actual incoming solar at each timestep, it would be:
% Sc = dq/dz = netflux./dy_p;

% Now let's say I wanted to let a portion, qsfactor, of the incoming radiation
% be assigned to surface heating, and the rest penetrate. I would need to
% account for that within xynet as well as Sc.
