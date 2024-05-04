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
   %
   %#codegen

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
   mie = load(fullfile(opts.pathinput, 'spectral', 'mie.mat'), 'mie');
   mie = mie.mie;

   % load the proto-typical spectral irradiance profile
   solar = load(fullfile(opts.pathinput, 'spectral', 'solar.mat'), 'solar');
   solar = solar.solar;

   % load the user-defined kabs/kice if provided
   if opts.kabs_user == true
      kabs = load(fullfile(opts.pathinput, 'spectral', 'kabs.mat'), 'kabs');
      kice = load(fullfile(opts.pathinput, 'spectral', 'kice.mat'), 'kice');
      kabs = kabs.kabs;
      kice = kice.kice;
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

% Notes
% (1) sum(netflux) = (1-albedo)*total_solar

% but sum(netflux) for the numerical calculation should be:
% (2) sum(netflux) = (1-albedo)*Qsi

% so multiply (1) by Qsi/total_solar:
% (3) (1-albedo)*total_solar * Qsi/total_solar = (1-albedo)*Qsi

% or as its expressed in the numerical setup:
% (4) Sc = (Qsi/total_solar .* netflux) ./dy_p;

% if total_solar was the actual incoming solar at each timestep, it would be:
% Sc = dq/dz = netflux./dy_p;

% To assign a portion of the incoming radiation, qsfactor, to surface heating,
% and the rest to subsurface layers, both xynet and Sc need to be adjusted for
% the surface portion.
