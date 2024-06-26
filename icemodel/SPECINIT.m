function [radii,mie,solar,kabs,kice] = SPECINIT(opts)
   %SPECINIT initialize the spectral coefficients
   %
   %#codegen

   % The grain radii that can be used for the two-stream spectral model
   radii = [ 0.040, 0.050, 0.065, 0.080, 0.100, ...
      0.120, 0.140, 0.170, 0.200, 0.240, ...
      0.290, 0.350, 0.420, 0.500, 0.570, ...
      0.660, 0.760, 0.870, 1.000, 1.100, ...
      1.250, 1.400, 1.600, 1.800, 2.000, ...
      2.250, 2.500, 2.750, 3.000, 3.500, ...
      4.000, 4.500, 5.000, 5.500, 6.000 ];

   % load the pre-defined Mie scattering values.
   mie = load(fullfile(opts.pathinput,'spectral','mie.mat'),'mie');
   mie = mie.mie;

   % load the proto-typical spectral irradiance profile
   solar = load(fullfile(opts.pathinput,'spectral','solar.mat'),'solar');
   solar = solar.solar;

   % load the user-defined kabs/kice if provided
   if opts.kabs_user == true
      kabs = load(fullfile(opts.pathinput,'spectral','kabs.mat'),'kabs');
      kice = load(fullfile(opts.pathinput,'spectral','kice.mat'),'kice');
      kabs = kabs.kabs;
      kice = kice.kice;
   else
      kabs = []; kice = [];
   end
end
