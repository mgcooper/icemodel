clean
%==========================================================================

% Dec 2020 this was named 'build_mie_library' suggesting it was the final
% one I used but in fact this seems to just be a comparison b/w old and new
% mie vs Mie function from matzler so I renamed it. I am nearly certain the
% sequence was: build_big_mie_library and then build_mie_08_v3

use_84      =   false;
save_data   =   true;

%==========================================================================
%% set paths
%==========================================================================
path.data   =   'GREENLAND/runoff/icemodel/data/preprocess/spectral/';
path.save   =   'GREENLAND/runoff/icemodel/data/preprocess/spectral/';
path        =   setpath(path);

%==========================================================================
%% read in the complex index of refraction of pure ice
%==========================================================================
if use_84 == true
    m       =   readtable([path.data 'm_warren_84.xlsx']);
else
    m       =   readtable([path.data 'm_warren.xlsx']);
end
%==========================================================================
%% extract the values for the solar spectrum
%==========================================================================
si          =   find(m.lambda == 0.25);
ei          =   find(m.lambda == 3.003);
m_im        =   m.m_imag(si:ei);
m_re        =   m.m_real(si:ei);
m_lambda    =   m.lambda(si:ei); clear m;
m_lambda    =   m_lambda ./ 1e6;    % convert wavelength from [um] to [m]           

% these are the radii I am using, which range from very fine 
% grained snow (0.04 mm = 40 um) to very large grained blue ice (6 mm)
radii       =   [   0.040, 0.050, 0.065, 0.080, 0.100, ...
                    0.120, 0.140, 0.170, 0.200, 0.240, ...
                    0.290, 0.350, 0.420, 0.500, 0.570, ...
                    0.660, 0.760, 0.870, 1.000, 1.100, ...
                    1.250, 1.400, 1.600, 1.800, 2.000, ...
                    2.250, 2.500, 2.750, 3.000, 3.500, ...
                    4.000, 4.500, 5.000, 5.500, 6.000 ];

nradii      =   length(radii);
radii       =   radii ./ 1000;  % convert to meters

%==========================================================================
%% Calculate mie parameters
%==========================================================================
for j = 1:nradii

    r_i     =   radii(j);       % Grain radius [m]
    N       =   3/4/(pi*r_i^3); % Number of grains per unit volume [m-3]
    k_num   =   2*pi./m_lambda; % wavenumber [m-1];
    
  % Single scattering properties g, w, and Qext for each wavelength
    for n = 1:length(m_lambda)

        k_n     =   k_num(n);   % wavenumber 
        mi      =   1i*m_im(n); % imaginary 
        mr      =   m_re(n);    % real
        m_n     =   mr + mi;
        
        % mie calculations
        mie_pars    =   Mie(m_n,k_n*r_i);
        qext(j,n)   =   mie_pars.Qext;          % extinction efficiency
        omeg(j,n)   =   mie_pars.Omega;         % single-scattering albedo
        asym(j,n)   =   mie_pars.Asy;           % assymetry factor
        sige(j,n)   =   N*pi*(r_i^2)*qext(j,n); % radiance extinction coefficient
        
        % activate this to compare the old 'mie' function to new 'Mie'
%         mie_pars    =   mie(m_n,k_n*r_i); 
%         qext(j,n)   =   mie_pars.Qext;          % extinction efficiency
%         omeg(j,n)   =   mie_pars.Omega;         % single-scattering albedo
%         asym(j,n)   =   mie_pars.Asy;           % assymetry factor
%         sige(j,n)   =   N*pi*(r_i^2)*qext(j,n); % radiance extinction coefficient
    end
end

%==========================================================================
%% plot old mie vs new Mie
%==========================================================================
% Qext
figure; 
plot(m_lambda,qext(j,:)); hold on;
plot(m_lambda,Qext2(j,:),':'); 
legend('v1','v2');
ylabel('Qext')

figure; 
plot(m_lambda,qext(j,:)-Qext2(j,:)); hold on; 
legend('v1-v2');
ylabel('Qext')

% omega
figure; 
plot(m_lambda,omeg(j,:)); hold on;
plot(m_lambda,w2(j,:),':'); 
legend('v1','v2');
ylabel('w')

figure; 
plot(m_lambda,omeg(j,:)-w2(j,:)); hold on; 
legend('v1-v2');
ylabel('w')

% asy
figure; 
plot(m_lambda,asym(j,:)); hold on;
plot(m_lambda,g2(j,:),':'); 
legend('v1','v2');
ylabel('g')

figure; 
plot(m_lambda,asym(j,:)-g2(j,:)); hold on; 
legend('v1-v2');
ylabel('g')

%==========================================================================
%% read in the 'mie.dat' data
%==========================================================================
mie_data            =   importdata([path.data 'mie.dat']);
nvalues             =   size(mie_data,2); % number of wavelength bands
nclasses            =   size(mie_data,1)/4; % number of grain radii
lambda              =   mie_data(142,:);

% these are the radii Glen originally included
% radii = [ 0.005, 0.007, 0.010, 0.015, 0.020, ...
%           0.030, 0.040, 0.050, 0.065, 0.080, ...
%           0.100, 0.120, 0.140, 0.170, 0.200, ...
%           0.240, 0.290, 0.350, 0.420, 0.500, ...
%           0.570, 0.660, 0.760, 0.870, 1.000, ...
%           1.100, 1.250, 1.400, 1.600, 1.800, ...
%           2.000, 2.500, 3.000, 3.500, 4.000, ...
%           4.500, 5.000, 5.500, 6.000, 6.500, ...
%           7.000, 7.500, 8.000, 8.500, 9.000, ...
%           9.500,10.000];
      
% each coefficient is a 47x118 vector, grain size x wavelength, where grain
% size ranges from 0.005 - 10.00 mm, lambda varies from 0.299 um to 2.999
% um with 
% g                   =   mie_data(1:47,:);
% Qext                =   mie_data(48:94,:);
% w                   =   mie_data(95:141,:);
% lambda              =   mie_data(142,:);
% dlambda             =   diff(lambda);

% these are the radii I am using, with more resolution at mid-range values
radii = [ 0.005, 0.007, 0.010, 0.015, 0.020, ...
          0.030, 0.040, 0.050, 0.065, 0.080, ...
          0.100, 0.120, 0.140, 0.170, 0.200, ...
          0.240, 0.290, 0.350, 0.420, 0.500, ...
          0.570, 0.660, 0.760, 0.870, 1.000, ...
          1.100, 1.250, 1.400, 1.600, 1.800, ...
          2.000, 2.250, 2.500, 2.750, 3.000, ...
          3.500, 4.000, 4.500, 5.000, 5.500, ...
          6.000, 7.000, 7.500, 8.000, 8.500, ...
          9.000, 10.000];
% The number of grain radii that can be used.

%==========================================================================
%%
%==========================================================================

load([path.data 'm_warren_full_spectrum.mat']);
si                  =   find(m.lambda == 299);
ei                  =   find(m.lambda == 3001);
lambda              =   m.lambda(si:ei);
mimag               =   m.m_img(si:ei);
mreal               =   m.m_real(si:ei);

%% Gardner recommends mu = r, sigma = 0.15*r, N = 1000
f = 'ice_wrn_0011.nc';
info = ncinfo(f);
wvl = ncread(f,'wvl');
rds = ncread(f,'rds');
bnd_nbr = ncread(f,'bnd_nbr');

% single scattering coalbedo
ss_alb = ncread(f,'ss_alb');
% asymmetry parameter
asym = ncread(f,'asm_prm');
% extinction coefficint
k = ncread(f,'ext_cff_mss');

plot(wvl(1:100),k(1:100))
plot(wvl,ss_alb)
plot(wvl,asym)

%% This could be used to interpolate to a finer/coarser grid
% interpolate the values to a slightly coarser grid at long wavelengths
% si                  =   find(m_lambda == 300);
% mi                  =   find(m_lambda == 900);
% ei                  =   find(m_lambda == 3003);
% 
% % interpolate to a 1 nm grid as per Steve's suggestion
% lambda_new          =   300:1:4000;
% m_im_log            =   log(m_im);
% lambda_log          =   log(m_lambda);
% lambda_new_log      =   log(lambda_new);
% m_im_log_int        =   interp1(lambda_log,m_im_log,lambda_new_log);
% 
% % back transform
% m_im_int            =   exp(m_im_log_int);
% 
% % repeat for mi_re
% m_re_int            =   interp1(lambda_log,m_re,lambda_new_log);
% 
% % transpose and covert to table
% m_imag              =   m_im_int';
% m_real              =   m_re_int';
% lambda              =   lambda_new';

% The first time I ran this, i used 10e6 instead of 1e6 to convert from um
% to m. So I processed from 25 nm to 303 nm instead of 250 nm to 3003. This
% means the library now has values from 25 nm to 3003 nm. Maybe the shorter
% wavelengths will be useful to someone. 

% correction. the data are useless because the mim and mre are indexed to
% the correct wavelength but the wavelength is wrong so the wavenumber is
% wrong too
