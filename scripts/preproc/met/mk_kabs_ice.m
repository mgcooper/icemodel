clean

use_warren84        =   false;
save_data           =   true;
%%

homepath            =   pwd;

if strcmp(homepath(2:6),'Users')
    path.data       =   ['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/' ...
                            'runoff/icemodel/data/preprocess/spectral/'];
    path.save       =   ['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/' ...
                            'runoff/icemodel/data/preprocess/spectral/'];
elseif strcmp(homepath(10:16),'mcooper')    
end

%% read in the complex index of refraction of pure ice

if use_warren84 == true
    m       =   readtable([path.data 'm_warren_84.xlsx']);
else
    m       =   readtable([path.data 'm_warren.xlsx']);
end

%% read in the 'mie.dat' data from the original model

mie_data    =   importdata([path.data 'mie.dat']); 

% convert the wavelengths to meters
lambda      =   mie_data(142,:)./1e6; 
nvalues     =   length(lambda); clear mie_data

%% extract the values for the solar spectrum

si          =   find(m.lambda == 0.25);
ei          =   find(m.lambda == 3.003);
m_im        =   (m.m_imag(si:ei))';
m_re        =   (m.m_real(si:ei))';
m_lambda    =   (m.lambda(si:ei))'; clear m;

% convert wavelength from [um] to [m]
m_lambda    =   m_lambda ./ 1e6;

% calculate the absorption coefficient of pure ice on the m_im grid
Kabs_m      =   4*pi./m_lambda.*m_im; 

% interpolate in log space
Kabs1       =   interp1(log(m_lambda),log(Kabs_m),log(lambda),'pchip');
Kabs1       =   exp(Kabs1);

% interpolate m_im first then compute Kabs
m_im        =   interp1(log(m_lambda),log(m_im),log(lambda),'pchip');
m_im        =   exp(m_im);

% calculate the absorption coefficient of pure ice
Kabs2       =   4*pi./lambda.*m_im; 

% plot the original and interpolated values
figure;
plot(m_lambda,Kabs_m); hold on;
plot(lambda,Kabs1);
plot(lambda,Kabs2);
set(gca,'YScale','log');
legend('original','kabs interpolated','m im interpolated');

%% use the second one and save it

kabs_theory         =   Kabs2;

%% process my field profile

load(['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/field/2018/data/' ...
        'processed/20july/f_coefficients/Kabs']);
kabs_user           =   Kabs.interp.Kabs;
lambda_user         =   Kabs.interp.Lambda;

% filter the noisy field values from 350-544
si                  =   1;
ei                  =   find(lambda_user == 600);
kabs_2              =   kabs_user(si:ei);
kabs_2              =   fillmissing(kabs_2,'movmean',3);
kabs_2_filt         =   sgolayfilt(kabs_2,3,21);
lambda_2            =   lambda_user(si:ei)./1e9;

% extend the field data to 350 nm using linear regression
si                  =   find(lambda_user == 360);
ei                  =   find(lambda_user == 390);
x                   =   (360:390)';
y                   =   kabs_2_filt(si:ei);
lm                  =   fitlm(x,y,'linear');
m                   =   lm.Coefficients.Estimate(2); % slope
b                   =   lm.Coefficients.Estimate(1); % intercept
xpred               =   (298:360)';
ypred               =   m.*xpred + b;
lambda_1            =   xpred./1e9;
kabs_1              =   ypred;

% combine the extrapolated values (300:360) with the filtered values (361:599)
si                  =   find(lambda_2 == 3.61e-7);
ei                  =   find(lambda_2 == 5.99e-7);
lambda_12           =   [lambda_1; lambda_2(si:ei)];
kabs_12             =   [kabs_1; kabs_2_filt(si:ei)];

% filter it again
kabs_12_filt        =   sgolayfilt(kabs_12,3,13);

% combine the filtered field values with the theoretical estimate
si                  =   find(m_lambda == 600/1e9);
kabs_3              =   (Kabs_m(si:end))';
lambda_3            =   (m_lambda(si:end))';

% Merge them into Kext_FullSpectrum_Experimental (Kext_FSE)
% kabs_FS             =   [kabs_1;kabs_2;kabs_3];
% lambda_FS           =   [lambda_1;lambda_2;lambda_3].*1e6;
kabs_FS             =   [kabs_12;kabs_3];
lambda_FS           =   [lambda_12;lambda_3].*1e6;

% plot the results
figure;
plot(m_lambda.*1e6,Kabs_m); hold on;
plot(lambda_FS,kabs_FS);
set(gca,'YScale','log')


%% interpolate the values to the 118 spectral bands used by IceModel
lambda_FS           =   lambda_FS./1e6;
kabs_user           =   interp1(log(lambda_FS),log(kabs_FS),log(lambda));
kabs_user           =   exp(kabs_user);

% plot the results
figure;
plot(m_lambda,Kabs_m); hold on;
plot(lambda_FS,kabs_FS);
plot(lambda,kabs_user);
legend('Theory','Field','Interpolated');
set(gca,'YScale','log')


%% save the data
kice                =   kabs_theory;
kabs                =   kabs_user;

if save_data == 1
    save([path.save 'kice'],'kice');
    save([path.save 'kabs'],'kabs');
end
