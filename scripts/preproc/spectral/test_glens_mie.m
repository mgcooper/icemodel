clean

%%

homepath            =   pwd;

if strcmp(homepath(2:6),'Users')
    path.data       =   ['/Users/mattcooper/Dropbox/CODE/MATLAB/'       ...
                        'GREENLAND/runoff/icemodel/data/preprocess/input/'];
    path.save       =   ['/Users/mattcooper/Dropbox/CODE/MATLAB/'       ...
                            'GREENLAND/icemodel/figs/'];
elseif strcmp(homepath(10:16),'mcooper')    
end

%% read in the 'mie.dat' data
mie_data            =   importdata([path.data 'mie.dat']);

% each coefficient is a 47x118 vector, grain size x wavelength, where grain
% size ranges from 0.005 - 10.00 mm, lambda varies from 0.299 um to 2.999
% um with 

g                   =   mie_data(1:47,:);
Qext                =   mie_data(48:94,:);
w                   =   mie_data(95:141,:);
lambda              =   mie_data(142,:); % just take one

% icemelt.f line 1507 gives r_eff array
r_eff               =   [   0.005, 0.007, 0.010, 0.015, 0.020, ...
                            0.030, 0.040, 0.050, 0.065, 0.080, ...
                            0.100, 0.120, 0.140, 0.170, 0.200, ...
                            0.240, 0.290, 0.350, 0.420, 0.500, ...
                            0.570, 0.660, 0.760, 0.870, 1.000, ...
                            1.100, 1.250, 1.400, 1.600, 1.800, ...
                            2.000, 2.500, 3.000, 3.500, 4.000, ...
                            4.500, 5.000, 5.500, 6.000, 6.500, ...
                            7.000, 7.500, 8.000, 8.500, 9.000, ...
                            9.500,10.000];

diff_lambda         =   diff(lambda);                        
%% make sample plots

% plot for r_eff = 0.5 - 10.0 as in Liston et al. 1999 (don't bother
% interpolating to 0.5 mm spacing)

si                  =   find(r_eff == 0.5);

figure(1);
p(1,:)              =   myplot(lambda,Qext(si:end,:)); ax = gca; 
ax.XLabel.String    =   'Wavelength (um)';
ax.YLabel.String    =   'Qext (-)';
offsetAxes

figure(2);
p(2,:)              =   myplot(lambda,1-w(si:end,:)); ax = gca; 
ax.XLabel.String    =   'Wavelength (um)';
ax.YLabel.String    =   '1-\omega (-)';
offsetAxes

figure(3);
p(3,:)              =   myplot(lambda,g(si:end,:)); ax = gca; 
ax.XLabel.String    =   'Wavelength (um)';
ax.YLabel.String    =   'g (-)';
offsetAxes

%% read in 'test_met_1hrly_dat'

% format is YYYY MM DD HH  

%% read in 'solar.dat'

% Glen's comments say "c Read Jerry Harrington's solar spectrum data."
solar_data          =   importdata([path.data 'solar.dat']);
min_lambda          =   solar_data(1,1);
max_lambda          =   solar_data(end,1);
% make a figure
figure(4)
p(4)                =   plot(solar_data(:,1),solar_data(:,2));

% the data is the solar spectrum from 289.99 to 380.99 nm in units um

%% this replicates GETSOLAR to confirm the units of solar.dat

solar_data          =   importdata([path.data 'solar.dat']);
% wavelength and solar radiation from solar.dat
wavel_tmp           =   solar_data(:,1);
solar_tmp           =   solar_data(:,2);

% wavelength from mie.dat
wavelength          =   lambda';
nvalues             =   length(wavelength);
isolarvals          =   250;

% interpolate solar.dat to mie.dat wavelengths
for k=1:nvalues
    x = wavelength(k,1);
    for i=1:isolarvals-1    
        if (x > wavel_tmp(i))
            icount = i;
        end
    end
    x1 = wavel_tmp(icount);
    x2 = wavel_tmp(icount+1);
    y1 = solar_tmp(icount);
    y2 = solar_tmp(icount+1);
    solar(k) = y1 + (x - x1) * (y2 - y1)/(x2 - x1);
end
solar = solar';

% plot for comparison
figure;
plot(wavel_tmp,solar_tmp); hold on;
plot(wavelength,solar,'--');
legend('solar.dat','interpolated to mie.dat');

% get dlambda. 
dwavelen(1)         =   2.0*(wavelength(2,1)-wavelength(1,1));
for k=2:nvalues-1
    dwavelen(k)     =   (wavelength(k+1,1)-wavelength(k-1,1))/2.0;
end
dwavelen(nvalues)   =   2.0*(wavelength(nvalues,1)-wavelength(nvalues-1,1));
dwavelen            =   dwavelen';

% compare with diffs
dtest               =   diff(wavelength);
dtest(1)            =   2.0*(wavelength(2,1)-wavelength(1,1));
dtest(nvalues)      =   2.0*(wavelength(nvalues,1)-wavelength(nvalues-1,1));
figure;
myscatter(dtest,dwavelen);
addOnetoOne

%
Qsi = 0.0;
for k=1:nvalues
    Qsi = Qsi + solar(k) * dwavelen(k);
end
total_solar = Qsi;

% compare with my method
sum(solar_tmp./sum(wavel_tmp))./sum(wavel_tmp)

trapz(wavelength,solar);

trapz(solar.*dtest)

sum(solar./wavelength);
sum(solar)/sum(wavelength);

sum(solar.*dtest)

figure
plot(solar)
            
      
solar(k) = y1 + (x - x1) * (y2 - y1)/(x2 - x1)


figure
plot(test(:,1),test(:,2));

figure
plot(test(:,1),test(:,2)./test(:,1)./sum(test(:,2)));


sum(test(:,1).*test(:,2))

