clean

%% 
save_data           =   1;
plot_figs           =   1;
year                =   '2015';

%% coart settings (also see screenshots)
% 2015 07 21 (Julian Day 202)
% UTC: 1600
% 67.5
% -48.9
% 
% 2016 07 09 (Julian Day 190)
% UTC: 1600
% 67.5
% -48.9

%%
homepath            =   pwd;

if strcmp(homepath(2:6),'Users')
    path.delim      =   '/';
    path.data       =   ['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/' ...
                            'runoff/icemodel/data/preprocess/input/'];
    path.save       =   ['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/' ...
                            'runoff/icemodel/data/preprocess/output/' year '/'];
elseif strcmp(homepath(10:16),'mcooper')
end

%% read in the coart data
fname               =   [path.data 'Qi_coart_' year '_clearsky.txt'];
coart               =   readtable(fname);

% read in the liston data
mie_data            =   importdata([path.data 'mie.dat']);
solar_data          =   importdata([path.data 'solar_in.dat']);

%% interpolate the coart data to the liston spectral grid
lambda              =   coart.WVLS_um_;
Qi                  =   coart.TotalDown;
lambda_new          =   solar_data(:,1);
Qi_interp           =   interp1(lambda,Qi,lambda_new);
Qi_coart            =   table(lambda_new,Qi_interp, ...
                            'VariableNames',{'Lambda','Qi'});
%% create a readme
Qi_coart.Properties.Description     =   ...
['These data are the clear sky spectral incoming total irradiance ' newline ...
'(direct + diffuse) at the surface for a sub Arctic summer atmosphere ' newline...
'with a standard background stratospheric aerosol and OPAC Arctic mixed ' newline...
'layer aerosol at 67.5N, -48.9E (center of Rio Behar) at 1600 GMT (~solar ' newline...
'noon local time) calculated using the COART model at ' newline...
'https://satcorps.larc.nasa.gov/jin/coart.html'];
%% save the data

if save_data == 1
    save([path.save 'Qi_coart_' year],'Qi_coart');
    dlmwrite([path.save 'solar.dat'],[lambda_new Qi_interp], ...
        'delimiter',' ','precision','%12.6f');
end

%% this stuff confirms that I need to write out solar.dat in units of W/m2/um

% now compare with the data in solar.dat to confirm that units are W/m2/um

% at this point, Qi is W/m2/um, lambda is um, so dlambda is um
dlambda             =   diff(lambda);
dlambda(end+1)      =   dlambda(end);

% integrate my incoming solar irradiance
total_solar         =   trapz(lambda,Qi); % = sum(Qi.*dlambda)
Qi_scaled           =   Qi./total_solar;

% wavelength and solar radiation from solar.dat
wavel_tmp           =   solar_data(:,1);
solar_tmp           =   solar_data(:,2);

% integrate Liston's solar irradiance
total_solar_tmp     =   trapz(wavel_tmp,solar_tmp);
solar_scaled        =   solar_tmp./total_solar_tmp;

% now reproduce Liston's interpolation/integration to confirm method

% wavelength from mie.dat
wavelength          =   (mie_data(142,:))'; % just take one
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
total_solar_check = Qsi
total_solar_check2 = trapz(wavelength,solar)
total_solar_check3 = sum(solar.*dwavelen)

%% plot things for comparison
if plot_figs == 1
figure
plot(coart.WVLS_um_,coart.TotalDown); hold on;
plot(coart.WVLS_um_,coart.Dir_Down); 
plot(coart.WVLS_um_,coart.Dif_Down); 
xlabel('\lambda [um]');
ylabel('Qi [W m^{-2} um^{-1}]')
legend('Total Downward','Direct Downward','Diffuse Downward'); 

% this is what I will pass into icemodel
figure; 
plot(lambda,Qi); hold on;
plot(wavel_tmp,solar_tmp)
xlabel('\lambda [um]');
ylabel('Qi [W m^{-2} um^{-1}]');
legend(['My data \Sigma Qi d\lambda = ' printf(sum(Qi.*dlambda),0) ' W m^{-2}'], ...
    ['Liston \Sigma Qi d\lambda = ' printf(trapz(wavel_tmp,solar_tmp),0) ' W m^{-2}'])

% this replicates Liston et al. 1999 Figure 4
figure; 
plot(lambda,Qi_scaled); hold on;
plot(wavel_tmp,solar_scaled)
xlabel('\lambda [um]');
ylabel('Q\lambda / Q_{total} [um^{-1}]');
legend(['\Sigma Qi d\lambda = ' printf(sum(Qi./sum(Qi)),0) ' W m^{-2}'])

% plot for comparison
figure;
plot(wavel_tmp,solar_tmp); hold on;
plot(wavelength,solar,'--');
legend('solar.dat','interpolated to mie.dat');

end

