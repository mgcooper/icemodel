clean

%% 
savedata = false;
plotfigs = false;
thisyear = '2015';

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
pathdata = setpath('runoff/data/icemodel/spectral/','project');
pathsave = setpath('runoff/data/icemodel/spectral/','project');

%% read in the coart data
fname = fullfile(pathdata,'Qi_coart_',thisyear,'_clearsky.txt');
coart = readtable(fname);

% read in the liston data
miedat   = importdata([pathdata 'mie.dat']);
solardat = importdata([pathdata 'solar_in.dat']);

%% interpolate the coart data to the liston spectral grid
wavl     = coart.WVLS_um_;
Qi       = coart.TotalDown;
wavlnew  = solardat(:,1);
Qinew    = interp1(wavl,Qi,wavlnew);
Qicoart  = table(wavlnew,Qinew,'VariableNames',{'Lambda','Qi'});

%% create a readme
Qicoart.Properties.Description     =   ...
['These data are the clear sky spectral incoming total irradiance ' newline ...
'(direct + diffuse) at the surface for a sub Arctic summer atmosphere ' newline...
'with a standard background stratospheric aerosol and OPAC Arctic mixed ' newline...
'layer aerosol at 67.5N, -48.9E (center of Rio Behar) at 1600 GMT (~solar ' newline...
'noon local time) calculated using the COART model at ' newline...
'https://satcorps.larc.nasa.gov/jin/coart.html'];
%% save the data

if savedata == true
    save([path.save 'Qi_coart_' thisyear],'Qicoart');
    dlmwrite([path.save 'solar.dat'],[wavlnew Qinew], ...
       'delimiter',' ','precision','%12.6f');
end

% %% this stuff confirms that I need to write out solar.dat in units of W/m2/um
% 
% % now compare with the data in solar.dat to confirm that units are W/m2/um
% 
% % at this point, Qi is W/m2/um, lambda is um, so dlambda is um
% dlambda             =   diff(wavl);
% dlambda(end+1)      =   dlambda(end);
% 
% % integrate my incoming solar irradiance
% total_solar         =   trapz(wavl,Qi); % = sum(Qi.*dlambda)
% Qi_scaled           =   Qi./total_solar;
% 
% % wavelength and solar radiation from solar.dat
% wavel_tmp           =   solardat(:,1);
% solar_tmp           =   solardat(:,2);
% 
% % integrate Liston's solar irradiance
% total_solar_tmp     =   trapz(wavel_tmp,solar_tmp);
% solar_scaled        =   solar_tmp./total_solar_tmp;
% 
% % now reproduce Liston's interpolation/integration to confirm method
% 
% % wavelength from mie.dat
% wavelength          =   (miedat(142,:))'; % just take one
% nvalues             =   length(wavelength);
% isolarvals          =   250;
% 
% % interpolate solar.dat to mie.dat wavelengths
% for k=1:nvalues
%     x = wavelength(k,1);
%     for i=1:isolarvals-1    
%         if (x > wavel_tmp(i))
%             icount = i;
%         end
%     end
%     x1 = wavel_tmp(icount);
%     x2 = wavel_tmp(icount+1);
%     y1 = solar_tmp(icount);
%     y2 = solar_tmp(icount+1);
%     solar(k) = y1 + (x - x1) * (y2 - y1)/(x2 - x1);
% end
% solar = solar';
% 
% % get dlambda. 
% dwavelen(1)         =   2.0*(wavelength(2,1)-wavelength(1,1));
% for k=2:nvalues-1
%     dwavelen(k)     =   (wavelength(k+1,1)-wavelength(k-1,1))/2.0;
% end
% dwavelen(nvalues)   =   2.0*(wavelength(nvalues,1)-wavelength(nvalues-1,1));
% dwavelen            =   dwavelen';
% 
% % compare with diffs
% dtest               =   diff(wavelength);
% dtest(1)            =   2.0*(wavelength(2,1)-wavelength(1,1));
% dtest(nvalues)      =   2.0*(wavelength(nvalues,1)-wavelength(nvalues-1,1));
% figure;
% myscatter(dtest,dwavelen);
% addOnetoOne
% 
% %
% Qsi = 0.0;
% for k=1:nvalues
%     Qsi = Qsi + solar(k) * dwavelen(k);
% end
% total_solar_check = Qsi
% total_solar_check2 = trapz(wavelength,solar)
% total_solar_check3 = sum(solar.*dwavelen)
% 
% %% plot things for comparison
% if plotfigs == 1
% figure
% plot(coart.WVLS_um_,coart.TotalDown); hold on;
% plot(coart.WVLS_um_,coart.Dir_Down); 
% plot(coart.WVLS_um_,coart.Dif_Down); 
% xlabel('\lambda [um]');
% ylabel('Qi [W m^{-2} um^{-1}]')
% legend('Total Downward','Direct Downward','Diffuse Downward'); 
% 
% % this is what I will pass into icemodel
% figure; 
% plot(wavl,Qi); hold on;
% plot(wavel_tmp,solar_tmp)
% xlabel('\lambda [um]');
% ylabel('Qi [W m^{-2} um^{-1}]');
% legend(['My data \Sigma Qi d\lambda = ' printf(sum(Qi.*dlambda),0) ' W m^{-2}'], ...
%     ['Liston \Sigma Qi d\lambda = ' printf(trapz(wavel_tmp,solar_tmp),0) ' W m^{-2}'])
% 
% % this replicates Liston et al. 1999 Figure 4
% figure; 
% plot(wavl,Qi_scaled); hold on;
% plot(wavel_tmp,solar_scaled)
% xlabel('\lambda [um]');
% ylabel('Q\lambda / Q_{total} [um^{-1}]');
% legend(['\Sigma Qi d\lambda = ' printf(sum(Qi./sum(Qi)),0) ' W m^{-2}'])
% 
% % plot for comparison
% figure;
% plot(wavel_tmp,solar_tmp); hold on;
% plot(wavelength,solar,'--');
% legend('solar.dat','interpolated to mie.dat');
% 
% end
% 
