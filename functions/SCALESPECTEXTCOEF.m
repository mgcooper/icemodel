%--------------------------------------------------------------------------
%   scale the theoretical spectral extinction coefficients by a
%   user-provided absorption coefficient profile from 0.3-0.9 um 
%--------------------------------------------------------------------------

function spect_extcoefs = SCALESPECTEXTCOEF(spect_extcoefs,wavelength, ...
                          kice,kabs)

% Use Eq. 15 of Warren et al. 2006 to scale the spectral extinction
% coefficients calculated using the asymptotic two-stream solution by a
% user-defined absorption spectrum for ice or snow.

% Stephen G. Warren, Richard E. Brandt, and Thomas C. Grenfell Visible and
% near-ultraviolet absorption spectrum of ice from transmission of solar
% radiation into snow, Applied Optics, Vol. 45 No. 21
	 
    wavel           =   wavelength(1,:);  % [um] (micrometers)
    
% pull out the reference values for kice (the value at 0.60 um)
% the spectral bands used here do not include 0.60 um, so interpolate
% between 0.597 and 0.616 um (indices 14 and 15). I could just as easily
% hard-code the value from Warren et al. 2008 (0.122 m-1)
    dwavel          =   wavel(15)-wavel(14);
    dk              =   kice(15)-kice(14);
    wavel_0         =   0.60;
    kice_0          =   kabs(14)+(dk/dwavel)*(wavel_0-wavel(14));
    
% get the reference value for kext
    dk              =   spect_extcoefs(15)-spect_extcoefs(14);
    kext_0          =   spect_extcoefs(14)+(dk/dwavel)*(wavel_0-wavel(14));

% apply Eq. 15:
    kext            =   kext_0 .* sqrt(kabs ./ kice_0);

% merge the scaled values from 0.3-0.9 um with the rest of the spectrum 
% spect_extcoefs_scaled  =   [kext(1:24) spect_extcoefs(25:end)];
    spect_extcoefs  =   horzcat(kext(1:24),spect_extcoefs(25:end));
    
% use this to see the result
%     figure;
%     plot(wavel,spect_extcoefs); hold on; set(gca,'YScale','log');
%     plot(wavel,spect_extcoefs_scaled);
%     legend('theory','user defined');
