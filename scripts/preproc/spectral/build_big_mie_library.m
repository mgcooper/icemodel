clean

use_warren84        =   false;
save_data           =   true;
post_process        =   false;
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
    m               =   readtable([path.data 'm_warren_84.xlsx']);
else
    m               =   readtable([path.data 'm_warren.xlsx']);
end

%% extract the values for the solar spectrum

si                  =   find(m.lambda == 0.25);
ei                  =   find(m.lambda == 3.003);
m_im                =   m.m_imag(si:ei);
m_re                =   m.m_real(si:ei);
m_lambda            =   m.lambda(si:ei); clear m;

% convert wavelength from [um] to [m]
m_lambda            =   m_lambda ./ 1e6;

% NOTE: I should have rounded the wavelengths before running mie, some of
% them had weird precision out to many decimal points

% these are the radii I am using, which range from very fine 
% grained snow (0.04 mm = 40 um) to very large grained blue ice (6 mm)
radii               =   [ 0.040, 0.050, 0.065, 0.080, 0.100, ...
                          0.120, 0.140, 0.170, 0.200, 0.240, ...
                          0.290, 0.350, 0.420, 0.500, 0.570, ...
                          0.660, 0.760, 0.870, 1.000, 1.100, ...
                          1.250, 1.400, 1.600, 1.800, 2.000, ...
                          2.250, 2.500, 2.750, 3.000, 3.500, ...
                          4.000, 4.500, 5.000, 5.500, 6.000 ];

nradii              =   length(radii);

% convert to meters
radii               =   radii ./ 1000;

%% Build the library in units of solid ice thickness
rho_i               =   917; 
rho_ice             =   917;

%% Calculate mie parameters. It takes ~five days to average over 100 radii 
n_avg                   =   1000;
r_sigma                 =   0.15; 

for j = 1:nradii

  % Grain radius [m]
    r_j                 =   radii(j);
    
  % 100 normally distributed values around r_i
    r_jj                =   normrnd(r_j,r_sigma*r_j,n_avg,1);
    r_dist(j,:)         =   r_jj;
    
    for k = 1:n_avg
        
      % the grain radius
        r_i             =   r_jj(k);
        
      % Number of ice grains per unit volume [m-3]
        N               =   (3/4)*(rho_i/rho_ice)*(1/(pi*r_i^3));

      % Wavenumber [m-1]
        k_num           =   2*pi./m_lambda;

      % Single scattering properties g, w, and Qext for each wavelength
        for n = 1:length(m_lambda)

          % wavenumber and complex index of refraction
            k_n         =   k_num(n);
            mi          =   1i*m_im(n);
            mr          =   m_re(n);
            m_n         =   mr + mi;
          % mie calculations
            mie_params  =   Mie(m_n,k_n*r_i);
            Qext(j,n,k) =   mie_params.Qext;            % extinction efficiency
            w(j,n,k)    =   mie_params.Omega;           % single-scattering albedo
            g(j,n,k)    =   mie_params.Asy;             % assymetry factor
            Rext(j,n,k) =   N*pi*(r_i^2)*Qext(j,n,k);   % radiance extinction coefficient

        end
    end
end

%%
if save_data == 1
    save([path.save 'mie_library_big_all_vars']); % save it all
end

%% post process the data

if post_process == 1

Qext_med            =   nanmedian(Qext,3);
Rext_med            =   nanmedian(Rext,3);
w_med               =   nanmedian(w,3);
g_med               =   nanmedian(g,3);

Qext_avg            =   nanmean(Qext,3);
Rext_avg            =   nanmean(Rext,3);
w_avg               =   nanmean(w,3);
g_avg               =   nanmean(g,3);

lambda              =   m_lambda; clear m_lambda;
lambda              =   roundn(10^9.*lambda,0);             % convert to nm


%% first plot each varible

figure                      % g                        
for n = 1:length(radii)
    r_n             =   radii(n);
    plot(lambda,g_avg(n,:)); hold on;
    plot(lambda,g_med(n,:)); 
    legend('avg','med');
    ylabel('g');
    pause
    hold off
end

figure                      % w
for n = 1:length(radii)
    r_n             =   radii(n);
    plot(lambda,1-w_avg(n,:)); hold on;
    plot(lambda,1-w_med(n,:)); set(gca,'YScale','log');
    legend('avg','med');
    ylabel('1-\omega');
    pause
    hold off
end


figure                      % Qext
for n = 1:length(radii)
    r_n             =   radii(n);
    plot(lambda,Qext_avg(n,:)); hold on;
    plot(lambda,Qext_med(n,:)); %set(gca,'YScale','log');
    legend('avg','med');
    ylabel('Qext');
    pause
    hold off
end

figure                      % Rext
for n = 1:length(radii)
    r_n             =   radii(n);
    plot(lambda,Rext_avg(n,:)); hold on;
    plot(lambda,Rext_med(n,:)); %set(gca,'YScale','log');
    legend('avg','med');
    ylabel('Rext');
    pause
    hold off
end


%% overlap with the prior compilation
figure;     % omega
for n = 1:length(radii_overlap)
    r_n             =   radii_overlap(n);
    i84             =   find(radii_84 == r_n);
    i08             =   find(radii_08 == r_n);
    
    plot(lambda08,1-w08(i08,:)); hold on;
    plot(lambda84,w84(i84,:)); 
    legend('08','84');
    ylabel('1-\omega');
    set(gca,'YScale','log');
    pause
    hold off
end

figure;     % Qext
for n = 1:length(radii_overlap)
    r_n             =   radii_overlap(n);
    i84             =   find(radii_84 == r_n);
    i08             =   find(radii_08 == r_n);    

    plot(lambda08,Qext08(i08,:)); hold on;
    plot(lambda84,Qext84(i84,:)); 
    legend('08','84');
    ylabel('Qext');
    pause
    hold off
end


%% Smooth Qext

% Qext - smooth up to 2725, then concatenate
nlambda             =   length(lambda);
ismooth0            =   find(lambda == 350);
ismooth1            =   find(lambda == 550);
ismooth2            =   find(lambda == 2725);
i00                 =   ismooth0;
i0                  =   ismooth1;
i1                  =   ismooth2;
i2                  =   i1+1;
i3                  =   nlambda;

order               =   3;
window              =   31;

figure;     
for n = 1:length(radii_overlap)
    
    r_n             =   radii_overlap(n);
    i84             =   find(radii_84 == r_n);
    i08             =   find(radii_08 == r_n);   

% smooth the data with a Savitsky Golay filter    
    Qext08smooth1   =   sgolayfilt(Qext(i08,1:i1),order,window,[],2);
    
% deal with the ends
    x               =   Qext08smooth1;
    lx              =   length(x);
    m               =   (window-1)/2;
    B               =   sgolay(order,window);
    steady          =   conv(x,B(m+1,:),'same');
    ybeg            =   B(1:m,:)*(x(1:window))';
    yend            =   B(window-m+1:window,:)*(x(lx-window+1:lx))';
    ybeg            =   ybeg';
    yend            =   yend';
    
    Qext08smooth1               =   steady;
    Qext08smooth1(1:m)          =   ybeg;
    Qext08smooth1(lx-m+1:lx)    =   yend;
    
% put the unaltered end back
    Qext08smooth1(i2:i3)    =   Qext(i08,i2:i3);    
        
% plot the results    
    plot(lambda84,Qext84(i84,:)); hold on; 
    plot(lambda,Qext(i08,:));
    plot(lambda,Qext08smooth1,'-','LineWidth',4)
    legend('84','08','08 - Smooth');
    ylabel('Qext');
    title([int2str(r_n) ' mm']);
    pause
    hold off
    clear Qext08smooth Qext08smooth1 ybeg yend x 
end


end





