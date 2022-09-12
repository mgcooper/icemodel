clean
%==========================================================================

% DON'T RENAME!
A
% Oct 2020: I am nearly certain this is the script that was used to build
% the final library. This is because 1) it was the most recently modified
% when i checked, 2) the final library is 'v3' and should be '08' as well

% Dec 2020: I copied this to my 2018 field folder to re-run for grain radii
% up to 12 mm, and to simplify the naming I renamed it build_mie_library,
% but I am leaviong this one here for reference with the same name

% I started to re-write this hence the line brreaks then switched over to
% the field experiment folder version
%==========================================================================
save_data   =   0;
plot_figs   =   0;

%==========================================================================
%% set paths
%==========================================================================
path.data   =   ['/Volumes/Samsung_T5/matlab/GREENLAND/runoff/icemodel/' ...
                    'data/spectral/mie_library_v3/'];
path.save   =   ['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/' ...
                'runoff/icemodel/data/preprocess/spectral/mie_library_v3/'];
%==========================================================================
%% load 
%==========================================================================

load([path.data 'mie_library_big_all_vars'],'Qext','g','w','m_lambda','radii');

% take the average (appears to work better than the median)
Qext08              =   nanmean(Qext,3);
w08                 =   nanmean(w,3);
g08                 =   nanmean(g,3);

% convert to nm
lambda08            =   m_lambda; clear m_lambda;
lambda08            =   roundn(10^9.*lambda08,0);

% convert omega to 1-omega
w08                 =   1-w08;
%%

% load glen's original data for comparison
mie_data            =   importdata([path.data 'mie.dat']);
g84                 =   mie_data(1:47,:);
Qext84              =   mie_data(48:94,:);
w84                 =   mie_data(95:141,:);
lambda84            =   1000.*mie_data(142,:); % just take one

%% make a couple of plots to be sure things are identical

if plot_figs == 1
% the first index in 08 is 0.04 mm. This is index 7 in 84
    figure; 
    plot(lambda84,g84(7,:)); hold on;
    plot(lambda08,g08(1,:));
    legend('84','08');

    figure; 
    plot(lambda84,w84(7,:)); hold on;
    plot(lambda08,w08(1,:));
    legend('84','08');

    figure; 
    plot(lambda84,Qext84(7,:)); hold on;
    plot(lambda08,Qext08(1,:));
    legend('84','08');
end

%% interpolate the new data to the old spectral bads

method              =   'spline';
g08new              =   interp1(lambda08,g08',lambda84,method);
g08new              =   g08new';

% omega
w08new              =   interp1(lambda08,w08',lambda84,method);
w08new              =   w08new';

% Qext
Qext08new           =   interp1(lambda08,Qext08',lambda84,method);
Qext08new           =   Qext08new';

%% plot the results

if plot_figs == 1
figure
for n = 1:size(g08new,1)
    plot(lambda08,g08(n,:)); hold on;
    plot(lambda84,g08new(n,:));
    legend('original','interpolated');
    pause
    hold off
end

figure
for n = 1:size(w08new,1)
    plot(lambda08,w08(n,:)); hold on;
    plot(lambda84,w08new(n,:));
    legend('original','interpolated');
    pause
    hold off
end


figure
for n = 1:size(Qext08new,1)
    plot(lambda08,Qext08(n,:)); hold on;
    plot(lambda84,Qext08new(n,:));
    legend('original','interpolated');
    pause
    hold off
end
end

%%

% Qext - smooth up to 2725, then concatenate
nlambda             =   length(lambda08);
ismooth0            =   find(lambda08 == 350);
ismooth1            =   find(lambda08 == 550);
ismooth2            =   find(lambda08 == 2725);
i00                 =   ismooth0;
i0                  =   ismooth1;
i1                  =   ismooth2;
i2                  =   i1+1;
i3                  =   nlambda;

order               =   3;
window              =   13;

figure;     
for n = 1:length(radii)
    
    r_n             =   radii(n);
% smooth the data with a Savitsky Golay filter    
    Qext08smooth1   =   sgolayfilt(Qext08(n,1:i1),order,window,[],2);
    
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
    Qext08smooth1(i2:i3)    =   Qext08(n,i2:i3);
    
    Qext08smooth(n,:)       =   Qext08smooth1;
        
% plot the results    
%     plot(lambda08,Qext08(n,:)); hold on;
%     plot(lambda08,Qext08smooth1,'-','LineWidth',2)
%     legend('08','08 - Smooth');
%     ylabel('Qext');
%     title([int2str(r_n) ' mm']);
%     pause
%     hold off
    clear Qext08smooth1 ybeg yend x 
end                        

%% Repeat the interpolation to the 118 bands using smoothed Qext

Qext08new           =   interp1(lambda08,Qext08smooth',lambda84,method);
Qext08new           =   Qext08new';
                        
%% build a new array with my values

nclasses            =   size(g08,1);
lambda84            =   mie_data(142:end,:); % just take one
lambda84            =   lambda84(1:nclasses,:);
    
%
mie                 =   [   g08new;
                            Qext08new;
                            w08new;
                            lambda84  ];

%% plot for comparison to be sure

% these are the radii used in Glen's files with the 84 compilation
radii_84            =   [   0.005, 0.007, 0.010, 0.015, 0.020, ...
                            0.030, 0.040, 0.050, 0.065, 0.080, ...
                            0.100, 0.120, 0.140, 0.170, 0.200, ...
                            0.240, 0.290, 0.350, 0.420, 0.500, ...
                            0.570, 0.660, 0.760, 0.870, 1.000, ...
                            1.100, 1.250, 1.400, 1.600, 1.800, ...
                            2.000, 2.500, 3.000, 3.500, 4.000, ...
                            4.500, 5.000, 5.500, 6.000, 6.500, ...
                            7.000, 7.500, 8.000, 8.500, 9.000, ...
                            9.500,10.000];
                        
                        
% these are the radii I used with the 08 compilation
radii_08            =   [   0.040, 0.050, 0.065, 0.080, 0.100, ...
                            0.120, 0.140, 0.170, 0.200, 0.240, ...
                            0.290, 0.350, 0.420, 0.500, 0.570, ...
                            0.660, 0.760, 0.870, 1.000, 1.100, ...
                            1.250, 1.400, 1.600, 1.800, 2.000, ...
                            2.250, 2.500, 2.750, 3.000, 3.500, ...
                            4.000, 4.500, 5.000, 5.500, 6.000 ];

% the radii overlap at radii_84(7:31) = radii_08(1:25) and then 
% radii_84(32) == radii_08(27)
% radii_84(33) == radii_08(29)
% radii_84(34:39) == radii_08(30:end)

radii_overlap       =   [   0.040, 0.050, 0.065, 0.080, 0.100, ...
                            0.120, 0.140, 0.170, 0.200, 0.240, ...
                            0.290, 0.350, 0.420, 0.500, 0.570, ...
                            0.660, 0.760, 0.870, 1.000, 1.100, ...
                            1.250, 1.400, 1.600, 1.800, 2.000, ...
                            2.500, 3.000, 3.500, 4.000, 4.500, ...
                            5.000, 5.500, 6.000] ;

figure
for n = 1:length(radii_overlap)
    r_n             =   radii_overlap(n);
    i84             =   find(radii_84 == r_n);
    i08             =   find(radii_08 == r_n);   
    
%     plot(lambda84,g84(i84,:)); hold on;
%     plot(lambda84,g08new(i08,:),'--');
%     legend('84','08');
%     pause
%     hold off
    
%     plot(lambda84,w84(i84,:)); hold on;
%     plot(lambda84,w08new(i08,:),'--');
%     legend('84','08');
%     pause
%     hold off
    
    plot(lambda84,Qext84(i84,:)); hold on;
    plot(lambda84,Qext08new(i08,:),'--');
    legend('84','08');
    pause
    hold off
    
end
%% save the data

if save_data == 1
    save([path.save 'mie.mat'],'mie');
end

%% This shows the smoothing and compares with the '84 compilation
                        
% Qext - smooth up to 2725, then concatenate
nlambda             =   length(lambda08);
ismooth0            =   find(lambda08 == 350);
ismooth1            =   find(lambda08 == 550);
ismooth2            =   find(lambda08 == 2725);
i00                 =   ismooth0;
i0                  =   ismooth1;
i1                  =   ismooth2;
i2                  =   i1+1;
i3                  =   nlambda;

order               =   3;
window              =   13;

figure;     
for n = 1:length(radii_overlap)
    
    r_n             =   radii_overlap(n);
    i84             =   find(radii_84 == r_n);
    i08             =   find(radii_08 == r_n);   

% smooth the data with a Savitsky Golay filter    
    Qext08smooth1   =   sgolayfilt(Qext08(i08,1:i1),order,window,[],2);
    
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
    Qext08smooth1(i2:i3)    =   Qext08(i08,i2:i3);
    
    Qext08smooth(n,:)       =   Qext08smooth1;
        
% plot the results    
%     plot(lambda84,Qext84(i84,:)); hold on; 
%     plot(lambda08,Qext08(i08,:));
%     plot(lambda08,Qext08smooth1,'-','LineWidth',2)
%     legend('84','08','08 - Smooth');
%     ylabel('Qext');
%     title([int2str(r_n) ' mm']);
%     pause
%     hold off
%     clear Qext08smooth Qext08smooth1 ybeg yend x 
end
