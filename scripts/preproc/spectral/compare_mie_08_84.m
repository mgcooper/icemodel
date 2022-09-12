clean

path.data           =   ['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/' ...
                            'runoff/icemodel/data/preprocess/spectral/'];
path.save           =   ['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/' ...
                            'runoff/icemodel/data/preprocess/spectral/'];                        
                        
%%

% % 
% m                   =   readtable([path.data 'm_warren.xlsx']);
% si                  =   find(m.lambda == 0.25);
% ei                  =   find(m.lambda == 3.003);
% m_im                =   m.m_imag(si:ei);
% m_re                =   m.m_real(si:ei);
% m_lambda            =   m.lambda(si:ei); 
% 
% % convert wavelength from [um] to [m]
% m_lambda            =   m_lambda ./ 1e6;


% read in Glen's files
mie_data            =   importdata([path.data 'mie.dat']);

% each coefficient is a 47x118 vector, grain size x wavelength, where grain
% size ranges from 0.005 - 10.00 mm, lambda varies from 0.299 um to 2.999
% um with 

g84                 =   mie_data(1:47,:);
Qext84              =   mie_data(48:94,:);
w84                 =   mie_data(95:141,:);
lambda84            =   mie_data(142,:); % just take one
lambda84            =   10^3.*lambda84; % convert from um to nm

% read in my results 
load([path.data 'mie_library_big_all_vars.mat']);
Qext08              =   nanmedian(Qext,3);
w08                 =   nanmedian(w,3);
g08                 =   nanmedian(g,3);
lambda08            =   m_lambda'; clear m_lambda;
lambda08            =   10^9.*lambda08;             % convert to nm
lambda08            =   roundn(lambda08,0);
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
                      
%%

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


%% plot the data


% smooth the data
order               =   3;
window              =   21;
dim                 =   2;

i0g                 =   find(lambda08 == 350);
i1g                 =   find(lambda08 == 1410);

figure      % g
for n = 1:length(radii_overlap)
    r_n             =   radii_overlap(n);
    i84             =   find(radii_84 == r_n);
    i08             =   find(radii_08 == r_n);

    g08smooth       =   sgolayfilt(g08(i08,i0g+1:i1g),order,window,[],dim);
    g08smooth       =   [g08(i08,1:i0g) g08smooth g08(i08,i1g+1:end)];
    
  % smooth just a tad around the kinks
%     g08smooth       =   smoothdata(g08smooth,'loess',13);
    
%     g08smooth       =   sgolayfilt(g08(i08,:),order,window,[],dim);
    
    plot(lambda84,g84(i84,:)); hold on;
    plot(lambda08,g08(i08,:));
    plot(lambda08,g08smooth,'LineWidth',4);
    
    legend('84','08','08 smoothed');
    ylabel('g');
    set(gca,'YScale','log');
    pause
    hold off
end

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

%% Smooth g 

% smooth the data
order               =   3;
window              =   21;
dim                 =   2;

i0g                 =   find(lambda08 == 350);
i1g                 =   find(lambda08 == 1410);


for n = 1:length(radii)

    g08smooth       =   sgolayfilt(g08(n,i0g+1:i1g),order,window,[],dim);
    g08smooth       =   [g08(n,1:i0g) g08smooth g08(n,i1g+1:end)];
    g08(n,:)        =   g08smooth;
    
end

% look at it
for n = 1:length(radii)
    plot(lambda08,g08(n,:));
    pause
end

%% (1-w) requires no smoothing but look at it

for n = 1:length(radii)
    plot(lambda08,w08(n,:));
    pause
end

%% Smooth Qext - smooth up to 2725, then concatenate

nlambda             =   length(lambda08);
ismooth0            =   find(lambda08 == 350);
ismooth1            =   find(lambda08 == 550);
ismooth2            =   find(lambda08 == 2725);
i00                 =   ismooth0;
i0                  =   ismooth1;
i1                  =   ismooth2;
i2                  =   i1+1;
i3                  =   nlambda;
  
for n = 1:length(radii)
    
% smooth the data with a Savitsky Golay filter    
    Qext08smooth    =   sgolayfilt(Qext08(n,1:i1),order,window,[],2);
    
% deal with the ends
    x               =   Qext08smooth;
    lx              =   length(x);
    m               =   (window-1)/2;
    B               =   sgolay(order,window);
    steady          =   conv(x,B(m+1,:),'same');
    ybeg            =   B(1:m,:)*(x(1:window))';
    yend            =   B(window-m+1:window,:)*(x(lx-window+1:lx))';
    ybeg            =   ybeg';
    yend            =   yend';
    
    Qext08smooth               =   steady;
    Qext08smooth(1:m)          =   ybeg;
    Qext08smooth(lx-m+1:lx)    =   yend;
    
% put the unaltered end back
    Qext08smooth(i2:i3)    =   Qext08(n,i2:i3);
    
% replace the original with the smoothed    
    Qext08(n,:)     =   Qext08smooth;

    clear Qext08smooth x lx m B steady ybeg yend 
end

% plot the results
figure;   
for n = 1:length(radii)
    plot(lambda08,Qext08(n,:),'-','LineWidth',2)
    legend('Qext');
    pause
end

%% save the output

if save_data == 1
    save([path.save 'mie_library_filtered'],'Qext08','w08','g08','lambda08','radii_08');
end

%% this was trial and error
% figure;     
% for n = 1:length(radii_overlap)
%     
%     r_n             =   radii_overlap(n);
%     i84             =   find(radii_84 == r_n);
%     i08             =   find(radii_08 == r_n);   
%     
% % smooth the data with a Savitsky Golay filter    
%     Qext08smooth1        =   smoothdata(Qext08(i08,1:i1),2,'sgolay',21);
%     Qext08smooth1(i2:i3) =   Qext08(i08,i2:i3);
% 
% % replace the first five with the original 
%     Qext08smooth1(1:i0)  =   Qext08(i08,1:i0);
%     
% % repeat the smoothing, this time with the first smoothed curve
%     Qext08smooth2        =   smoothdata(Qext08smooth1(1:i1),2,'sgolay',21);    
%     Qext08smooth2(i2:i3) =   Qext08(i08,i2:i3);
%     
% % replace the first three with the first smoothed 
%     Qext08smooth2(1:i0)  =   Qext08smooth1(1:i0);
%     
% % plot the results    
%     plot(lambda84,Qext84(i84,:)); hold on; 
%     plot(lambda08,Qext08(i08,:));
%     plot(lambda08,Qext08smooth1,'-','LineWidth',4)
%     plot(lambda08,Qext08smooth2,':','LineWidth',4)
%     legend('84','08','08 - Smooth 1','08 - Smooth 2');
%     ylabel('Qext');
%     title([int2str(r_n) ' mm']);
%     pause
%     hold off
%     clear Qext08smooth
%   
% end
% 
% figure
% plot(lambda08,g); hold on;
% plot(lambda84,g84(i84,:),':'); 
% legend('08','84');
% ylabel('g');
% 
% figure
% plot(lambda08,Qext); hold on;
% plot(lambda84,Qext84(i84,:),':'); 
% legend('08','84');
% ylabel('Qext');
% 
% figure
% plot(lambda08,g); hold on;
% plot(lambda84,g84(i84,:),':'); 
% legend('08','84');
% ylabel('g');
%     

