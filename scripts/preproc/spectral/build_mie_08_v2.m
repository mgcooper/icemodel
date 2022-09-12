clean

path.data   =   ['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/' ...
                'runoff/icemodel/data/preprocess/spectral/mie_library_v2/'];
path.save   =   ['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/' ...
                'runoff/icemodel/data/preprocess/spectral/mie_library_v2/'];
%%

load([path.data 'mie_library_filtered']);

% convert wavelngth back to um
lambda08            =   lambda08./1000;

% convert omega to 1-omega
w08                 =   1-w08;
%%

% load glen's original data for comparison
mie_data            =   importdata([path.data 'mie.dat']);


g84                 =   mie_data(1:47,:);
Qext84              =   mie_data(48:94,:);
w84                 =   mie_data(95:141,:);
lambda84            =   mie_data(142,:); % just take one

%% make a couple of plots to be sure things are identical

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


%% build a new array with my values

nclasses            =   size(g08,1);
lambda84            =   mie_data(142:end,:); % just take one
lambda84            =   lambda84(1:nclasses,:);
    
%
mie                 =   [   g08new;
                            Qext08new;
                            w08new;
                            lambda84  ];

%% save the data

if save_data == 1
    save([path.save 'mie.mat'],'mie');
end

