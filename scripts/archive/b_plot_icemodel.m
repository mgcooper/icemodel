clean

% notes on where to pick up - need to figure out the plot at end of this
% script - what was being plotted in the fortran version vs. here, seems
% weird that i plotted melt_frac and then on top the cumulative meltflux

% after that, i need to define Cp according to Clark et al., and create
% liq, air, ice fraction variables. FIRST CREATE A NEW FOLDER. Then I can
% explore other things like qsfactor, changing density, etc. Then I can
% move on to spatial once the point model is tits

% at some point, I should make the spectral model separate, and compute
% kbulk as a function of ice density, and then update kbulk at each
% timestep as the ice density evolves

%%
save_figs           =   0;
plot_surface        =   1;
plot_subsurface     =   1;

%%
homepath            =   pwd;

if strcmp(homepath(2:6),'Users')
    path.data       =   ['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/' ...
                            'runoff/icemodel/model/v9/output/'];
    path.save       =   ['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/' ...
                            'runoff/icemodel/model/v9/figs/'];
elseif strcmp(homepath(10:16),'mcooper')    
end

%% 

% load the test data
load([path.data 'icemodel']);

hourly              =   icemodel.hourly;
daily               =   icemodel.daily;

% pull out the var names and units
vars                =   daily.vars;
varnames            =   daily.varnames;
units               =   daily.units;

% 
xhourly             =   hourly.dates(:,4);
xdaily              =   daily.dates(:,4);

%% plot some stuff

% no need to plot the meteorological stuff
metinds             =   [1,3,4,5,6,13]; % tair, rh, wspd, qsi, qli, albedo
surfinds            =   [2,7,8,9,10,11,12,14,15,16,17,18]; % tsfc, qle, qh, qe, qc, qm, balance, qsip, waterdepth, waterflux, freezedepth, freezeflux

if plot_surface == 1
for n = 1:length(vars)
    
%     if ~any(n == surfinds(:))
%         continue
%     end
    
    var_n           =   vars{n};
    varname_n       =   varnames{n};
    units_n         =   units{n};
    
    yhourly         =   hourly.(var_n);
    ydaily          =   daily.(var_n);
    
    figure('Units','inches','Position',[1 1 1.5*12 1.5*4])
    
    subplot(1,2,1)
    plot(xhourly,yhourly); datetick('x','mmm');
    ylabel(['hourly ' varname_n ' ' units_n]);
    
    subplot(1,2,2)
    plot(xdaily,ydaily); datetick('x','mmm');
    ylabel(['daily ' varname_n ' ' units_n ' hr^{-1}']);
    
%     if save_figs == 1
%         export_fig([path.save varname_n '.png'],'-r400');
%     end
end
% close all
end
%% plot the subsurface

% build a grid
z                   =   0:0.03:15;
zi                  =   60; % 1.8 m
si                  =   100; % doy 200
ei                  =   300; % end of year
t                   =   icemodel.daily.dates(si:ei,4);
[X,Y]               =   meshgrid(t,z(1:zi));

% if plot_subsurface == 1
for n = 1:length(ice2.units)
    var_n               =   ice2.vars{n};
    varname_n           =   ice2.varnames{n};
    units_n             =   ice2.units{n};
    Z                   =   ice2.(var_n);
    Z                   =   Z(1:zi,si:ei);
    f(n)                =   figure;
    h(n)                =   pcolor(X,Y,Z); ax(n) = gca; datetick('x');
    h(n).EdgeColor      =   'none';
    ax(n).YDir          =   'reverse';
    ax(n).YLabel.String =   'Depth (m)';
    ax(n).YLim          =   [0 2];
    c                   =   colorbar;
    set(get(c,'title'),'String',units_n);
    title(varname_n);
    
%     if save_figs == 1
%         export_fig([path.save 'z_' varname_n '.png'],'-r400');
%     end
end
% close all
% end

% cd(path.save)
%% these came from another script I wrote , might need to be added

% [nlayers,ndays]     =   size(ice2.melt_frac);
% X                   =   1:1:ndays;
% Y                   =   0.03:0.03:nlayers*0.03;
% 
% % plot melt_frac
% Z                   =   ice2.melt_frac;
% f                   =   figure; 
% h1                  =   pcolor(X,Y,Z); ax = gca; hold on;
% h2                  =   plot(X,daily.meltdepth,'m');
% h1.EdgeColor        =   'none';
% ax.YDir             =   'reverse';
% ax.YLim             =   [0 2];
% c                   =   colorbar;
% ax.XLabel.String    =   'Day of Year';
% ax.YLabel.String    =   'Melt Fraction [-]';
% legend('Melt Fraction','Cumulative Water Flux','Location','best')
% if save_figs == 1
%     export_fig([path.save 'z_melt_fraction_cumulative_overlay.png'],'-r400');
% end
    
%% plot freeze_frac
% if str2num(run(4)) >= 6
% Z                   =   ice2.freeze_frac;
% f                   =   figure; 
% h1                  =   pcolor(X,Y,Z); ax = gca; hold on;
% h2                  =   plot(X,cumsum(daily.freezeflux),'m');
% h1.EdgeColor        =   'none';
% ax.YDir             =   'reverse';
% c                   =   colorbar;
% ax.XLabel.String    =   'Day of Year';
% ax.YLabel.String    =   'Freeze Fraction [-]';
% legend('Refreeze Fraction','Cumulative Refreeze Flux','Location','best')
% end

% % water depth is meaningless, I think - it must be adding every day 
% f                   =   figure;
% h                   =   plot(X,daily.waterdepth); ax = gca; 
% ax.XLabel.String    =   'Day of Year';
% ax.YLabel.String    =   'Water Depth [m]';
% 
% % plot water flux
% f                   =   figure;
% h                   =   plot(X,daily.waterflux); ax = gca; 
% ax.XLabel.String    =   'Day of Year';
% ax.YLabel.String    =   'Water Flux [m]';
% 
% % plot cumulative water flux
% f                   =   figure;
% h                   =   plot(X,cumsum(daily.waterflux)); ax = gca; 
% ax.XLabel.String    =   'Day of Year';
% ax.YLabel.String    =   'Cumulative Water Flux [m]';
% ax.YDir             =   'reverse';


