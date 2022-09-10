clean

comp_v2     =   true;
yyyy        =   '2016';
metfname    =   'met.mat';
model       =   'icemodel';
albedo      =   'kanm';

%--------------------------------------------------------------------------
%% set paths
%--------------------------------------------------------------------------
path.met    =   ['input/' yyyy '/met/'];
path.data   =   ['output/' yyyy '/' model '/'];
path.eval   =   ['eval/' yyyy '/'];

%--------------------------------------------------------------------------                
%% load model and field data
%--------------------------------------------------------------------------
load([path.met metfname]);
load([path.eval 'mar_' yyyy '.mat']);
load([path.eval 'merra_' yyyy '.mat']);
load([path.eval 'racmo_' yyyy '.mat']);
load([path.eval 'seb_' yyyy '.mat']);
load([path.data 'ice1_' albedo '.mat']);
load([path.data 'opts_' albedo '.mat']);
load([path.eval 'discharge_' yyyy '.mat']);
load([path.eval 'catchment_' yyyy '.mat']);

ice1keep = ice1;

% read in the comp
if comp_v2 == true
    ice1_v2 = load([path.data 'ice1_' albedo '_v2.mat']);
    opts_v2 = load([path.data 'opts_' albedo '_v2.mat']);
    ice1_v2 = ice1_v2.ice1;
    opts_v2 = opts_v2.opts;
end


% eeal with 
if size(ice1,1) == 2928 || size(ice1,1) == 4392
    ice1    = retime(ice1,datetime(mar.dates,'ConvertFrom','datenum'),  ...
                    'linear');
end

% FIGURE 1 = RUNOFF (options: 'raw','mean','members','sensitivity','surf')
if opts.skin_model == true
	h = plotrunoff(mar,merra,racmo,seb,ice1,discharge,catchment,'raw','surf');
else
    if comp_v2 == true
        plotrunoff(mar,merra,racmo,ice1_v2,ice1,discharge,catchment,'comp');
    else
        plotrunoff(mar,merra,racmo,seb,ice1,discharge,catchment,'raw');
    end
end
title(['albedo = ' albedo])

    

T   = ice1.Time(3649:end); 
R1  = ice1.runoff(3649:end); R1 = R1-R1(1);
R2  = cumsum(racmo.runoff(3649:end)); R2 = R2-R2(1);
R3  = cumsum(mar.runoff(3649:end)); R3 = R3-R3(1);
R4  = cumsum(merra.runoff(3649:end)); R4 = R4-R4(1);

figure; 
plot([T T T T],[R1 R2 R3 R4]); hold on;
legend('icemodel','racmo','mar','merra','location','best')
ylabel('Runoff'); title(['albedo = ' albedo])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([path.data 'enbal_' albedo '.mat']);
figure; plot(enbal.Time,enbal.Qg1); set(gca,'YLim',[-200 200])
for n = 1:height(enbal)
    if enbal.Tflag(n)==true
        vline(enbal.Time(n));
    end
end
ylabel('Qg'); title(['albedo = ' albedo])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([path.data 'ice2_' albedo '.mat']);

% pull out the field dates
si      =   discharge.plotUTC(1);
ei      =   discharge.plotUTC(end);
[s,e]   =   dateInds(si,ei,merra.dates);
T       =   datetime(merra.dates,'ConvertFrom','datenum');
t       =   datetime(merra.dates(s:e),'ConvertFrom','datenum');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check the upper layer gamma0 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure; plot(T,ice2.gamma(1,:)./ice2.gamma0(1,:));
% ylabel('top layer thermal k')
% legend('10*gamma/gamma');
% 
% % during the field study
% figure; plot(t,ice2.gamma(1,s:e)./ice2.gamma0(1,s:e));
% ylabel('top layer thermal k')
% legend('10*gamma/gamma');
% 
% % top 20 cm during the field study
% figure; 
% plot(t,mean(ice2.gamma(1:10,s:e))); hold on;
% plot(t,mean(ice2.gamma0(1:10,s:e)));
% ylabel('top 20 cm thermal k')
% legend('10*gamma','gamma');
% 
% figure; 
% plot(t,mean(ice2.gamma(1:10,s:e))./mean(ice2.gamma0(1:10,s:e)));
% ylabel('top 20 cm thermal k')
% legend('10*gamma/gamma');


% figure; plot(t,mean(ice2.gamma0(1:20,s:e))); ylabel('top 40 cm thermal k')


% figure; myscatter(enbal.Qc,-enbal.Qg1); addOnetoOne
% xlabel('Qc'); ylabel('Qg')
% 
% % figure; plot(enbal.Time,enbal.Qg1)
% % figure; plot(enbal.Time,enbal.Qm)
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check the upper layer thermal k to see if a smaller hmin smooths it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure; plot(T,ice2.gamma(1,:)); ylabel('top node thermal k')
% figure; plot(t,ice2.gamma(1,s:e)); ylabel('top node thermal k')
% figure; plot(t,mean(ice2.gamma(1:10,s:e))); ylabel('top 20 cm thermal k')
% figure; plot(t,mean(ice2.gamma(1:20,s:e))); ylabel('top 40 cm thermal k')
% 
% % plot the ice density
% figure; plot(T,ice2.ro_snow(1,:)); ylabel('top node density')
% figure; plot(t,ice2.ro_snow(1,s:e)); ylabel('top node density')
% figure; plot(t,mean(ice2.ro_snow(1:10,s:e))); ylabel('top 20 cm density')
% figure; plot(t,mean(ice2.ro_snow(1:20,s:e))); ylabel('top 40 cm density')
% figure; plot(t,mean(ice2.ro_snow(1:40,s:e))); ylabel('top 80 cm density')
% 
