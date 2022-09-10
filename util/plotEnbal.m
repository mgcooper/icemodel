function [h1,h2,h3] = plotEnbal(enbal,met,varargin)
   
p = MipInputParser;
p.FunctionName = 'plotEnbal.m';
p.addRequired('enbal',@(x)istimetable(x));
p.addRequired('met',@(x)istimetable(x));
p.parseMagically('caller');


% % plot icemodel enbal
% h1 = plotIceEnbal(enbal); axis tight;
% 
% % plot forcing data enbal
% h2 = plotMetEnbal(met); axis tight;

% plot the comparison 
h3 = plotCompEnbal(enbal,met);

end



function h = plotCompEnbal(enbal,met)
   
   % moving mean window size
   win         =  7;

   % retime to daily and monthly
   enbald      =  retime(enbal,'daily','mean');
   enbalm      =  retime(enbal,'monthly','mean');
   metd        =  retime(met,'daily','mean');
   metm        =  retime(met,'monthly','mean');
   
   % replicate figure 12 in van As
   vars        =  {'swd','swu','lwd','lwu','shf','lhf','netr'};
   col         =  {'dark blue','blue','red','orange','dark green','green','purple'};
   h.line.fig  =  macfig('size','full'); hold on;
   % subtight(4,1,1:3,'style','fitted'); 
   for n = 1:numel(vars)
      hh(n) = plot(metd.Time,movmean(metd.(vars{n}),win),'Color',rgb(col{n})); 
      plot(enbald.Time,movmean(enbald.(vars{n}),win),':','Color',rgb(col{n}));
   end
   title('icemodel versus forcing');
   set(gca,'XAxisLocation','origin'); axis tight;
   ltext = {'SW\downarrow','SW\uparrow','LW\downarrow','LW\uparrow','SHF','LHF','Rnet'};
   legend(hh,ltext,'Orientation','horizontal','Location','northeast')
   title('icemodel (dotted) versus forcing (solid)');

   
% % this isn't actually very helpful   
%    % in a smaller subplot plot the surface temp and albedo
%    subtight(4,1,4,'style','fitted');
%    yyaxis left
%    plot(metd.Time,movmean(metd.tsfc,win)); hold on;
%    plot(enbald.Time,movmean(enbald.tsfc,win),':');
%    ylabel('tsfc');
%    
%    yyaxis right
%    plot(metd.Time,movmean(metd.albedo,win)); hold on;
%    plot(enbald.Time,movmean(enbald.albedo,win),':');
%    ylabel('albedo');
   
   
%    mticks = datetime(simYears,1,1):calmonths(1):datetime(simYears+1,1,1);
%    mticks.TimeZone = met.Time.TimeZone;
% 
%    figure; 
%    plot(met.Time,met.albedo); ylabel('albedo'); xtickformat('MMM'); 
%    xticks(mticks)


   
%%%%%%%%%%%%%%%%%%%%%

   macfig('size','large');
   subplot(2,2,1);myscatter(enbal.tsfc,met.tsfc); 
   ltext = sprintf('\nNSE Tsfc: %.2f\n',nse(enbal.tsfc,met.tsfc)); 
   addOnetoOne; xylabel('icemodel','forcing'); title('Tsfc');
   legend(ltext,'Location','nw');
   
   subplot(2,2,2);myscatter(enbal.netr,met.netr);
   addOnetoOne; xylabel('icemodel','forcing'); title('Net Rad');
   ltext = sprintf('\nNSE Rnet: %.2f\n',nse(enbal.netr,met.netr));
   legend(ltext,'Location','nw');

   subplot(2,2,3);myscatter(enbal.lhf,met.lhf);
   addOnetoOne; xylabel('icemodel','forcing'); title('LHF');
   ltext = sprintf('\nNSE LHF: %.2f\n',nse(enbal.lhf,met.lhf));
   legend(ltext,'Location','nw');

   subplot(2,2,4);myscatter(enbal.shf,met.shf);
   addOnetoOne; xylabel('icemodel','forcing'); title('SHF');
   ltext = sprintf('\nNSE SHF: %.2f\n',nse(enbal.shf,met.shf));
   legend(ltext,'Location','nw');

% 
% % % these are daily or hourly plots which are very noisy   
% 
%    % timeseries of upper grid cell diffusion length scale
%    figure; plot(enbal.Time,enbal.zD)
%    ylabel('zD')
%    
%    % timeseries of energy balance errors
%    figure; 
%    subplot(3,3,3); plot(enbal.Time,enbal.balance); ylabel('balance');
%    subplot(1,2,2); plot(enbal.Time,enbal.dt); ylabel('dt');
% 
%    % timeseries of latent heat flux and sensible heat flux
%    macfig;
%    subplot(2,1,1); plot(met.Time,met.lhf); hold on;
%    plot(enbal.Time,enbal.lhf,':'); ylabel('LHF'); legend('met','icemodel') 
% 
%    subplot(2,1,2); plot(met.Time,met.shf); hold on; 
%    plot(enbal.Time,enbal.shf,':'); ylabel('SHF'); legend('met','icemodel')
% 
%    % timeseries of tair vs tsfc
   
   metd     = retime(met,'daily','mean');
   enbald   = retime(enbal,'daily','mean');

   macfig('size','horizontal');
   subplot(2,1,1);
   plot(met.Time,met.tsfc); hold on; 
   plot(enbal.Time,enbal.tsfc,':'); ylabel('tsfc'); legend('met','icemodel')
   
   subplot(2,1,2);
   plot(metd.Time,metd.tsfc); hold on; 
   plot(enbald.Time,enbald.tsfc,':'); ylabel('tsfc'); legend('met','icemodel')
   
end




