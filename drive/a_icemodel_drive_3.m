clean

%%   run dependent information
savedata    =  true;
sitename    =  'region';
forcingdata =  'mar';
userdata    =  'modis';
uservars    =  'albedo';
meltmodel   =  'icemodel';

startyear   =   2008;
endyear     =   2018;
npts        =   1487;
si          =   1;              % 298, 671, 1046,  1420 
ei          =   743;             % 372, 744, 1116, 1487

% to complete the modis runs, need:
% 1:297
% 373:670
% 745:1045
% 1117:1419

% and re-save the ones I did after:
% 1. removing runoff (since it's just M-F) and surf_runoff
% 2. rename depth_melt and depth_freeze to melt and freeze
% 3. 

% 223+(372-223)/2
% 597+(744-597)/2
% 975+(1116-975)/2
% 1352+(1487-1352)/2

% MAKE SURE TO SET THE RIGHT PATH
drive          =  '/Users/mattcooper/data/';
% opts.pathmet   =  [drive 'mar3.11/matfiles/region/level2/old/'];
opts.pathmet   =  [drive 'mar3.11/matfiles/region/level2/lists/'];
opts.pathsave  =  ['/Volumes/Samsung_T5b/icemodel/output/region/v10/' userdata '/'];
opts.pathuser  =  '/Users/mattcooper/projects/icemodel/input/';

opts.savedata  =  savedata;

% this is unique to this version, which uses the new met files that have 
% both mar and modis albedo
opts.userdata  =  userdata;

%% run the model

disp([meltmodel ' ' userdata ' ' int2str(startyear) ':'   ...
        int2str(endyear) ' ' int2str(si) ':' int2str(ei)])

% initialize the model options
opts = a_opts_region(opts,sitename,meltmodel,startyear,endyear);

% run the model for each point
for ipt = si:ei % 1045 = rio behar
    
   opts.ipt    = ipt;
    
   % run the model
   [ice1,ice2] = icemodel_region_lists(opts);

end

% save the opts
if ~exist([opts.pathsave 'opts/'],'dir'); mkdir([opts.pathsave 'opts/']); end
save([opts.pathsave 'opts/opts_' meltmodel '_' userdata '.mat'],'opts');
        
% use this to check rio behar       
        
% if iyr == 8 && check_pt && strcmp(alb_model,'mar')
%     % v8 same point, mar albedo, most important check
%     tmp1 = load(['/Volumes/Samsung_T5/matlab/GREENLAND/icemodel/' ...
%                 'model/experiment/v8/region/skinmodel/output/mar/' ...
%                 'reference/2016/ice1/ice1_1045.mat']);
%     figure; plot(ice1.runoff); hold on; plot(tmp1.ice1.runoff,':')
%     legend('v8c, pt 1045, mar albedo','v8, pt 1045, mar albedo')
% 
%     % v8c rio behar kan-m forcings w/mar albedo
%     tmp2 = load(setpath(['GREENLAND/icemodel/model/versions/v8c/'   ...
%                     'output/2016/skinmodel/ice1_mar.mat']));
%     tmp3 = load(setpath(['GREENLAND/icemodel/model/versions/v8c/'   ...
%                     'eval/2016/mar_2016.mat']));
%     figure; plot(ice1.runoff); hold on; plot(tmp2.ice1.runoff,':')
%         plot(cumsum(tmp3.mar.runoff)); 
%     legend('v8c, pt 1045, mar albedo','v8c, behar kanm-mar albedo','mar, behar')
%                          
% % same as above, but if running modis albedo, so use kanm
% elseif iyr == 8 && check_pt && strcmp(alb_model,'modis')
%     
%     tmp1 = load(['/Volumes/Samsung_T5/matlab/GREENLAND/icemodel/' ...
%             'model/experiment/v8/region/skinmodel/output/modis/' ...
%             'reference/2016/ice1/ice1_1045.mat']);
%     figure; plot(ice1.runoff); hold on; plot(tmp1.ice1.runoff,':')
%     legend('v8c, pt 1045, modis','v8, pt 1045, modis')
%     
%     tmp2 = load(setpath(['GREENLAND/icemodel/model/versions/v8c/'   ...
%                     'output/2016/skinmodel/ice1_kanm.mat']));
%     tmp3 = load(setpath(['GREENLAND/icemodel/model/versions/v8c/'   ...
%                     'eval/2016/racmo_2016.mat']));
% 
%     figure; plot(ice1.runoff); hold on; plot(tmp2.ice1.runoff,':')
%         plot(cumsum(tmp3.racmo.runoff)); 
%     legend('v8c, pt 1045, modis alb','v8c, behar kanm alb','racmo, behar')
% 
%     end