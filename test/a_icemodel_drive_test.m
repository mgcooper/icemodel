clean
%--------------------------------------------------------------------------
%   run dependent information
%--------------------------------------------------------------------------

% this is to run the POSTPROC2_test and _test2 to see how much performance
% and/or storage can be improved relative to POSTPROC

% note: mar forcing, modis albedo, 2008, point 1117 was used for the tests

savedata    =  true;
sitename    =  'region';
forcingdata =  'mar';
userdata    =  'modis';
uservars    =  'albedo';
meltmodel   =  'icemodel';

testtype    =  'proc'; % 'save';


startyear   =   2008;
endyear     =   2008;
npts        =   1487;
si          =   1117;              % 1,   373, 745,  1117 
ei          =   1117;             % 372, 744, 1116, 1487

% MAKE SURE TO SET THE RIGHT PATH
drive          =  '/Users/coop558/mydata/';
opts.pathmet   =  [drive 'mar3.11/matfiles/region/level2/'];
opts.pathsave  =  [pwd '/test/'];
opts.pathuser  =  '/Users/coop558/MATLAB/GREENLAND/icemodel/input/';

opts.savedata  =  savedata;

%--------------------------------------------------------------------------
%% run the model
%--------------------------------------------------------------------------
disp([meltmodel ' ' userdata ' ' int2str(startyear) ':'   ...
        int2str(endyear) ' ' int2str(si) ':' int2str(ei)])

% initialize the model options
opts = a_opts_region(opts,sitename,meltmodel,startyear,endyear);
      
for ipt = si:ei % 1045 = rio behar

   opts.ipt    = ipt;

   % run the model
   switch testtype
      case 'save'
         [ice1,ice2,times] = icemodel_region_test_save(opts);
      case 'proc'
         [ice1,ice2,times] = icemodel_region_test_proc(opts);
   end
   
end

switch testtype
   case 'save'
      times.v4(1:2) = nan;
      save([opts.pathsave 'times_save'],'times');
   case 'proc'
      save([opts.pathsave 'times_proc'],'times');
end

% save the opts
% save([opts.pathsave 'opts/opts_' meltmodel '_' userdata '.mat'],'opts');
        
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