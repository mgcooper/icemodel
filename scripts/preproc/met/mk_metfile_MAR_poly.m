clean

% this script uses the new 'data' timetables to build catchment-scale
% metfiles. the only difference is that the metfile renames the variables
% for ease of use in the model, and there are options here to interpolate
% to higher resolution and save separate files with mar vs modis albedo 

% mk_metfile_MAR.m script calls makeMarMetfile.m function to build the
% point-scale metfiles. 


savedata =  true;
sitename =  'behar';
startyr  =  2015;
endyear  =  2016;
nyears   =  endyear-startyr+1;
interpT  =  true;
dtnew    =  '15m';

%--------------------------------------------------------------------------
% set paths
%--------------------------------------------------------------------------
pathdata =  ['/Users/coop558/mydata/mar3.11/matfiles/' sitename '/data/'];
pathsave =  setpath('GREENLAND/icemodel/input/met/');

%--------------------------------------------------------------------------
% read in the catchment to get the x,y coordinate 
%--------------------------------------------------------------------------

for n = 1:nyears
    
    thisyr  = num2str(startyr+n-1);
    fdata   = [pathdata 'MAR_' sitename '_' thisyr '.mat'];
    
%--------------------------------------------------------------------------
%     met     = marData2Met(fdata);
    metcopy = met;
%--------------------------------------------------------------------------

    if savedata == true
        
        fsave = ['met_' sitename '_MAR_' thisyr];

        save([pathsave fsave '_1hr.mat'],'met');

        % interpolate, re-check out-of-bound values, and save
        if interpT == true
            met      = interpMet(met,dtnew);
            met.date = datenum(met.Time);
            met      = metchecks(met,false); % false = don't plot

            save([pathsave fsave '_' dtnew '.mat'],'met');
        end
        
    end
end