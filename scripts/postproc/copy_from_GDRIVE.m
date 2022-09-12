clean
recycle off

save_data   =   true;
y1          =   [2009,2009,2009,2009];
y2          =   [2018,2018,2018,2018];
region      =   {'region','region','region','region'};
model       =   {'skinmodel','skinmodel','skinmodel','skinmodel'};
climate     =   {'t1','t2','t3','t4'};
albedo      =   {'modis','modis','modis','modis'};
folds       =   {'enbal','ice1','ice2','opts'};
%==========================================================================
%% loop through each run and copy/paste
%==========================================================================

for ii = 1:length(model)
    y1_ii   =   y1(ii);
    y2_ii   =   y2(ii);
    reg_ii  =   region{ii};
    mod_ii  =   model{ii};
    clim_ii =   climate{ii};
    alb_ii  =   albedo{ii};
    
    % set paths
%==========================================================================    
    p.data  = ['/Volumes/GDRIVE/matlab/GREENLAND/runoff/icemodel/model/' ...
                'experiment/v8/' reg_ii '/' mod_ii '/output/' alb_ii    ...
                '/' clim_ii '/'];
    p.save  = ['/Volumes/Samsung_T5/matlab/GREENLAND/runoff/icemodel/model/' ...
                'experiment/v8/' reg_ii '/' mod_ii '/output/' alb_ii    ...
                '/' clim_ii '/'];

    % copy the data
%==========================================================================
    nyrs    =   y2_ii-y1_ii+1;
    for i = 1:nyrs
        yr      = y1_ii+i-1;
        for j = 1:length(folds)
            pcopy   =   [p.data int2str(yr) '/' folds{j} '/'];
            ppaste  =   [p.save int2str(yr) '/' folds{j} '/'];
            list    =   getlist(pcopy,'*.mat');
            for k = 1:length(list)
                fcopy   =   [pcopy list(k).name]; 
                fpaste  =   [ppaste list(k).name];

                copyfile(fcopy,fpaste);
                delete(fcopy)
            end
        end
    end
end

