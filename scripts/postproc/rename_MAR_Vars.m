clean

% this is based on mk_metfile_MAR_poly. Instead that script, this one
% renames the vars in 'Data' and re-saves it. That should have been done in
% saveMarData. BE CAREFUL!

% update jul 2022, this is now iincorporated into saveMarData, so can prob
% be deleted eventually

saveData    =   true;
siteName    =   'L41A';
startYear   =   2012;
endYear     =   2012;
nYears      =   endYear-startYear+1;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% set paths
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p.data  =   ['/Users/coop558/mydata/mar3.11/matfiles/' siteName '/data/'];
p.save  =   ['/Users/coop558/mydata/mar3.11/matfiles/' siteName '/data/'];
p.copy  =   setpath(['GREENLAND/icemodel/model/input/userData/' siteName '/']);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% read in the catchment to get the x,y coordinate 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for n = 1:nYears
    
    thisYr  = num2str(startYear+n-1);
    fData   = [p.data 'mar_' siteName '_' thisYr '.mat'];
    fCopy   = [p.copy 'mar_' siteName '_' thisYr '.mat'];
    fSave   = fData;
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Data    = marData2Met(fData);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if saveData == true
        save(fSave,'Data');
        save(fCopy,'Data');        
    end

end
