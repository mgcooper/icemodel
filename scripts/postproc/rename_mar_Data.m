clean

% this loads the data files and renames them 'Data' from 'data'

% dataPath    = '/Users/coop558/mydata/mar3.11/matfiles/';
dataPath    = '/Users/coop558/MATLAB/GREENLAND/icemodel/model/input/external/';

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
siteName    = 'ak4';
firstYear   = 2009;
lastYear    = 2016;

% cd([dataPath siteName '/data/']);
cd([dataPath siteName '/']);

for n = firstYear:lastYear
    
    load(['mar_' siteName '_' int2str(n) '.mat'],'data')
    Data = data;
    save(['mar_' siteName '_' int2str(n) '.mat'],'Data')
    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
siteName    = 'slv1';
firstYear   = 2015;
lastYear    = 2015;

% cd([dataPath siteName '/data/']);
cd([dataPath siteName '/']);

for n = firstYear:lastYear
    
    load(['mar_' siteName '_' int2str(n) '.mat'],'data')
    Data = data;
    save(['mar_' siteName '_' int2str(n) '.mat'],'Data')
    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
siteName    = 'slv2';
firstYear   = 2015;
lastYear    = 2015;

% cd([dataPath siteName '/data/']);
cd([dataPath siteName '/']);

for n = firstYear:lastYear
    
    load(['mar_' siteName '_' int2str(n) '.mat'],'data')
    Data = data;
    save(['mar_' siteName '_' int2str(n) '.mat'],'Data')
    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% siteName    = 'L41A';
% firstYear   = 2012;
% lastYear    = 2012;
% 
% % cd([dataPath siteName '/data/']);
% cd([dataPath siteName '/']);
% 
% for n = firstYear:lastYear
%     
%     load(['mar_' siteName '_' int2str(n) '.mat'],'data')
%     Data = data;
%     save(['mar_' siteName '_' int2str(n) '.mat'],'Data')
%     
% end



