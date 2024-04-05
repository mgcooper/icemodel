clean

%% March 2024

% Below here are notes which have been at the bottom of a_save_icemodel
% forever, so I moved them here for posterity.

% NOTES from the older ice mask

% this save the 1487 icemodel runs that cover the SW sector as matfiles.
% NOTE: because i used list=dir() to read in the data, the values are
% not stored in the order 1:1487, instead they are roughly 1, 10, 101,
% 1000, 200, 201, 2000, etc. The x,y coordinates and the data are correct,
% so as long as the x,y coordinates saved in the output files from this
% script are used for interpolation, everythign is correct, but if the x,y
% coordinates from the rcm data, or the modis ice mask, both of which are
% ordered correctly as 1:1487 then the results will be wrong. I created
% 'PATCH_ICEMASK' to deal with this (i think superceded by mask.xmodel,
% mask.ymodel coordinates in the mask file

%%% test

% load('modis_ice_mask.mat')
% elev  =  icemask.elevmodel;
% x     =  icemask.xmodel;
% y     =  icemask.ymodel;

% % this confirms I can go from x/y/elev to xmodel/ymodel/elevmodel
% xtest1 = icemask.x(logical(icemask.mask));
% ytest1 = icemask.y(logical(icemask.mask));
% xtest1 = xtest1(icemask.idxmodel);
% ytest1 = ytest1(icemask.idxmodel);
%
% xtest2 = icemask.xmask(icemask.idxmodel);
% ytest2 = icemask.ymask(icemask.idxmodel);
%
% figure; plot(xtest1,x,'o'); addOnetoOne;
% figure; plot(ytest1,y,'o'); addOnetoOne;

%% UPDATE aug 2022
% I should have also added the indices for the og mask 1
% ... 1487 relative to the met 1 ... 1487

% to explain ... if we apply the icemask.mask to the icemask.x/y we get
% 1487 points, and the met files 1...1487 are numbered in the same order.
% But when I later read in the model ouptut i used list = dir() which
% returned the files in the order 1,10,100,1000,1001,1002, and so on. So I
% made this script which collected the x,y coordinates in that order and I
% call those coordinates xmodel,ymodel. But then in aug 2022 I made the new
% gridded met files and I used xmodel/ymodel, but I should have used the
% original x/y and fixed it on the other end, i.e., not using dir(),
% because when I wanted to compare the new gridded metfiles to the old
% ones, I wasn't able to do that easily b/c they are in a different order.
% So I came back here and added the mapping b/w the dir() ordring and the
% original met ordering, see second section below

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % this is the original section where I added xmodel/ymodel, which are the
% % x/y coords in the dir() ordering
%
% load(setpath('GREENLAND/runoff/data/region/modis_ice_mask.mat'));
%
% f = 'icemodel_reference_mar_albedo_2009.mat';
% load(['/Volumes/Samsung_T5/matlab/GREENLAND/runoff/icemodel/output/region/' f]);
%
% x = data.x';
% y = data.y';
%
% icemask.xmodel = x;
% icemask.ymodel = y;
%
% figure;
% scatter(x,y,6,data.elev,'filled')
% colorbar
%
% save(setpath('GREENLAND/runoff/data/region/modis_ice_mask.mat'),'icemask');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % this is the new section where I add the linear indices that map the
% dir() ordering onto the mask (and thereofre the og met file) numbering,
% so I can build met files in list format from the new gridded met files.
% When time permits i should delete the new grids adn resave them using the
% mask ordering, or just delete them altogether since it ended up making
% more sense to save one point at a time as a ten-year met file


% the new grids are ordered in the same order as xmodel/ymodel, whereas the
% old met files are ordered in the same order as the mask, because they
% were created using x(mask), like below. The issue was that when I saved
% the output, I used list = dir() to read-in ice1, ice2, ... , ice1487

pathmask = '/Users/coop558/MATLAB/GREENLAND/runoff/data/region/';

load([pathmask 'modis_ice_mask.mat']);

% this is the syntax used in the new mk_metfile_MAR_grids but it is
% misleading b/c Xmet/Ymet should follow the icemask.mask numbering,
% whereas I assigned xmodel/ymodel (which follow dir() numbering) to
% Xmet/Ymet here and passed it to the function, but keeping the syntax for
% reference
Xmet     = icemask.xmodel;
Ymet     = icemask.ymodel;

% load the mask and get the x/y coords following the og met file numbering
mask     = logical(icemask.mask);
xmask    = icemask.x(mask);
ymask    = icemask.y(mask);

% find the mapping between them:
for n = 1:numel(Xmet)
   xn       = Xmet(n);
   yn       = Ymet(n);
   idx(n)   = find(xmask==xn & ymask==yn);
end

% add that to the mask and save it:
icemask.idxmodel  = idx(:);
icemask.xmask     = xmask;
icemask.ymask     = ymask;

icemask.readme    = ...
   ['xmodel/ymodel follow dir() numbering, xmask/ymask follow ' newline ...
   'the numbering of icemask.x(icemask.mask), which is the ' newline ...
   'numbering that the original met files follow.'];

save([pathmask 'modis_ice_mask.mat'],'icemask');

% % this is how I figured it out:
%
% idx = 22;
% p = '/Users/coop558/mydata/mar3.11/matfiles/region/level2/modis/2000/';
% load([p 'met_' int2str(idx) '.mat'])
% [Xmet(idx) met.Properties.CustomProperties.X xmask(idx)]
% [Ymet(idx) met.Properties.CustomProperties.Y ymask(idx)]
