clean

% this finds the point numbers for the points that have not been run, to
% read in met_X, where X is the point number

% note that in lists3, met_X is not the same as met_x in the other lists
% because the numbering goes from 1:npts, and there were 1487 points in the
% old lists, and 2479 points in the new mask

movedata =  true;
albedo   =  'modis';
pathold  =  '/Volumes/Samsung_T5b/icemodel/output/region/v10/';
pathnew  =  '/Volumes/Samsung_T5b/icemodel/output/region/v10_b/';

pathmetold = '/Volumes/Samsung_T5b/mar3.11/matfiles/region/level2/lists/';
pathmetnew = '/Volumes/Samsung_T5b/mar3.11/matfiles/region/level2/lists3/';

years = 2008:2018;
nyrs = numel(years);

load('modis_ice_mask.mat','icemask'); oldmask = icemask; clear icemask;
load('IceMask.mat','IceMask'); newmask = IceMask; clear IceMask;

% find the points that were in the old mask and new mask
xold     = oldmask.xmodel;
yold     = oldmask.ymodel;
xyold    = [xold yold];
xall     = newmask.X(newmask.IceMaskSLA);
yall     = newmask.Y(newmask.IceMaskSLA);
xyall    = [xall yall];
isnew    = ~ismember(xyall,xyold,'rows');
isold    = ismember(xyall,xyold,'rows');
isout    = ~ismember(xyold,xyall,'rows');
xnew     = xall(isnew);
ynew     = yall(isnew);
idxold   = find(isold);
idxnew   = find(isnew);
idxout   = find(isout);

noldpts  = numel(xold);
nnewpts  = numel(xnew);
nallpts  = numel(xall);
ngridpts = numel(newmask.X);

% the idxnew skips should be: 56, 57, 75, 98, 99, 100, 140, 184, 185, 245
% final file count should be 1486*2=2972 (the 20 are included in the 2972)

% n = 1; m = 1; % for testing

for m = 2:nyrs
% for m = 1:1

   fprintf('\ndoing %d \n\n',years(m))

   % load the saved data
   load([pathold 'icemodel_' albedo '_' num2str(years(m)) '.mat']);

   for n = 1:noldpts
      
      idxold = n;
      idxnew = findnewidxfromold(idxold);

      if isnan(idxnew)
         fprintf('\nskipping idxold %d (not in idxnew)\n\n',idxold)
         continue
      end
   
      % DATA CHECK - check if the x/y mapping is correct
      xyold = [data.x(idxold) data.y(idxold)];
      xynew = [xall(idxnew) yall(idxnew)];
   
      if ~isequal(xyold,xynew)
         error('xy not equal')
      end

      % MET CHECK - check that the forcing is identical

      % load the met files to double check the idx is correct
      load([pathmetold 'met_' num2str(idxold)]); metold = met; 
      load([pathmetnew 'met_' num2str(idxnew)]); metnew = met;
      metnew = rmleapinds(metnew);

      r2check(1) = corr(metold.tair,metnew.tair);
      r2check(2) = corr(metold.swd,metnew.swd);

%       figure; plot(metold.tair,metnew.tair,'o'); addOnetoOne;
%       figure; plot(metold.swd,metnew.swd,'o'); addOnetoOne;
%       figure; plot(metold.lwd,metnew.lwd,'o'); addOnetoOne;
%       figure; plot(metold.albedo,metnew.albedo,'o'); addOnetoOne;
%       metold.modis(metold.modis>1) = metold.albedo(metold.modis>1);
%       figure; plot(metold.modis,metnew.modis,'o'); addOnetoOne;
%       figure; scatter(xall,yall); hold on;
%       scatter(xyold(1),xyold(2),100,'filled');

      % at this point we've already checked that x/y are equal so relax the
      % correlation check b/c some points e.g. near the edges will not be
      % exactly equal
      if any(round(r2check,3) ~= 1)
         if n ~= 361 % 361 was checked
            error('met not equal');
         end
      end

      % create old and new file names
      fold = ['ice1_' num2str(idxold) '.mat'];
      fnew = ['ice1_' num2str(idxnew) '.mat'];
      fold = fullfile(pathold, 'icemodel', albedo, num2str(years(m)), fold);
      fnew = fullfile(pathnew, 'icemodel', albedo, num2str(years(m)), fnew);

      % check if the new file already exists - should be 10 of them
      if exist(fnew,'file')
         fprintf('skipping idxold %d (idxnew %d)\n',idxold,idxnew)
         continue;
      end

      % otherwise copy the file
      if movedata == true
         cmd = ['rsync -a --remove-source-files ' fold ' ' fnew];
         status = system(cmd);
      end

      % repeat for ice2
      fold = ['ice2_' num2str(idxold) '.mat'];
      fnew = ['ice2_' num2str(idxnew) '.mat'];
      fold = fullfile(pathold, 'icemodel', albedo, num2str(years(m)), fold);
      fnew = fullfile(pathnew, 'icemodel', albedo, num2str(years(m)), fnew);

      if movedata == true
         cmd = ['rsync -a --remove-source-files ' fold ' ' fnew];
         status = system(cmd);
      end
   end

end





