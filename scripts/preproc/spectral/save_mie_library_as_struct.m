% i just wrote this to save the mie library in a more accessible structure
% format

clean

% gett the radii, they aren't in the cleaned final icemodel version
path.mie    =   ['/Volumes/Samsung_T5/matlab/GREENLAND/runoff/icemodel/' ...
                    'data/spectral/mie_library_v3/'];
load([path.mie 'mie_library_big_all_vars.mat'],'radii');

% this one has the final values
path.mie    =   setpath('GREENLAND/runoff/icemodel/data/preprocess/spectral/');
load([path.mie 'mie.mat']);
lambda      =   mie(end,:);

% make a new empty 'mie' structure
mie_old     =   mie; clear mie
rs          =   length(radii);          % starting index

mie.radii   = radii';
mie.lambda  = lambda.*1000;
mie.Qext    = mie_old(rs+1:2*rs,:);
mie.g       = mie_old(1:rs,:);
mie.w       = 1-mie_old(2*rs+1:3*rs,:);

%%
figure; 
plot(mie.lambda,mie.Qext(1,:)); hold on; 
plot(mie.lambda,mie.Qext(end,:))

figure; 
plot(mie.lambda,mie.g(1,:)); hold on; 
plot(mie.lambda,mie.g(end,:))

figure; 
plot(mie.lambda,mie.w(1,:)); hold on; 
plot(mie.lambda,mie.w(end,:))

%%
save([path.mie 'mie_library_v3'],'mie');
%% 

% this is how the indexing was done in plot_mie_libary, so i followed the
% pattern to pull out the data above
% for n = 1:length(radii)    
%     mie.g(:,n)      = mie(n,:);
%     mie.Qext(:,n)   = mie(1*rs+n,:); % qext
%     mie.w(:,n)      = 1-mie(2*rs+n,:); % 1-w
% end
    
    

