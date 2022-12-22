clean

fname    = 'solar.dat'; % output file name
savedata = false;
plotdata = false;

% set paths
pathdata = setpath('runoff/data/icemodel/spectral/','project');
pathsave = setpath('runoff/data/icemodel/spectral/','project');

%% read in the solar spectrum
load(fullfile(pathdata,'Qi_coart.mat'));
QiCoart = Qi_coart.Qi;
wavlCoart = Qi_coart.Lambda;

% read in Glen's spectrum
solar = importdata(fullfile(pathdata,'solar.dat'));
wavl = 1000.*solar(:,1);
Qitmp = solar(:,2);

% interpolate my spectrum to his lambda
Qi = interp1(wavlCoart,QiCoart,wavl);
Qi(isnan(Qi)) = Qitmp(isnan(Qi)); % fill in missing lambda values 

%% make figure

if plotdata == true
    figure(1)
    p(1) = plot(wavl,solar(:,2)); hold on;
    p(2) = plot(wavlCoart,QiCoart);
    p(3) = plot(wavl,Qi);
    l    = legend('Antarctica','Greenland - OG','Greenland - Interp');
    if savefigs == true
        exportgraphics(gcf,fullfile(pathsave,'Qi_incoming.png'),'Resolution','400')
    end
end

%% save the output
data  = transpose(roundn([wavl./1000 Qi],-5));
pwavl = '%8.5f';
pQi   = '%12.5f\n';

if savedata == 1
   fid = fopen(fullfile(pathsave,fname),'w'); 
   fprintf(fid,[pwavl pQi],data);
   fclose(fid);
end
