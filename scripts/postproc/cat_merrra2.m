
% because cdo mergetime changes the dimension order, I am gonna use nco
% instead, alhtough it might be prefereable to change the dimorder, since i
% managed to figure out the ncrowcol for merra i am leaving it for now

% not sure why the shell script wouldn't work, but this works. Note that
% ncrcat "glues" the files along the record dimension, and since I first
% ran makeTimeRecDimLoop.sh, whihc specified time, it should cat correctly

clean

p.data  =   '/Users/coop558/mydata/merra2/1hrly/ncfiles/';
p.save  =   '/Users/coop558/mydata/merra2/1hrly/ncfiles/annual/';
dirList =   {'flx','glc','rad','slv'};

startYear   = 2012;
endYear     = 2018;
numYears    = endYear-startYear+1;
numDirs     = numel(dirList);
fileParts   = { 'MERRA2_400.tavg1_2d_flx_Nx.', ...
                'MERRA2_400.tavg3_2d_glc_Nx.', ...
                'MERRA2_400.tavg1_2d_rad_Nx.', ...
                'MERRA2_400.tavg1_2d_slv_Nx.'   };
            
cp1stat     = nan(numDirs,numYears);
catstat     = nan(numDirs,numYears);
cp2stat     = nan(numDirs,numYears);

for n = 1:numDirs
    
    dirData   = [p.data dirList{n} '/'];
    dirTemp   = [p.save 'scratch/'];
    
    filePrefix = fileParts{n};
    
    for m = 1:numYears
        
        thisYear = num2str(startYear + m - 1);
    
        % build the copy command
        cmdCopy  = ['cp ' dirData filePrefix thisYear '* ' dirTemp];
        
        % copy the files
        cp1stat(n,m) = system(cmdCopy);
        
        % prep for catting
        fSave   = [p.save filePrefix thisYear '.nc4'];
        
        % build the cat command
        cmdCat = ['ncrcat ' dirTemp '*.nc4 ' fSave];
        
%         cmdCat = ['cdo mergetime ' dirTemp '*.nc4 ' fSave];
        
        catstat(n,m) = system(cmdCat);
        
        % if catting successful, delete the scratch files
        if catstat(n,m) == false
            cmdRemove   = ['rm ' dirTemp '*.nc4'];
            cp2stat(n,m) = system(cmdRemove);
        end
        
    end
end


if test_data
    T       = annualCalendar(2012,hours(3));

    load('merra_grid');
    Lat     = Grid.LAT;
    Lon     = Grid.LON;

    pold    = '/Users/coop558/mydata/merra2/1hrly/ncfiles/glc/';
    pnew    = '/Users/coop558/mydata/merra2/1hrly/ncfiles/annual/';

    fnew    = [pnew 'MERRA2_400.tavg3_2d_glc_Nx.2012.nc4'];
    V       = ncread(fnew,'SNOMAS_GL');
    V       = rot90_3D(V,3,1);

    testday = '20121021';
    testidx = find(T==datetime(2012,str2num(testday(5:6)),str2num(testday(7:8))));

    fold    = [pold 'MERRA2_400.tavg3_2d_glc_Nx.' testday '.nc4.nc4'];
    vold    = ncread(fold,'SNOMAS_GL');
    vold    = rot90_3D(vold,3,1);
    vold    = vold(:,:,1);
    vnew    = V(:,:,testidx);
    dif     = vnew-vold;

    max(dif(:))
    min(dif(:))

    figure;
    subplot(1,2,1);
    geoshow(Lat,Lon,vold,'DisplayType','texturemap');
    subplot(1,2,2);
    figure; geoshow(Lat,Lon,vnew,'DisplayType','texturemap');

end


% check results:
% fn      = 'MERRA2_400.tavg1_2d_flx_Nx.2012.nc4';
% pcdo    = '/Users/coop558/mydata/merra2/1hrly/ncfiles/annual/bk/';
% pnco    = '/Users/coop558/mydata/merra2/1hrly/ncfiles/annual/';
% 
% fcdo    = [pcdo fn];
% fnco    = [pnco fn];
% 
% icdo    = ncinfo(fcdo);
% inco    = ncinfo(fnco);
% 
% vcdo    = ncread(fcdo,'EVAP');
% vnco    = ncread(fnco,'EVAP');
% 
% figure; 
% for n = 1:size(vcdo,1)
%     for m = 1:size(vcdo,2)
%         myscatter(squeeze(vcdo(n,m,:)),squeeze(vnco(n,m,:))); addOnetoOne;
%         pause;
%     end
% end
% 
% dataNew = ncreaddata(fSave);
% 
% 
% 
% % % % % % % % % % % 
% %% initial tests
% % % % % % % % % % % 
% 
% 
% fcdo = '/Users/coop558/mydata/merra2/1hrly/ncfiles/test/test_cdo.nc4';
% fnco = '/Users/coop558/mydata/merra2/1hrly/ncfiles/test/test_nco.nc4';
% 
% infocdo = ncinfo(fcdo);
% infonco = ncinfo(fnco);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
