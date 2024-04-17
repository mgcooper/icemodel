
%% runscript settings

% If sitename is 'sector', specify which grid point to run
runpoint = 239; % 2195;

% % These were the settings I had to rerun
% sitename = 'sector';        % options: 'kanm', 'behar'
% forcings = 'mar';         % options: 'mar','kanm'
% userdata = 'modis';          % options: 'modis','racmo','merra','mar','kanm','none'
% simyears = 2013:2013;
% runpoint = 1794;

%% rungrid settings

testname = 'zobs2';
% Option to run specific grid points:
gridnums = 2195;
% gridnums = [368, 2337]; % bad mar pixels
% gridnums = [2320]; % bad modis pixel

icemodel.run.grid(gridnums=gridnums, testname=testname)





