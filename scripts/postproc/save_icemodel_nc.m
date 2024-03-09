clean

% There was a version of this script called _test that tested only writing the
% 2d vars (and maybe only the 0-d vars based on the note below but the note may
% be a mistake, it seems it tested only writing 2-d vars). In the test version,
% comprsz was set to 9, so I added that to the if test_compression == true
% section, but I think that was just a one-off setting, separate from the actual
% test of including 0-d, 1-d, or 2-d vars.

savedata = true;
sitename = 'sector';
simmodel = 'icemodel'; % {'icemodel', 'skinmodel'};
userdata = 'mar'; % {'mar', 'modis'};
siteopts = setBasinOpts('sitename', sitename, 'simmodel', simmodel, 'userdata', userdata);

make_backups = true;

% I also tested if compression is better w/o the 1-d and 2-d vars - it isn't
test_compression = true;
if test_compression == true
   comprsz = 9;           % compression size
else
   comprsz = 1;           % compression size
end
% file size in mb/cell for min/med/max compression
% 1: 10.9 (10.8 w/o 1-d or 0-d vars)
% 5: 9.4
% 9: 9.2 (9.1 w/o 1-d or 0-d vars)


% Set path to data
pathdata = fullfile(getenv('ICEMODELOUTPUTPATH'), sitename);
pathsave = fullfile(getenv('ICEMODELOUTPUTPATH'), sitename, 'nc');
pathmet  = fullfile(getenv('ICEMODELINPUTPATH'), 'met', sitename);

%%

icemodel.netcdf.makencfile(pathdata, pathsave, siteopts)
