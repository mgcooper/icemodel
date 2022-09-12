% later that same day ... moved icemodel/ out of runoff/ and renamed
% icemodel/data/preprocess/ to icemodel/preproc/data and moved
% icemodel/scripts/ to icemodel/preproc/scripts

% JUNE 2021, I need a way to immediately build a new model run, including
% met files. for ak4 I want to use kanL, which makes it somewhat more
% complicated. 

% looking around, I found build_directories.m in
% icemodel/model/versions/v10/scripts, and also:
% icemodel/model/experiment/v8/behar/scripts, but nowhere else

% I moved it here.

% I checked all other dirs there e.g. v8, v9, and there weren't any other
% obviously usefel (general) scripts that shouldn't be hidden away in
% random folders 

% I then checked everywhere else in GREENLAND/runoff/ and confirmed there
% aren't any scripts squirelled away that should be somehwere obvious, like
% here, for example, I though I might find a script to create new icemodel
% metfiles but the only one in /GREENLAND/runoff/icemodel/ is the one saved
% here

% NOTE THAT IT CALLS FUNCTIONS THAT ARE SAVED IN myFunctions/ 