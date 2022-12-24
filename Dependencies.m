clean

%%
refpath = fullfile(getuserenvs('MATLABPROJECTPATH'),'icemodel','functions');
Info = getFunctionDependencies('icemodel','refpath',refpath);

%% cycle through the dependent functions and copy them to util/

% NOTE: this is for private use, it won't work if you don't have the functions
% on your local computer. Please contact me at matt.cooper@pnnl.gov if you have
% any trouble running this toolbox or if any function dependencies are missing
% from the toolbox. Alternatively, look for the missing functions in
% https://github.com/mgcooper/matfunclib (I suggest the dev branch). Thank you. 

% TODO: add method to clone from https://github.com/mgcooper/matfunclib

for n = 1:numel(Info.fmissing)
   [~,fname,ext] = fileparts(Info.fmissing{n});
   copyfile(Info.fmissing{n},['scripts/util/',fname,ext]);
end


% where to pick up 
% i think skinmodel is integrated into icemodel repo and I started to confirm
% that skinmodel replicates the emulators, at first i got it wrong b/c dt = 900
% and it has to = 3600 to replicate. I still need to check which functions are
% in skinmodel python feature and copy into icemodel then move to private, I
% planned to do that but forgot and committed getnewt, thenI need to check
% skinmodel main and dev and decide if i want to delete skinmodel altogether,
% also check the differences in the opts and run scripts, in particular the
% ensemble vs normal runs
