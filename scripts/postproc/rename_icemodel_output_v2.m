clean

% p = '/Volumes/Samsung_T5b/icemodel/output/leverett/v8c/';
p = '/Users/coop558/mydata/icemodel/output/upperBasin/v8/';

list = getlist(p,'*.mat');

for n = 1:numel(list)
   f = [p list(n).name];
   fnew = strrep(f,'albedo_','');
   copyfile(f,fnew);
   delete(f);
end

% fnew = strrep(fnew,'reference_','');