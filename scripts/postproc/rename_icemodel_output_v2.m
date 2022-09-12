clean

p = '/Users/coop558/mydata/icemodel/output/upperBasin/v8/';

list = getlist(p,'*.mat');

for n = 1:numel(list)
   
   f = [p list(n).name];
   
   fnew = strrep(f,'albedo_','');
%    fnew = strrep(fnew,'reference_','');
   
   copyfile(f,fnew);
   
   delete(f);
end