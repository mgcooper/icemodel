clean

checkFirst  = false;

siteName = 'upperBasin';
cd(['/Users/coop558/MATLAB/GREENLAND/icemodel/output/' siteName]);

list     = getlist([pwd '/'],'mat');
dirName  = [list(1).folder '/'];
dirSave  = strrep(dirName,[siteName '/'],'');

uniqueVals  = {'swap_KANL','swap_MERRA','swap_MODIS','swap_NONE','swap_RACMO'};
uniqueYears = 2016;
numYears    = numel(uniqueYears);

sumTotal    = numel(list);
sumFiles    = 0;


for n = 1:numel(uniqueVals)
   
   idx      = find(contains({list.name},uniqueVals{n}));
   thisList = string({list(idx).name})';
   numFiles = numel(thisList);
   
   sumFiles = sumFiles + numFiles;
   
   if sumFiles > sumTotal
      break;
   end
   
   for m = 1:numYears
         
      thisYear    = string(uniqueYears(m));
      iyear       = find(contains(thisList,thisYear));
      thisYearList = thisList(iyear);
      
      % for each year, there should be four files
      if numel(iyear) == 0
         continue;
      end
      
      if numel(iyear) ~= 4
         error('check file list');
      end
      
      % build a new filename to save the output
      thisYear = char(thisYear);
      newfile  = char(thisYearList(contains(thisYearList,'enbal')));
      newfile  = strrep(newfile,'enbal_','');
      newfile  = strrep(newfile,['_' thisYear],'');
      newfile  = strrep(newfile,'icemodel_',['icemodel_' siteName '_' thisYear '_']);

      newfile  = [dirSave newfile];
      
      
      % USE THIS TO CONFIRM BEFORE SAVEING
      disp(thisList(iyear(1)))
      disp(thisList(iyear(2)))
      disp(thisList(iyear(3)))
      disp(thisList(iyear(4)))
      disp(newfile); 
         
      if checkFirst == true
      	pause;
      else
         % load the four files and resave
         load(thisList(iyear(1)))
         load(thisList(iyear(2)))
         load(thisList(iyear(3)))
         load(thisList(iyear(4)))
         save(newfile,'enbal','ice1','ice2','opts');
         clear newfile enbal ice1 ice2 opts
      end
               
   end
     
end

% 
% % rmStrs   = {'MARalbedo_','KANLalbedo_'};
% rmStrs   = {'forcing_MARalbedo'};
% 
% for n = 1:numel(list)
%       
%    oldName  = list(n).name;
%    
%    for m = 1:numel(rmStrs)
%          
%       if contains(oldName,rmStrs{m})
%          newName  = strrep(oldName,rmStrs{m},'');
%       end
%    end
%    
%    movefile([dirName oldName],[dirName newName]);
%    
% end
% 
% 
% % % % % % % % % % % 
% 
% list     = getlist([pwd '/'],'mat');
% dirName  = [list(1).folder '/'];
% 
% % rmStrs   = {'MARalbedo_','KANLalbedo_'};
% rmStrs   = {'forcing_MARalbedo'};
% 
% for n = 1:numel(list)
%       
%    oldName  = list(n).name;
%    
%    for m = 1:numel(rmStrs)
%          
%       if contains(oldName,rmStrs{m})
%          newName  = strrep(oldName,rmStrs{m},'');
%       end
%    end
%    
%    movefile([dirName oldName],[dirName newName]);
%    
% end
% 
% % update the list
% list     = getlist([pwd '/'],'mat');
% varNames    = {'enbal','ice1','ice2','opts'};
% 
% for n = 1:numel(list)
%    
%    oldName  = list(n).name;
%    
%    for m = 1:numel(varNames)
%    
%       thisVar  = varNames{m};
%          
%       if contains(oldName,thisVar)
%          newName = strrep(oldName,[thisVar '_icemodel_'],'');
%          newName = ['icemodel_' thisVar '_' newName];
%       end
%       
%    end
%    
%    movefile([dirName oldName],[dirName newName]);
%    
% end
% 
% % the original:
% % varNames    = {'enbal','ice1','ice2','opts'};
% % 
% % for n = 1:numel(list)
% %    
% %    oldName  = list(n).name;
% %    
% %    for m = 1:numel(varNames)
% %    
% %       thisVar  = varNames{m};
% %          
% %       if contains(oldName,thisVar)
% %          newName = strrep(oldName,['_' thisVar],'');
% %          newName = [thisVar '_' newName];
% %       end
% %       
% %    end
% %    
% %    movefile([dirName oldName],[dirName newName]);
% %    
% % end
%    