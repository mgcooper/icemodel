function path = setpath(path)

% nov 2021 learned about userpath, using this instead of old method
    
if isstruct(path)
    fields = fieldnames(path);

    for n = 1:length(fields)
        pathname = [getenv('MATLABUSERPATH') path.(fields{n})];
        path.(fields{n}) = pathname;
    end
    
elseif ischar(path)
    path = [getenv('MATLABUSERPATH') path];
end

    
% the old way that required knowing which computer I was on:
    
% homepath = pwd;
% 
% if isstruct(path)
%     fields = fieldnames(path);
% 
%     for n = 1:length(fields)
% 
%         pathname = path.(fields{n});
% 
%         if strcmp(homepath(8:14),'coop558')
%             pathname = ['/Users/coop558/Dropbox/MATLAB/' pathname];
%         elseif strcmp(homepath(8:17),'mattcooper')
%             pathname = ['/Users/mattcooper/Dropbox/MATLAB/' pathname];
%         end
% 
%         path.(fields{n}) = pathname;
%     end
%     
% elseif ischar(path)
%     
%     if strcmp(homepath(8:14),'coop558')
%             path = ['/Users/coop558/Dropbox/MATLAB/' path];
%     elseif strcmp(homepath(8:17),'mattcooper')        
%             path = ['/Users/mattcooper/Dropbox/MATLAB/' path];
%     end
%         
% end
% 
%     
% 
