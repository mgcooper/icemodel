function tf = isfullfile(filename)
    %ISFULLFILE Check if a string represents a full file path.
    %
    % tf = isfullfile(filename) returns true if the string filename
    % appears to represent a full file path on the system,
    % and false otherwise.
    %
    % Example:
    %   tf = isfullfile('C:\Users\file.txt'); % returns true on Windows
    %   tf = isfullfile('/home/user/file.txt'); % returns true on Linux/Mac
    %   tf = isfullfile('file.txt'); % returns false
    %
    % See also: ismfile, isonpath, isfullpath

    [allparts{1:3}] = fileparts(filename);
    tf = all(~cellfun(@isempty, allparts));

    if ispc
        tf = tf && ~isempty(regexp(allparts{1}, '^[a-zA-Z]:\\', 'once'));
    else
        tf = tf && startsWith(allparts{1}, '/');
    end
end
