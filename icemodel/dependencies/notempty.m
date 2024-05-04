function tf = notempty(x)
   %NOTEMPTY Determine whether array X contains any non-empty elements
   %
   %  TF = NOTEMPTY(X)
   %
   % Description
   %  TF = NOTEMPTY(X) returns true if ~all(isempty(X)) is true
   %
   % Matt Cooper, 27-Oct-2022, https://github.com/mgcooper
   %
   % See also any, isempty, notempty
   %
   %#codegen

   tf = ~all(isempty(x));
end
