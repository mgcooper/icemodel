function out = tocellstr(in)
   %TOCELLSTR Convert char and string arrays to cellstr, or return cellstr input
   %
   %  OUT = TOCELLSTR(IN) Converts string, cell, or char array IN to cellstr
   %  array OUT
   %
   % See also: tocolumn, torow, rowvec
   %
   %#codegen

   %#ok<*ISCLSTR>
   validateattributes(in, {'string', 'cell', 'char'}, {'vector'}, ...
      mfilename, 'IN', 1)
   out = ifelse(iscellstr(in), in, cellstr(in));
end
