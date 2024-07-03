function [tf, vi] = isvariable(varname, tbl)
   %ISVARIABLE Determine if VARNAME is a variable in table TBL.
   %
   %  [TF, VI] = ISVARIABLE(VARNAME, TBL) returns TF = true if VARNAME is a
   %  variable in table TBL, and the variable (column) index VI.
   %
   % See also: ismember, table
   %
   %#codegen

   % arguments
   %    varname (:,1) string
   %    T (:,:) tabular
   % end
   validateattributes(varname, {'char'}, {'row'}, mfilename, 'VARNAME', 1)
   validateattributes(tbl, {'tabular'}, {'nonempty'}, mfilename, 'TBL', 1)

   tf = any(strcmp(varname, tbl.Properties.VariableNames));
   vi = find(strcmp(varname, tbl.Properties.VariableNames));
end
