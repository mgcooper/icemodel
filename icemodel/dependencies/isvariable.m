function [tf, vi] = isvariable(varname, T)
   %ISVARIABLE Determine if VARNAME is a variable in table T.
   %
   %  [TF, VI] = ISVARIABLE(VARNAME,T) returns TF = true if VARNAME is a
   %  variable in table T, and the variable (column) index VI.
   %
   % See also: ismember, table
   %
   %#codegen

   arguments
      varname (:,1) string
      T (:,:) tabular
   end
   tf = any(strcmp(varname, T.Properties.VariableNames));
   vi = find(strcmp(varname, T.Properties.VariableNames));
end
