function mustBeStruct(S, function_name, variable_name, argument_position)
   %MUSTBESTRUCT Verifies that the input is a structure
   %
   %  MUSTBESTRUCT(S, FUNCTION_NAME, VARIABLE_NAME, ARGUMENT_POSITION)
   %
   % Description
   %   MUSTBESTRUCT(S, FUNCTION_NAME, VARIABLE_NAME, ARGUMENT_POSITION)
   %   verifies that S is a structure.  If it isn't, MUSTBESTRUCT issues an
   %   error message using FUNCTION_NAME, VARIABLE_NAME, and
   %   ARGUMENT_POSITION.  FUNCTION_NAME is the name of the user-level
   %   function that is checking the struct, VARIABLE_NAME is the name of the
   %   struct variable in the documentation for that function, and
   %   ARGUMENT_POSITION is the position of the input argument to that
   %   function.
   %
   %#codegen

   % Based on checkstruct
   %   Copyright 1996-2011 The MathWorks, Inc.

   if ~isstruct(S)
      error(['MATFUNCLIB:LIBSTRUCT' function_name, ':expectedStructInput'], ...
         'Function %s expected input number %d, %s, to be a struct.', ...
         upper(function_name), argument_position, upper(variable_name))
   end
end
