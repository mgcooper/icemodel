function tf = iscomplex(X,varargin)
   %ISCOMPLEX Return true for complex number elements of X
   %
   %  Y = ISCOMPLEX(X)
   %
   % Description
   %
   %  Y = ISCOMPLEX(X) returns true/false for elements of X that have an
   %  imaginary component
   %
   % Matt Cooper, 30-Jan-2023, https://github.com/mgcooper
   %
   % See also
   %
   %#codegen

   tf = imag(X)~=0;
end
