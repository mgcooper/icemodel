function tf = iscomplex(X,varargin)
   %ISCOMPLEX Return true for complex number elements of X
   %
   % Syntax
   %
   %  Y = ISCOMPLEX(X) returns true/false for elements of X that have an
   %  imaginary component
   %
   %
   % Matt Cooper, 30-Jan-2023, https://github.com/mgcooper
   %
   % See also
   %
   %#codegen

   tf = imag(X)~=0;

   % % got confused and was thinking of sending back [tfx,tfy] versus tf(:), but
   % the parsing below is relevant to those cases or if find() is used

   % switch nargout
   %    case 1
   %       varargout{1} = tf;
   %    case 2
   % end
   %
   % if isvector(X)
   %    tf = imag(X)~=0;
   % elseif ismatrix(X) % vectors are matrices but this will catch mxn with m,n>1
   %    tf = imag(X)~=0;
   % end
end
