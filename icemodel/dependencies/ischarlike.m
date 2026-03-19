function tf = ischarlike(x, varargin)
   %ISCHARLIKE Compatibility wrapper for ISTEXTLIKE.
   %
   %  TF = ISCHARLIKE(X, ...) forwards to ISTEXTLIKE(X, ...).
   %
   %  The historical name ISCHARLIKE is preserved for backward compatibility.
   %  New code should prefer ISTEXTLIKE when the intent is "contains only text
   %  values" rather than specifically "is char-like".
   %
   % See also ISTEXTLIKE

   %#codegen

   tf = istextlike(x, varargin{:});
end
