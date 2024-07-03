function tf = isoctave()
   %ISOCTAVE Return logical true if the environment is Octave otherwise false.
   %
   %  tf = isoctave()
   %
   % See also: ishg2, isinteractive

   persistent cacheval;
   if isempty (cacheval)
      cacheval = (exist ("OCTAVE_VERSION", "builtin") > 0);
   end
   tf = cacheval;
end
