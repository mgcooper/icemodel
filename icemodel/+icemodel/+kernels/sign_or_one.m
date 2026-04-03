function s = sign_or_one(x)
   %SIGN_OR_ONE Return the sign of X, treating zero as positive one.

   %#codegen

   if x >= 0
      s = 1;
   else
      s = -1;
   end
end
