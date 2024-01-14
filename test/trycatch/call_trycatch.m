function call_trycatch(a, b, option)
   if nargin == 2
      option = 1;
   end
   if option == 1
      trycatch_w_handling(a, b);
   elseif option == 2
      trycatch_wo_handling(a, b);
   end
end
