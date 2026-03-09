function T = structToTable(s)
%STRUCTTOTABLE Robust struct-to-table conversion for saved test artifacts.

   if numel(s) == 1
      fn = fieldnames(s);
      scalarFields = true;
      n = [];
      for i = 1:numel(fn)
         v = s.(fn{i});
         if iscell(v) || isnumeric(v) || isstring(v) || isdatetime(v)
            if isscalar(v)
               continue
            end
            if isempty(n)
               n = numel(v);
            end
            if numel(v) ~= n
               scalarFields = false;
               break
            end
         else
            scalarFields = false;
            break
         end
      end
      if scalarFields && ~isempty(n)
         out = struct();
         for i = 1:numel(fn)
            v = s.(fn{i});
            if isscalar(v)
               if isstring(v) || ischar(v)
                  out.(fn{i}) = repmat(string(v), n, 1);
               elseif isnumeric(v)
                  out.(fn{i}) = repmat(v, n, 1);
               elseif isdatetime(v)
                  out.(fn{i}) = repmat(v, n, 1);
               else
                  out.(fn{i}) = repmat({v}, n, 1);
               end
            else
               out.(fn{i}) = tocol(v);
            end
         end
         T = struct2table(out);
         return
      end
   end
   T = struct2table(s);
end

function y = tocol(x)
   if ischar(x)
      y = string({x});
   elseif isstring(x) || isnumeric(x) || isdatetime(x)
      y = x(:);
   elseif iscell(x)
      y = x(:);
   else
      y = x;
   end
end
