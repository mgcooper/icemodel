function [ice1, ice2] = updateoutput(step, ice1, ice2, vars1, vars2, data1, data2)
   %UPDATEOUTPUT Store one timestep of model output into the ice1/ice2 structs.
   %
   %  [ice1, ice2] = icemodel.updateoutput(step, ice1, ice2, ...
   %     vars1, vars2, data1, data2)
   %
   % Writes each element of data1 (scalar) into ice1.(vars1{n})(step, 1) and
   % each element of data2 (column vector) into ice2.(vars2{n})(:, step).
   % The order of vars1/vars2 must match the order of data1/data2 as assembled
   % by icemodel.buildOutputPayload.
   %
   %#codegen

   for n = 1:numel(data1)
      ice1.(vars1{n})(step, 1) = data1{n};
   end

   for n = 1:numel(data2)
      ice2.(vars2{n})(:, step) = data2{n};
   end
end
