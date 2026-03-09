function S = summarizeIce1Metrics(ice1)
%SUMMARIZEICE1METRICS Extract core regression metrics from ice1 output.
%
%  S = test.helpers.summarizeIce1Metrics(ice1)
%
% Output fields:
%  runoff_final, melt_final, mean_Tice_numiter, max_Tice_numiter,
%  n_not_converged

   S = struct( ...
      'runoff_final', nan, ...
      'melt_final', nan, ...
      'mean_Tice_numiter', nan, ...
      'max_Tice_numiter', nan, ...
      'n_not_converged', nan);

   if istimetable(ice1)
      vn = string(ice1.Properties.VariableNames);
   elseif isstruct(ice1)
      vn = string(fieldnames(ice1));
   else
      error('ice1 must be timetable or struct')
   end

   if any(vn == "runoff")
      r = getField(ice1, 'runoff');
      if ~isempty(r)
         if size(r, 2) > 1
            r = r(:, 1);
         end
         S.runoff_final = r(end);
      end
   end

   if any(vn == "melt")
      m = getField(ice1, 'melt');
      if ~isempty(m)
         if size(m, 2) > 1
            m = m(:, 1);
         end
         S.melt_final = m(end);
      end
   end

   if any(vn == "Tice_numiter")
      niter = getField(ice1, 'Tice_numiter');
      niter = niter(isfinite(niter));
      if ~isempty(niter)
         S.mean_Tice_numiter = mean(niter);
         S.max_Tice_numiter = max(niter);
      end
   end

   if any(vn == "Tice_converged")
      conv = logical(getField(ice1, 'Tice_converged'));
      S.n_not_converged = sum(~conv);
   end
end

function x = getField(s, name)
   x = s.(name);
end
