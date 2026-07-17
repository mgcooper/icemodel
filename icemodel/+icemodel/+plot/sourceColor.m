function colors = sourceColor(names)
   %SOURCECOLOR Return stable verification colors keyed by source identity.
   %
   %  colors = icemodel.plot.sourceColor(names)
   %
   % Names may include plot roles such as "promice met" or
   % "mar3.11 userdata". MAR, RACMO, and MERRA-2 use the established
   % runoff/RunoffPlot palette. Other known verification sources use durable
   % colors across roles, while unknown labels receive a deterministic fallback.

   names = reshape(string(names), [], 1);
   colors = zeros(numel(names), 3);
   fallback = [ ...
      0.0000 0.4470 0.7410
      0.6350 0.0780 0.1840
      0.3010 0.7450 0.9330
      0.4660 0.6740 0.1880
      0.8500 0.3250 0.0980];

   for k = 1:numel(names)
      key = lower(strtrim(names(k)));
      if contains(key, "merra")
         colors(k, :) = [0 0 0.172413793103448];
      elseif contains(key, "racmo")
         colors(k, :) = [0.25 0.80 0.54];
      elseif contains(key, "mar3") || key == "mar" ...
            || startsWith(key, "mar ")
         colors(k, :) = [0.866 0.329 0];
      elseif contains(key, "promice")
         colors(k, :) = [0.4940 0.1840 0.5560];
      elseif contains(key, "modis")
         % Keep the observational albedo overlay neutral but visibly distinct
         % from the near-black MERRA-2 source color.
         colors(k, :) = [0.55 0.55 0.55];
      elseif startsWith(key, "observations")
         colors(k, :) = [0 0 0];
      elseif contains(key, "sumup")
         % Keep sparse SUMup observations distinct from the fixed RCM colors.
         colors(k, :) = [0.0000 0.4470 0.7410];
      elseif contains(key, "esm_snowmip") || contains(key, "esm-snowmip") ...
            || contains(key, "esmsnowmip")
         colors(k, :) = [0.0000 0.4470 0.7410];
      elseif contains(key, "laugh") || contains(key, "colbeck")
         colors(k, :) = [0.6350 0.0780 0.1840];
      else
         % Unknown sources keep one color across roles and variable suffixes.
         base_key = regexprep(key, "\s+-\s+[^-]+$", "");
         base_key = regexprep(base_key, ...
            "\s+(forcing|observations|met|userdata)$", "");
         code = double(char(base_key));
         index = mod(sum(code .* (1:numel(code))), size(fallback, 1)) + 1;
         colors(k, :) = fallback(index, :);
      end
   end
end
