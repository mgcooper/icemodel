function forcing_sources = normalizeForcingSources(values, build_forcing)
   %NORMALIZEFORCINGSOURCES Normalize one public forcing-source selection.
   %
   %  forcing_sources = ...
   %     icemodel.verification.setup.normalizeForcingSources( ...
   %     values, build_forcing)

   % Public selectors are ordered sets; blank tokens do not name artifacts.
   forcing_sources = reshape(string(values), 1, []);
   forcing_sources = forcing_sources(strlength(strtrim(forcing_sources)) > 0);
   forcing_sources = unique(forcing_sources, 'stable');

   % A true build toggle without a requested artifact is contradictory.
   if build_forcing && isempty(forcing_sources)
      error('icemodel:verification:normalizeForcingSources:emptySelection', ...
         ['forcing_sources cannot be empty when build_forcing=true. ' ...
         'Use build_forcing=false for observation-only staging.']);
   end
end
