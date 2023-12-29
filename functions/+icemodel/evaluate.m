function results = evaluate(ice1, met)
   
   % varnames = commonfields()
   
   varnames = {'tsfc', 'shf', 'lhf'};
   
   for n = 1:numel(varnames)
      results.(varnames{n}) = nashsutcliffe(ice1.(varnames{n}), met.(varnames{n}));
   end
end
