function target = copyFields(target, source, only_matching)
   %COPYFIELDS Copy SOURCE struct fields onto TARGET.
   % ONLY_MATCHING keeps TARGET's field set.
   %
   % See also: icemodel.buildOutputPayload,
   %  icemodel.couplers.update_solver_diag,
   %  icemodel.forcing.helpers.promiceShortwave
   %
   %#codegen

   % Two-argument calls copy every source field.
   if nargin < 3
      only_matching = false;
   end
   fields = fieldnames(source);

   % ONLY_MATCHING keeps the target field set unchanged.
   for k = 1:numel(fields)
      if ~only_matching || isfield(target, fields{k})
         target.(fields{k}) = source.(fields{k});
      end
   end
end
