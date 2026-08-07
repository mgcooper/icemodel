function [ice1, ice2, opts] = runModelCase(c, varargin)
   %RUNMODELCASE Resolve, execute, and postprocess one supported model case.
   %
   %  [ice1, ice2, opts] = icemodel.test.helpers.runModelCase(c)
   %  [ice1, ice2, opts] = icemodel.test.helpers.runModelCase( ...
   %     case_manifest, startdate=..., enddate=...)
   %
   % C accepts the same formal-row or verification-manifest contract as
   % setModelOptsForCase. Additional inputs are forwarded unchanged.

   % Keep case resolution, production dispatch, and canonical postprocessing
   % together so verification and formal regression execute the same path.
   opts = icemodel.test.helpers.setModelOptsForCase(c, varargin{:});
   [ice1, ice2, opts] = icemodel.test.helpers.runSmbModel(opts);
   [ice1, ice2] = icemodel.postprocess( ...
      ice1, ice2, opts, opts.output_years);
end
