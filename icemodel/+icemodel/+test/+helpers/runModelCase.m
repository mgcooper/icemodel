function [ice1, ice2, opts] = runModelCase(c, varargin)
   %RUNMODELCASE Resolve, execute, and postprocess one supported model case.
   %
   %  [ice1, ice2, opts] = icemodel.test.helpers.runModelCase(c)
   %  [ice1, ice2, opts] = icemodel.test.helpers.runModelCase( ...
   %     case_manifest, startdate=..., enddate=...)
   %
   % C accepts the same formal row or verification manifest as
   % setModelOptsForCase. Additional inputs are forwarded unchanged.

   % Verification and formal regression both execute this path: case
   % resolution, production dispatch, then canonical postprocessing.
   opts = icemodel.test.helpers.setModelOptsForCase(c, varargin{:});
   [ice1, ice2, opts] = icemodel.test.helpers.runSmbModel(opts);
   [ice1, ice2] = icemodel.postprocess( ...
      ice1, ice2, opts, opts.output_years);
end
