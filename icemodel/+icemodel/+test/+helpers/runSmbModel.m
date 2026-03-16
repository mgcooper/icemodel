function varargout = runSmbModel(opts)
%RUNSMBMODEL Dispatch to the requested core SMB model kernel.
%
%  [ice1, ice2] = icemodel.test.helpers.runSmbModel(opts)
%  [ice1, ice2, opts] = icemodel.test.helpers.runSmbModel(opts)

   switch opts.smbmodel
      case 'icemodel'
         [varargout{1:nargout}] = icemodel(opts);
      case 'skinmodel'
         [varargout{1:nargout}] = skinmodel(opts);
      otherwise
         error('unsupported smbmodel: %s', opts.smbmodel)
   end
end
