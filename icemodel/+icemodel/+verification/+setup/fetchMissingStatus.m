function status = fetchMissingStatus(missing)
   %FETCHMISSINGSTATUS Convert missing source patterns to fetch status rows.
   %
   %  status = icemodel.verification.setup.fetchMissingStatus(missing)
   %
   % Missing-list fetchers use this adapter before delegating to the shared
   % finishFetchStatus completion contract.

   missing = reshape(string(missing), 1, []);
   row = struct('product', "", 'present', false);
   status = repmat(row, 1, numel(missing));
   for k = 1:numel(missing)
      status(k).product = missing(k);
   end
end
