function returned = captureExpectedWarning(testCase, fcn, expected_id)
   %CAPTUREEXPECTEDWARNING Run fcn once, verify its warning, and capture output.
   if ~isa(fcn, 'function_handle')
      error('icemodel:test:captureExpectedWarning:badInput', ...
         'fcn must be a function handle')
   end

   % Capture command-window output so expected warning text does not clutter
   % test logs while still leaving the warning id available for verification.
   lastwarn('');
   [~, returned] = evalc('fcn();');
   [~, actual_id] = lastwarn();
   testCase.verifyEqual(actual_id, expected_id);
end
