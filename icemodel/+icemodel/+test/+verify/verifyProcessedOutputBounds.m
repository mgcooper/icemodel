function verifyProcessedOutputBounds(testCase, ice1, ice2)
   %VERIFYPROCESSEDOUTPUTBOUNDS Verify basic physical bounds of processed output.
   %
   %  icemodel.test.verify.verifyProcessedOutputBounds(testCase, ice1, ice2)

   arguments
      testCase
      ice1 timetable
      ice2 struct
   end

   % Processed outputs should always contain finite retained timesteps.
   testCase.verifyGreaterThan(height(ice1), 0);
   testCase.verifyTrue(all(isfinite(ice1{:, :}), 'all'), ...
      'processed ice1 output contains non-finite values');

   if isfield(ice2, 'Time')
      testCase.verifyEqual(numel(ice2.Time), height(ice1));
   end

   % Check phase-fraction bounds using the same density convention as the
   % core column model.
   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   if isfield(ice2, 'f_ice')
      testCase.verifyGreaterThanOrEqual(min(ice2.f_ice(:)), -1e-9);
      testCase.verifyLessThanOrEqual(max(ice2.f_ice(:)), 1 + 1e-9);
   end
   if isfield(ice2, 'f_liq')
      testCase.verifyGreaterThanOrEqual(min(ice2.f_liq(:)), -1e-9);
      testCase.verifyLessThanOrEqual(max(ice2.f_liq(:)), 1 + 1e-9);
   end
   if isfield(ice2, 'f_ice') && isfield(ice2, 'f_liq')
      massfrac = ice2.f_ice + ice2.f_liq * ro_liq / ro_ice;
      testCase.verifyLessThanOrEqual(max(massfrac(:)), 1 + 1e-6);
   end

   % Retained cumulative diagnostics should remain monotone after processing.
   if ismember('runoff', ice1.Properties.VariableNames)
      testCase.verifyGreaterThanOrEqual(min(diff(ice1.runoff)), -1e-6);
   end
   if ismember('melt', ice1.Properties.VariableNames)
      testCase.verifyGreaterThanOrEqual(min(diff(ice1.melt)), -1e-6);
   end
   if ismember('tsfc', ice1.Properties.VariableNames)
      testCase.verifyLessThanOrEqual(max(ice1.tsfc), 1e-6);
   end
end
