function verifyEqualNested(testCase, A, B, tol)
%VERIFYEQUALNESTED Recursively compare nested structs, timetables, and arrays.
%
%  icemodel.test.verify.verifyEqualNested(testCase, A, B, tol)

   if istimetable(A) && istimetable(B)
      verifyTimetableEqual(testCase, A, B, tol);
   elseif isstruct(A) && isstruct(B)
      verifyStructEqual(testCase, A, B, tol);
   else
      verifyArrayEqual(testCase, A, B, tol);
   end
end

function verifyTimetableEqual(testCase, A, B, tol)

   testCase.verifyEqual(A.Time, B.Time);
   testCase.verifyEqual(string(A.Properties.VariableNames), ...
      string(B.Properties.VariableNames));

   names = string(A.Properties.VariableNames);
   for i = 1:numel(names)
      verifyArrayEqual(testCase, A.(char(names(i))), B.(char(names(i))), tol);
   end
end

function verifyStructEqual(testCase, A, B, tol)

   testCase.verifyEqual(string(fieldnames(A)), string(fieldnames(B)));

   names = string(fieldnames(A));
   for i = 1:numel(names)
      name = char(names(i));
      icemodel.test.verify.verifyEqualNested(testCase, ...
         A.(name), B.(name), tol);
   end
end

function verifyArrayEqual(testCase, A, B, tol)

   testCase.verifyEqual(size(A), size(B));

   if isdatetime(A)
      testCase.verifyEqual(A, B);
      return
   end

   if isnumeric(A) || islogical(A)
      nan_mask = isnan(A) & isnan(B);
      finite_mask = isfinite(A) == isfinite(B);
      testCase.verifyTrue(all(nan_mask(:) | finite_mask(:)), ...
         'finite/NaN pattern mismatch');
      A = double(A(~nan_mask));
      B = double(B(~nan_mask));
      if isempty(A)
         return
      end
      testCase.verifyLessThanOrEqual(max(abs(A - B)), tol);
   else
      testCase.verifyEqual(A, B);
   end
end
