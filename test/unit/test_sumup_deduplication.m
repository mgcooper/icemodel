function tests = test_sumup_deduplication
   %TEST_SUMUP_DEDUPLICATION Verify exact SUMup scientific row identities.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % Install the package path without creating any staged artifacts.
   [~, ~, ~, ~, cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
end

function teardownOnce(testCase)
   % Release the shared test environment after all direct helper procedures.
   testCase.TestData.cleanup = [];
end

function test_density_identity_missing_zero_signed_zero_and_provenance(testCase)
   % Every density identity field must distinguish a row, while repeated missing
   % uncertainty and signed zero compare equal and retain first-row provenance.
   [record, identity_names] = identityFixture("density");

   [actual, raw_rows, unique_rows, removed_rows] = ...
      icemodel.verification.setup.deduplicateSumupRecords(record, "density");

   verifyIdentityFixture(testCase, actual, identity_names, ...
      raw_rows, unique_rows, removed_rows);
   testCase.verifyEqual(actual.name_key(1), 101);
   testCase.verifyEqual(actual.name(1), "first_alias");
   testCase.verifyEqual(actual.elevation(1), 1001);
   testCase.verifyEqual(actual.measurement_id(1), 1);
   testCase.verifyEqual(actual.measurement_id, ...
      [1; (3:height(record))']);
   % The base identity and every non-error perturbation legitimately retain
   % missing uncertainty; only the error perturbation carries genuine zero.
   testCase.verifyEqual(nnz(ismissing(actual.error)), numel(identity_names));
   testCase.verifyEqual(nnz(actual.error == 0), 1);
end

function test_temperature_identity_fields(testCase)
   % Temperature depth/support/time/value fields all participate in identity.
   [record, identity_names] = identityFixture("temperature");

   [actual, raw_rows, unique_rows, removed_rows] = ...
      icemodel.verification.setup.deduplicateSumupRecords( ...
      record, "temperature");

   verifyIdentityFixture(testCase, actual, identity_names, ...
      raw_rows, unique_rows, removed_rows);
end

function test_smb_identity_fields(testCase)
   % SMB period endpoints/years/value/provenance keys all participate in identity.
   [record, identity_names] = identityFixture("SMB");

   [actual, raw_rows, unique_rows, removed_rows] = ...
      icemodel.verification.setup.deduplicateSumupRecords(record, "SMB");

   verifyIdentityFixture(testCase, actual, identity_names, ...
      raw_rows, unique_rows, removed_rows);
end

function test_empty_input_and_invalid_contracts(testCase)
   % A schema-valid empty selection reports zero counts; unknown variables and
   % missing identity columns fail before silently weakening equality.
   [record, ~] = identityFixture("density");
   empty_record = record([], :);

   [actual, raw_rows, unique_rows, removed_rows] = ...
      icemodel.verification.setup.deduplicateSumupRecords( ...
      empty_record, "density");
   testCase.verifyEqual(height(actual), 0);
   testCase.verifyEqual([raw_rows, unique_rows, removed_rows], [0, 0, 0]);

   testCase.verifyError(@() ...
      icemodel.verification.setup.deduplicateSumupRecords(record, "unknown"), ...
      'icemodel:verification:deduplicateSumupRecords:badVariable');
   missing_column = removevars(record, 'longitude');
   testCase.verifyError(@() ...
      icemodel.verification.setup.deduplicateSumupRecords( ...
      missing_column, "density"), ...
      'icemodel:verification:deduplicateSumupRecords:missingIdentity');
end

function verifyIdentityFixture(testCase, actual, identity_names, ...
      raw_rows, unique_rows, removed_rows)
   %VERIFYIDENTITYFIXTURE Check the common all-fields-and-one-duplicate contract.
   expected_raw = numel(identity_names) + 2;
   expected_unique = numel(identity_names) + 1;
   testCase.verifyEqual(raw_rows, expected_raw);
   testCase.verifyEqual(unique_rows, expected_unique);
   testCase.verifyEqual(removed_rows, 1);
   testCase.verifyEqual(height(actual), expected_unique);
end

function [record, identity_names] = identityFixture(variable)
   %IDENTITYFIXTURE Create one duplicate plus one perturbation per identity field.
   switch lower(variable)
      case "density"
         identity_names = ["latitude", "longitude", "start_depth", ...
            "stop_depth", "midpoint", "timestamp", "reference_key", ...
            "method_key", "density", "error"];
      case "temperature"
         identity_names = ["latitude", "longitude", "depth", "duration", ...
            "timestamp", "reference_key", "method_key", "temperature", ...
            "error"];
      case "smb"
         identity_names = ["latitude", "longitude", "start_date", ...
            "end_date", "start_year", "end_year", "reference_key", ...
            "method_key", "smb", "error"];
   end

   % Rows one and two are scientifically identical despite repeated missing
   % uncertainty and, for density, positive versus negative zero. Each remaining
   % row changes exactly one identity field and therefore must survive.
   n_identity = numel(identity_names);
   base = 10 .* (1:n_identity);
   error_index = find(identity_names == "error", 1);
   base(error_index) = NaN;
   value_index = find(identity_names == lower(variable), 1);
   if lower(variable) == "density"
      base(value_index) = 0;
   end
   values = repmat(base, n_identity + 2, 1);
   if lower(variable) == "density"
      values(2, value_index) = -0;
   end
   for k = 1:n_identity
      if ismissing(base(k))
         values(k + 2, k) = 0;
      else
         values(k + 2, k) = base(k) + 1;
      end
   end
   record = array2table(values, 'VariableNames', cellstr(identity_names));

   % These source fields are intentionally outside identity. Distinct duplicate-
   % row values prove the untouched first row supplies retained provenance.
   n_rows = height(record);
   record.elevation = (1001:1000 + n_rows)';
   record.name_key = (101:100 + n_rows)';
   record.name = "alias_" + string((1:n_rows)');
   record.name(1:2) = ["first_alias"; "second_alias"];
   record.measurement_id = (1:n_rows)';
end
