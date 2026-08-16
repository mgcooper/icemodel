function tf = isPhysicsFingerprint(value)
   %ISPHYSICSFINGERPRINT Return true for one well-formed physics stamp.
   %
   %  tf = icemodel.verification.helpers.isPhysicsFingerprint(value)
   %
   % A stamp is one scalar struct carrying opts_sha256 and icemodel_version,
   % each a nonempty text scalar. VALUE comes out of a saved results MAT file
   % and can hold anything. Every test below therefore short-circuits before
   % the one after it reads a field.
   %
   % Shape only. Neither field has a fixed format here: requiring 64 hex
   % characters would reject every version string.
   %
   % See also: icemodel.verification.helpers.physicsFingerprint,
   %  icemodel.verification.report.buildAblationEvaluationReport

   fields = ["opts_sha256", "icemodel_version"];
   tf = isstruct(value) && isscalar(value) && all(isfield(value, fields)) ...
      && isTextScalar(value.opts_sha256) ...
      && isTextScalar(value.icemodel_version);
end

function tf = isTextScalar(value)
   %ISTEXTSCALAR Return true for one nonempty text scalar.
   %
   % The type test comes first. string() therefore never runs on a numeric or
   % cell value, and cannot turn one into text that passes the later tests.

   tf = (isstring(value) || ischar(value)) && isscalar(string(value)) ...
      && ~ismissing(string(value)) && strlength(string(value)) > 0;
end
