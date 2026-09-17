% +HELPERS
%
%   Contents file for +HELPERS and its subfolders.
%
%   +HELPERS
%   icemodel.verification.helpers.ablationLedgerIncrements         - Per-interval ablation terms from the mass ledger
%   icemodel.verification.helpers.alignObservationSeries           - Align one observation/model series by its support
%   icemodel.verification.helpers.assertArtifactSha256             - Require current bytes to match a pinned SHA-256
%   icemodel.verification.helpers.assertRootRelativeArtifactSha256 - Verify one root-scoped artifact identity
%   icemodel.verification.helpers.classifyObservationSupport       - Apply the flag-support rules to observation rows
%   icemodel.verification.helpers.classifySnowDepth                - Classify snow support without accepting negative depth
%   icemodel.verification.helpers.esmRuntimeMetFiles               - Resolve an atomic ESM case's standard runtime met paths
%   icemodel.verification.helpers.esmSnowmipWaterYear              - Return one snow water year for an ESM-SnowMIP site
%   icemodel.verification.helpers.evaluationDataRoot               - Resolve the base evaluation-data root
%   icemodel.verification.helpers.evaluationSeason                 - Summertime display and evaluation bounds for one year
%   icemodel.verification.helpers.familyManifestFiles              - List staged verification family manifest files
%   icemodel.verification.helpers.fieldOr                          - Return a struct field or default value
%   icemodel.verification.helpers.inputDataRoot                    - Resolve the base icemodel input-data root
%   icemodel.verification.helpers.isPhysicsFingerprint             - Return true for one well-formed physics stamp
%   icemodel.verification.helpers.loadArtifact                     - Load one named staged verification artifact from a MAT file
%   icemodel.verification.helpers.loadColocatedData                - Assemble a timeseries bundle from staged per-source files
%   icemodel.verification.helpers.metricRowSchema                  - Return canonical comparison-metric field names/defaults
%   icemodel.verification.helpers.observationSupportFields         - Columns the observation support rules read
%   icemodel.verification.helpers.physicsFingerprint               - Fingerprint the model's default physics configuration
%   icemodel.verification.helpers.profileGroups                    - Split profile rows by stable source identity and UTC date
%   icemodel.verification.helpers.readFamilyManifest               - Read one verification family manifest JSON file
%   icemodel.verification.helpers.residualMetrics                  - Bias, MAE, RMSE, max error, and NSE for one paired series
%   icemodel.verification.helpers.resolveCandidateBundle           - Resolve the comparison bundle for one case
%   icemodel.verification.helpers.sampleQuantile                   - Linear-interpolated sample quantile, no extra toolboxes
%   icemodel.verification.helpers.sumupColocation                  - Flag a SUMup point as co-located with a mixed anchor
%   icemodel.verification.helpers.validateAblationModelSchema      - Check a saved cohort still supports the report
%   icemodel.verification.helpers.writeRunReport                   - Write a concise markdown report for a verification run
%
%   updatecontents.m generated this file on 14 Sep 2026 at 20:14:01.
