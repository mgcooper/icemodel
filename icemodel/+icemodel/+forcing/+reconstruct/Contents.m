% +RECONSTRUCT
%
%   Contents file for +RECONSTRUCT and its subfolders.
%
%   +RECONSTRUCT
%   icemodel.forcing.reconstruct.acceptanceWindow                 - Per-site forcing-ready policy window from staged proxies
%   icemodel.forcing.reconstruct.admissionGate                    - Apply the approved per-variable admission thresholds
%   icemodel.forcing.reconstruct.applyDonorTransfer               - Apply a fitted donor transfer to donor samples
%   icemodel.forcing.reconstruct.applyProxyCalibration            - Apply a fitted proxy calibration to model samples
%   icemodel.forcing.reconstruct.assertNotEvaluationDestination   - Refuse reconstruction writes under eval data
%   icemodel.forcing.reconstruct.assertPromiceFilledArtifact      - Prove product and station provenance
%   icemodel.forcing.reconstruct.auditSegments                    - Build one audit row per contiguous selected segment
%   icemodel.forcing.reconstruct.blendFallbackSeams               - Apply the policy seam taper to fallback fills
%   icemodel.forcing.reconstruct.blendSeams                       - Taper excess anchored boundary mismatch across one run
%   icemodel.forcing.reconstruct.bucketEdges                      - Return the gap-duration bucket edges in hours
%   icemodel.forcing.reconstruct.clearSkyIndex                    - Normalize shortwave flux by station-specific TOA irradiance
%   icemodel.forcing.reconstruct.climatologyFill                  - Day-of-year climatology estimate of one channel
%   icemodel.forcing.reconstruct.commonSupportSkill               - Compare candidate and baseline on identical samples
%   icemodel.forcing.reconstruct.deriveUpwardShortwave            - Fill missing swu from final albedo and swd
%   icemodel.forcing.reconstruct.elevationAdjust                  - Adjust a donor channel across an elevation difference
%   icemodel.forcing.reconstruct.fillPromiceStation               - Produce the gap-filled met product for one station
%   icemodel.forcing.reconstruct.fillShortGaps                    - Tier-1 bounded interior interpolation of one channel
%   icemodel.forcing.reconstruct.fillTwilightClimatology          - Fill one-posting SWD gaps beside known night
%   icemodel.forcing.reconstruct.fitDonorTransfer                 - Fit an overlap-calibrated donor-to-target transfer
%   icemodel.forcing.reconstruct.fitProxyCalibration              - Calibrate a model proxy on its observed overlap
%   icemodel.forcing.reconstruct.flatRunScreen                    - Flag multi-day buried/rime-encased sensor runs in met data
%   icemodel.forcing.reconstruct.gapCensus                        - Census contiguous missing runs in one role-contract series
%   icemodel.forcing.reconstruct.gapDurationBucket                - Assign positive durations to right-closed policy bins
%   icemodel.forcing.reconstruct.icemodelRequiredChannels         - The POLICY A5 seven-channel icemodel set
%   icemodel.forcing.reconstruct.interpolationCapHours            - Approved per-channel interpolation ceilings
%   icemodel.forcing.reconstruct.lastResortProxies                - Adopt aligned proxy values for residual gaps
%   icemodel.forcing.reconstruct.loadWidestTimetable              - Load the staged timetable with the widest time axis
%   icemodel.forcing.reconstruct.lwdEstimator                     - Empirical downward-longwave candidate from temperature and RH
%   icemodel.forcing.reconstruct.mustBeCapHours                   - Require a gap cap within any approved channel ceiling
%   icemodel.forcing.reconstruct.mustBeStationToken               - Require canonical lowercase alphanumeric station IDs
%   icemodel.forcing.reconstruct.partitionPrecipitation           - Split total precipitation by air temperature
%   icemodel.forcing.reconstruct.persistenceEstimate              - Hold pre-gap values without held-out-data leakage
%   icemodel.forcing.reconstruct.physicalBounds                   - Return the approved physical bounds for one channel
%   icemodel.forcing.reconstruct.physicalValidity                 - Enforce scalar and relational reconstruction bounds
%   POLICY.md
%   icemodel.forcing.reconstruct.policySha256                     - Return the SHA-256 fingerprint of reconstruction POLICY.md
%   icemodel.forcing.reconstruct.promiceFilledVerificationMatches - Match a prevalidated runtime identity
%   icemodel.forcing.reconstruct.provenanceCodes                  - Return the per-sample reconstruction provenance registry
%   icemodel.forcing.reconstruct.proxyArtifactIdentity            - Verify one staged proxy's target and producer
%   README.md
%   icemodel.forcing.reconstruct.reconstructSeries                - Compose admitted fill methods into one target series
%   icemodel.forcing.reconstruct.scalarValidity                   - Check finite samples against the A15 scalar registry
%   icemodel.forcing.reconstruct.seasonOf                         - Meteorological season label of each timestamp
%   icemodel.forcing.reconstruct.selectedDataRoot                 - Resolve one selected met path to its data and met roots
%   icemodel.forcing.reconstruct.setopts                          - Central options for the reconstruction pipeline
%   icemodel.forcing.reconstruct.smoothShortwaveSeams             - Repair empirical outlier boundaries in filled SWD
%   icemodel.forcing.reconstruct.solarElevationBands              - Solar-elevation thresholds for the swd science path
%   icemodel.forcing.reconstruct.stampGapfillIdentity             - Stamp the gapfill_* identity fields on a filled met
%   icemodel.forcing.reconstruct.stationMethodPlan                - Select and fit admitted fill methods for one station
%   icemodel.forcing.reconstruct.stepScale                        - Per-season median absolute step of one observed channel
%   icemodel.forcing.reconstruct.syntheticMissingness             - Draw blocked synthetic gaps into observed segments
%   icemodel.forcing.reconstruct.toaIrradiance                    - Top-of-atmosphere irradiance on a horizontal surface
%   icemodel.forcing.reconstruct.validationMetrics                - Grade reconstructed samples against withheld truth
%   icemodel.forcing.reconstruct.validationSplit                  - Partition station years into selection and evaluation sets
%   icemodel.forcing.reconstruct.verifyPromiceFilledReadiness     - Gate derived PROMICE forcing by coverage
%
%   updatecontents.m generated this file on 18 Sep 2026 at 01:12:19.
