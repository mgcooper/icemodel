% +FORCING
%
%   Contents file for +FORCING and its subfolders.
%
%   +FORCING
%   icemodel.forcing.buildGcnetVandecruxData                       - Build GC-Net/Vandecrux surface forcing Data
%   icemodel.forcing.buildGcnetVandecruxFirnTemperature            - Read GC-Net firn temperatures
%   icemodel.forcing.buildGcnetVandecruxMet                        - Build native met from GC-Net/Vandecrux surface data
%   icemodel.forcing.buildImauHourlyData                           - Build native IMAU hourly AWS Data
%   icemodel.forcing.buildImauHourlyMet                            - Build native met from IMAU hourly AWS data
%   icemodel.forcing.buildKtransectData                            - Build native K-transect 30-minute AWS Data
%   icemodel.forcing.buildKtransectMet                             - Build native met from K-transect 30-minute AWS data
%   icemodel.forcing.buildMarData                                  - Build a Data timetable from MAR v3.11 yearly NetCDF files
%   icemodel.forcing.buildMarMet                                   - Build an icemodel met timetable from MAR v3.11 data
%   icemodel.forcing.buildMerraData                                - Build a Data timetable from MERRA-2 daily NetCDF files
%   icemodel.forcing.buildMerraMet                                 - Build an icemodel met timetable from MERRA-2 data
%   icemodel.forcing.buildPromiceData                              - Build PROMICE data for model evaluation
%   icemodel.forcing.buildPromiceMet                               - Build an icemodel met timetable from PROMICE AWS data
%   icemodel.forcing.buildRacmoData                                - Build a Data timetable from RACMO2.3 NetCDF files
%   icemodel.forcing.buildSamimiDye2Data                           - Build RetMIP Dye-2 2016 native Data from Samimi AWS
%   icemodel.forcing.buildSamimiDye2Met                            - Build RetMIP Dye-2 2016 native met from Samimi AWS
%   icemodel.forcing.data2met                                      - Convert a Data timetable to an icemodel met timetable
%   icemodel.forcing.destepSurface                                 - Detect and optionally correct step-shifts in a surface
%   icemodel.forcing.fillPromiceAlbedo                             - Fill PROMICE albedo gaps with the winter-fill policy
%   icemodel.forcing.marGridInfo                                   - Read the MAR grid coordinates and static fields
%   icemodel.forcing.modisToMetCadence                             - Attach staged daily MODIS albedo onto a met time axis
%   icemodel.forcing.promiceAlbedoSourceValid                      - Identify finite native PROMICE albedo samples
%   icemodel.forcing.readGeusModis                                 - Read the GEUS MODIS daily albedo at points or a polygon
%   icemodel.forcing.readMar3p11                                   - Read one MAR v3.11 variable in standard units
%   README.md
%   icemodel.forcing.readMerra2                                    - Read one MERRA-2 variable in standard units
%   icemodel.forcing.readPromiceAws                                - Read a pypromice L3 AWS NetCDF into icemodel channels
%   icemodel.forcing.readRacmo2p3                                  - Read one RACMO 2.3 (FGRN11) variable in standard units
%   icemodel.forcing.stageModisAlbedo                              - Stage per-site GEUS MODIS daily albedo userdata artifacts
%
%   +FORCING/+HELPERS
%   icemodel.forcing.helpers.alignMarDailyMetadata                 - Align MAR per-day provenance to a retained time axis
%   icemodel.forcing.helpers.applyMarDailyQualityControl           - Constrain MAR hourly mass data by daily totals
%   icemodel.forcing.helpers.applyMarSnowDepthQualityControl       - Mask source-discontinuous SHSN2 years
%   icemodel.forcing.helpers.applyMerraTimeSupport                 - Apply the MERRA interval-start and support rules
%   icemodel.forcing.helpers.applyRacmoPrecipitationQualityControl - Enforce nonnegative RACMO ppt
%   icemodel.forcing.helpers.artifactCadenceMatches                - Check a saved timetable for the requested cadence
%   icemodel.forcing.helpers.artifactIdentityMatches               - Reject reuse across concrete provenance conflicts
%   icemodel.forcing.helpers.artifactMetadata                      - Build a source-light top-level artifact metadata record
%   icemodel.forcing.helpers.artifactScalarIdentityMatches         - Compare concrete scalar provenance facts
%   icemodel.forcing.helpers.attachLocationMetadata                - Add location CustomProperties to a Data timetable
%   icemodel.forcing.helpers.columnizeMetadata                     - Store metadata vectors as columns for inspection
%   icemodel.forcing.helpers.completeMetVariables                  - Add missing met-contract variables as NaN placeholders
%   icemodel.forcing.helpers.dailyAlbedoAnomalyFlags               - Flag transient reflected-shortwave collapses
%   icemodel.forcing.helpers.dailyToHourly                         - Interpolate daily data onto an hourly (or finer) time axis
%   icemodel.forcing.helpers.data2metCollection                    - Convert one Data timetable or a cell collection to met
%   icemodel.forcing.helpers.findEnclosingWindowFile               - Name of a staged window file bracketing a query
%   icemodel.forcing.helpers.gcnetHourlyAxis                       - Use the documented hourly row-index time convention
%   icemodel.forcing.helpers.gcnetTime                             - Convert Vandecrux/GC-Net numeric time to UTC datetimes
%   icemodel.forcing.helpers.gcnetVandecruxCatalog                 - The nine Vandecrux/GC-Net stations, defined once
%   icemodel.forcing.helpers.gcnetVandecruxInputs                  - Resolve shared Vandecrux/GC-Net roots and aliases
%   icemodel.forcing.helpers.gcnetVandecruxStation                 - Normalize Vandecrux/GC-Net station aliases
%   icemodel.forcing.helpers.gcnetVandecruxStationMetadata         - Return station aliases and location
%   icemodel.forcing.helpers.geusModisCoverageMetadata             - Build GEUS MODIS coverage provenance
%   icemodel.forcing.helpers.geusModisProjection                   - Native polar-stereographic projection of the GEUS grid
%   icemodel.forcing.helpers.gridLocation                          - Map a point or polygon onto a grid hyperslab + collapse rule
%   icemodel.forcing.helpers.hasCanonicalMerraTimeSupport          - True for the complete MERRA time contract
%   icemodel.forcing.helpers.hasConstantMerraTavg3Support          - True when glacier channels hold each UTC block
%   icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid         - True for an exact native tavg3 inventory
%   icemodel.forcing.helpers.interpRcm                             - Interpolate RCM grid-cell data to query points
%   icemodel.forcing.helpers.intervalMaximumSolarElevation         - Maximum sun angle over each time support
%   icemodel.forcing.helpers.listSourceFiles                       - Return flat and nested files under a source/cache root
%   icemodel.forcing.helpers.locateImauHourlyFile                  - Find a flat or hourly-subfolder PANGAEA tab file
%   icemodel.forcing.helpers.locateKtransectFiles                  - Find every cached K-transect annual file for a station
%   icemodel.forcing.helpers.locateNormalizedFile                  - Find a file by normalized path/name token
%   icemodel.forcing.helpers.locateSamimiDye2Workbook              - Find the Dye-2 summer 2016 AWS workbook
%   icemodel.forcing.helpers.marDiagnosticMetadata                 - Describe optional MAR mass-balance diagnostics
%   icemodel.forcing.helpers.marDynamicProfileQa                   - Diagnose MAR dynamic snow/firn layer consistency
%   icemodel.forcing.helpers.marRefreezeMetadata                   - Describe signed native MAR RZ values in one artifact
%   icemodel.forcing.helpers.metchecks                             - Gap-fill and clamp met variables to physically valid ranges
%   icemodel.forcing.helpers.metfilename                           - Build a standard icemodel met-file name
%   icemodel.forcing.helpers.metTimestepSuffix                     - Return the file tag for a model-met cadence
%   icemodel.forcing.helpers.metvariables                          - Met-file variable names for the forcing builders
%   icemodel.forcing.helpers.modisAlbedoChannel                    - GEUS MODIS daily albedo on a time axis, per location
%   icemodel.forcing.helpers.normalizedFileToken                   - Compare filenames case-insensitively across separators
%   icemodel.forcing.helpers.normalizeGeusModisAlbedo              - Mask undocumented GEUS albedo sentinels
%   icemodel.forcing.helpers.normalizeLocations                    - Normalize one location or a point list to a row cell
%   icemodel.forcing.helpers.optionalDate                          - Normalize optional date-like inputs to UTC datetimes
%   icemodel.forcing.helpers.precipitationConsistency              - Check finite, nonnegative phase mass balance
%   icemodel.forcing.helpers.precipitationValidity                 - Validate complete and partial precipitation splits
%   icemodel.forcing.helpers.precipitationVariables                - Total and partitioned precipitation names
%   icemodel.forcing.helpers.projectLocation                       - Ensure a location struct carries EPSG:3413 coordinates
%   icemodel.forcing.helpers.promiceShortwave                      - Select physical public PROMICE shortwave channels
%   icemodel.forcing.helpers.pruneSupersededWindowFiles            - Remove shorter windows contained by a new file
%   icemodel.forcing.helpers.psnProjection                         - Polar stereographic north projection used by the builders
%   icemodel.forcing.helpers.readGcnetDonor                        - Load one GC-Net surface file as an observed-only donor
%   icemodel.forcing.helpers.readImauHourlyTable                   - Parse one IMAU hourly PANGAEA AWS table
%   icemodel.forcing.helpers.readKtransectHeights                  - Parse the K-transect sensor-height workbook
%   icemodel.forcing.helpers.readKtransectTable                    - Parse one K-transect annual PANGAEA AWS table
%   icemodel.forcing.helpers.readMarDensitySnapshots               - Read requested MAR RO1 density profiles
%   icemodel.forcing.helpers.readMerra2Time                        - Decode only a MERRA-2 file's native UTC time coordinate
%   icemodel.forcing.helpers.readNetcdfAttribute                   - Return one NetCDF attribute or empty string when absent
%   icemodel.forcing.helpers.readPangaeaTab                        - Read one PANGAEA tab-delimited dataset export
%   icemodel.forcing.helpers.regexpOnce                            - Return one stripped regexp token or an empty string
%   icemodel.forcing.helpers.remapPolygon                          - Conservative area-weighted remap of a grid block to a polygon
%   icemodel.forcing.helpers.resampleMetTimestep                   - Resample model met without crossing source outages
%   icemodel.forcing.helpers.slabMean                              - Mean of a static grid field over hyperslab target cells
%   icemodel.forcing.helpers.solarElevation                        - Approximate NOAA geometric solar elevation in degrees
%   icemodel.forcing.helpers.sourceAlbedo                          - Derive broadband albedo where shortwave input is valid
%   icemodel.forcing.helpers.sourceSearchDirs                      - Candidate dirs for a staged file: per-source subfolder first
%   icemodel.forcing.helpers.sourceVariableAttributes              - Preserve NetCDF variable attributes by variable
%   icemodel.forcing.helpers.stampMetadata                         - Embed CF-ish metadata in a timetable's properties
%   icemodel.forcing.helpers.stationTransitionTimes                - Within-record AWS handover times for a PROMICE site
%   icemodel.forcing.helpers.surfaceFlags                          - Per-sample quality flags for a PROMICE surface-height series
%   icemodel.forcing.helpers.timeWindowMask                        - Select an optional datetime window from a source time axis
%   icemodel.forcing.helpers.uniformCadenceSeconds                 - Derive one timetable's exact regular cadence
%   icemodel.forcing.helpers.validatemet                           - Assert that MET satisfies the icemodel met-file contract
%   icemodel.forcing.helpers.variableUnits                         - Unit string for each forcing-builder channel
%   icemodel.forcing.helpers.verificationSourceDir                 - Resolve repo-local verification source data roots
%   icemodel.forcing.helpers.windFromComponents                    - Wind speed and direction from zonal/meridional components
%   icemodel.forcing.helpers.writemet                              - Validate and save an icemodel met file
%   icemodel.forcing.helpers.writeuserdata                         - Save a Data timetable as met-swap userdata files
%
%   +FORCING/+RECONSTRUCT
%   icemodel.forcing.reconstruct.acceptanceWindow                  - Per-site forcing-ready policy window from staged proxies
%   icemodel.forcing.reconstruct.admissionGate                     - Apply the approved per-variable admission thresholds
%   icemodel.forcing.reconstruct.applyDonorTransfer                - Apply a fitted donor transfer to donor samples
%   icemodel.forcing.reconstruct.applyProxyCalibration             - Apply a fitted proxy calibration to model samples
%   icemodel.forcing.reconstruct.assertNotEvaluationDestination    - Refuse reconstruction writes under eval data
%   icemodel.forcing.reconstruct.assertPromiceFilledArtifact       - Prove product and station provenance
%   icemodel.forcing.reconstruct.auditSegments                     - Build one audit row per contiguous selected segment
%   icemodel.forcing.reconstruct.blendFallbackSeams                - Apply the policy seam taper to fallback fills
%   icemodel.forcing.reconstruct.blendSeams                        - Taper excess anchored boundary mismatch across one run
%   icemodel.forcing.reconstruct.bucketEdges                       - Return the gap-duration bucket edges in hours
%   icemodel.forcing.reconstruct.clearSkyIndex                     - Normalize shortwave flux by station-specific TOA irradiance
%   icemodel.forcing.reconstruct.climatologyFill                   - Day-of-year climatology estimate of one channel
%   icemodel.forcing.reconstruct.commonSupportSkill                - Compare candidate and baseline on identical samples
%   icemodel.forcing.reconstruct.deriveUpwardShortwave             - Fill missing swu from final albedo and swd
%   icemodel.forcing.reconstruct.elevationAdjust                   - Adjust a donor channel across an elevation difference
%   icemodel.forcing.reconstruct.fillPromiceStation                - Produce the gap-filled met product for one station
%   icemodel.forcing.reconstruct.fillShortGaps                     - Tier-1 bounded interior interpolation of one channel
%   icemodel.forcing.reconstruct.fillTwilightClimatology           - Fill one-posting SWD gaps beside known night
%   icemodel.forcing.reconstruct.fitDonorTransfer                  - Fit an overlap-calibrated donor-to-target transfer
%   icemodel.forcing.reconstruct.fitProxyCalibration               - Calibrate a model proxy on its observed overlap
%   icemodel.forcing.reconstruct.flatRunScreen                     - Flag multi-day buried/rime-encased sensor runs in met data
%   icemodel.forcing.reconstruct.gapCensus                         - Census contiguous missing runs in one role-contract series
%   icemodel.forcing.reconstruct.gapDurationBucket                 - Assign positive durations to right-closed policy bins
%   icemodel.forcing.reconstruct.icemodelRequiredChannels          - The POLICY A5 seven-channel icemodel set
%   icemodel.forcing.reconstruct.interpolationCapHours             - Approved per-channel interpolation ceilings
%   icemodel.forcing.reconstruct.lastResortProxies                 - Adopt aligned proxy values for residual gaps
%   icemodel.forcing.reconstruct.loadWidestTimetable               - Load the staged timetable with the widest time axis
%   icemodel.forcing.reconstruct.lwdEstimator                      - Empirical downward-longwave candidate from temperature and RH
%   icemodel.forcing.reconstruct.mustBeCapHours                    - Require a gap cap within any approved channel ceiling
%   icemodel.forcing.reconstruct.mustBeStationToken                - Require canonical lowercase alphanumeric station IDs
%   icemodel.forcing.reconstruct.partitionPrecipitation            - Split total precipitation by air temperature
%   icemodel.forcing.reconstruct.persistenceEstimate               - Hold pre-gap values without held-out-data leakage
%   icemodel.forcing.reconstruct.physicalBounds                    - Return the approved physical bounds for one channel
%   icemodel.forcing.reconstruct.physicalValidity                  - Enforce scalar and relational reconstruction bounds
%   POLICY.md
%   icemodel.forcing.reconstruct.policySha256                      - Return the SHA-256 fingerprint of reconstruction POLICY.md
%   icemodel.forcing.reconstruct.promiceFilledVerificationMatches  - Match a prevalidated runtime identity
%   icemodel.forcing.reconstruct.provenanceCodes                   - Return the per-sample reconstruction provenance registry
%   icemodel.forcing.reconstruct.proxyArtifactIdentity             - Verify one staged proxy's target and producer
%   README.md
%   icemodel.forcing.reconstruct.reconstructSeries                 - Compose admitted fill methods into one target series
%   icemodel.forcing.reconstruct.scalarValidity                    - Check finite samples against the A15 scalar registry
%   icemodel.forcing.reconstruct.seasonOf                          - Meteorological season label of each timestamp
%   icemodel.forcing.reconstruct.selectedDataRoot                  - Resolve one selected met path to its data and met roots
%   icemodel.forcing.reconstruct.setopts                           - Central options for the reconstruction pipeline
%   icemodel.forcing.reconstruct.smoothShortwaveSeams              - Repair empirical outlier boundaries in filled SWD
%   icemodel.forcing.reconstruct.solarElevationBands               - Solar-elevation thresholds for the swd science path
%   icemodel.forcing.reconstruct.stampGapfillIdentity              - Stamp the gapfill_* identity fields on a filled met
%   icemodel.forcing.reconstruct.stationMethodPlan                 - Select and fit admitted fill methods for one station
%   icemodel.forcing.reconstruct.stepScale                         - Per-season median absolute step of one observed channel
%   icemodel.forcing.reconstruct.syntheticMissingness              - Draw blocked synthetic gaps into observed segments
%   icemodel.forcing.reconstruct.toaIrradiance                     - Top-of-atmosphere irradiance on a horizontal surface
%   icemodel.forcing.reconstruct.validationMetrics                 - Grade reconstructed samples against withheld truth
%   icemodel.forcing.reconstruct.validationSplit                   - Partition station years into selection and evaluation sets
%   icemodel.forcing.reconstruct.verifyPromiceFilledReadiness      - Gate derived PROMICE forcing by coverage
%
%   updatecontents.m generated this file on 18 Sep 2026 at 01:12:18.
