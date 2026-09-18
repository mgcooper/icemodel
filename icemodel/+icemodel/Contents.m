% +ICEMODEL
%
%   Contents file for +ICEMODEL and its subfolders.
%
%   +ICEMODEL
%   icemodel.buildOutputPayload                                             - Return icemodel.updateoutputs cell arrays in OPTS order
%   icemodel.chunkgridcell                                                  - This function is incomplete. rungrid_1 would call it like this:
%   icemodel.completions                                                    - Function completions
%   icemodel.concatoutput                                                   - Concatenate yearly icemodel output structures
%   icemodel.config                                                         - Configure icemodel project paths
%   icemodel.configureRun                                                   - Fill derived run settings from OPTS
%   icemodel.createMetFileNames                                             - Create icemodel met file names for model OPTS
%   icemodel.cvconvert                                                      - Convert between control volume properties
%   icemodel.dependencies                                                   - Add external project dependencies to the MATLAB path
%   icemodel.extractice2                                                    - Extract 2-d icemodel data
%   icemodel.getopts                                                        - Return model options by name
%   icemodel.getpath                                                        - Return canonical icemodel data and run paths
%   icemodel.interpmet                                                      - Interpolate met data
%   icemodel.isIncrementChannel                                             - True for per-step increment channels
%   icemodel.isPathInside                                                   - True when a canonical path resolves within a selected root
%   icemodel.loadmet                                                        - Load one or more icemodel met files as a timetable
%   icemodel.loadRestartState                                               - Load a saved year-boundary restart state
%   icemodel.loadresults                                                    - Load saved yearly output files and concatenate them if needed
%   icemodel.mkfolders
%   icemodel.outputYears                                                    - Return the simulation years retained in saved output
%   icemodel.pairedWindow                                                   - Normalize one optional start/end pair to UTC datetimes
%   icemodel.parameterLookup                                                - Return the value of a model parameter
%   icemodel.physicalConstant                                               - Return the value of a physical constant
%   icemodel.postprocess                                                    - Calculate diagnostic outputs and format simulation output
%   icemodel.prepareRunOutput                                               - Prepare output folders and optional opts logging
%   icemodel.processmet                                                     - Post-process an icemodel met timetable
%   README.md
%   icemodel.resetopts                                                      - Override existing OPTS fields by name
%   icemodel.resolvePrecipPhase                                             - Select the runtime rainf/snowf split
%   icemodel.restartfile                                                    - Return the canonical restart-state file path for a run/year
%   icemodel.retimeHourlyFixedStep                                          - Aggregate fixed-step data to hourly values
%   icemodel.saveRestartState                                               - Save the year-boundary state needed for a restart
%   icemodel.saveRunOpts                                                    - Save the resolved OPTS struct for a run
%   icemodel.setcase                                                        - Input parsing
%   icemodel.setopts                                                        - Set model options
%   icemodel.shellQuote                                                     - Quote one value as a literal argument for the host shell
%   icemodel.time2iter
%   icemodel.updateoutput                                                   - Store one timestep of model output into the ice1/ice2 structs
%   icemodel.writeoutput                                                    - Post-process and save ice1/ice2 output to disk
%
%   +ICEMODEL/+COLUMN
%   icemodel.column.accumulate_phase_budget                                 - Add one substep's phase-change storage increments
%   icemodel.column.accumulate_remesh_budget                                - Add one substep's remesh events to the budget
%   icemodel.column.accumulate_vapor_exchange                               - Budget one substep's surface vapor exchange
%   icemodel.column.accumulate_vapor_transport                              - Budget one substep's interior vapor transport
%   icemodel.column.apply_vapor_transfer                                    - Apply vapor phase change increments to the column
%   icemodel.column.assemble_enthalpy_system                                - Compute the general equation coefficients
%   icemodel.column.assert_max_water                                        - Assert that water fraction does not exceed the maximum
%   icemodel.column.available_liquid_water                                  - Compute available liquid water in one control volume
%   icemodel.column.budget_surface_mass_balance                             - Apply and budget the surface mass balance
%   icemodel.column.bulk_density                                            - Bulk density of ice, liquid, and air mixture
%   icemodel.column.bulk_enthalpy                                           - Compute the solver-state bulk enthalpy [J m-3]
%   icemodel.column.bulk_specific_heat_capacity                             - Bulk cp of ice, liquid, and air mixture
%   icemodel.column.bulk_thermal_conductivity                               - Compute bulk effective thermal conductivity
%   icemodel.column.control_volume_mesh                                     - Compute cell edges and nodes for the column mesh
%   icemodel.column.couple_vapor_step                                       - Apply subsurface vapor transport across cells
%   icemodel.column.diagnose_column_runoff                                  - Diagnose cumulative runoff from column mass changes
%   icemodel.column.enforce_control_volume_balance                          - Enforce the total-volume constraint
%   icemodel.column.finalize_budget_state                                   - Record the storage end endpoints for one forcing step
%   icemodel.column.firn_thermal_conductivity                               - Compute porous ice thermal conductivity
%   icemodel.column.infiltration                                            - Snow column liquid mass + cold-content + conduction update
%   icemodel.column.initialize_budget_state                                 - Return the zeroed budget for one forcing step
%   icemodel.column.initialize_column_state                                 - Initialize the 1-d ice column state
%   icemodel.column.initialize_remesh_ledger                                - Zeroed remesh event ledger
%   icemodel.column.integrate_column_budget                                 - Return column-integrated mass and enthalpy storage
%   icemodel.column.liquid_flux                                             - Compute the liquid water flux between snowpack layers
%   icemodel.column.liquid_fraction_derivative                              - Liquid fraction derivative wrt temperature
%   icemodel.column.liquid_fraction_function                                - Project state onto the liquid-fraction function
%   icemodel.column.max_liquid_fraction_change                              - Largest f_liq increase a control volume takes
%   icemodel.column.meltzone_bounds                                         - Return the canonical mushy-zone liquid fraction bounds
%   icemodel.column.meltzone_transform                                      - Apply the melt-zone temperature-enthalpy transform
%   icemodel.column.merge_layer_indices                                     - Choose the pair of layers to merge
%   icemodel.column.merge_layers                                            - Combine two control volumes conserving state and sources
%   icemodel.column.merge_thin_layers                                       - Merge layers that fall below the minimum ice fraction
%   icemodel.column.potential_sublimation                                   - Convert surface vapor demand to ice fraction
%   README.md
%   icemodel.column.residual_water_fraction                                 - Volumetric residual-water floor for control volumes
%   icemodel.column.residual_water_pore_fraction                            - Residual liquid fraction per pore volume
%   icemodel.column.saturated_hydraulic_conductivity                        - Snow saturated hydraulic conductivity
%   icemodel.column.shortwave_source_term                                   - Solve the spectral shortwave source term
%   icemodel.column.solve_column_enthalpy                                   - Solve the column enthalpy balance
%   icemodel.column.solve_column_temperature                                - Solve the 1-dimensional column conduction equation
%   icemodel.column.subsurface_linearization_error                          - Diagnose the top-node enthalpy
%   icemodel.column.surface_linearization_error                             - Diagnose the Robin surface linearization error
%   icemodel.column.update_grain_radius                                     - Grow thermal grains from the substep vapor exchange
%   icemodel.column.updatestate                                             - Update column thermodynamic state variables
%   icemodel.column.vapor_exchange_is_wet                                   - Decide which phase a cell exchanges vapor with
%   icemodel.column.vapor_transport_terms                                   - Build the coupled vapor face transport terms
%   icemodel.column.water_fraction                                          - Compute the total volumetric water fraction
%
%   +ICEMODEL/+COUPLERS
%   icemodel.couplers.accelerate_coupler_iterate                            - Accelerate one surface-temperature Picard step
%   icemodel.couplers.initialize_coupler_history                            - Initialize the coupler-iteration history
%   icemodel.couplers.initialize_solver_diag                                - Initialize forcing-step solver diagnostics
%   icemodel.couplers.initialize_solver_settings                            - Initialize solver and timestep settings
%   README.md
%   icemodel.couplers.solve_skin_surface_column                             - Coupled skin-subsurface T_sfc-T_ice solve
%   icemodel.couplers.solve_surface_column_dirichlet                        - Coupled surface-subsurface solve for
%   icemodel.couplers.solve_surface_column_robin                            - Coupled surface-subsurface solve for Robin-type
%   icemodel.couplers.update_solver_diag                                    - Record the accepted substep solver result
%
%   +ICEMODEL/+FORCING
%   icemodel.forcing.buildGcnetVandecruxData                                - Build GC-Net/Vandecrux surface forcing Data
%   icemodel.forcing.buildGcnetVandecruxFirnTemperature                     - Read GC-Net firn temperatures
%   icemodel.forcing.buildGcnetVandecruxMet                                 - Build native met from GC-Net/Vandecrux surface data
%   icemodel.forcing.buildImauHourlyData                                    - Build native IMAU hourly AWS Data
%   icemodel.forcing.buildImauHourlyMet                                     - Build native met from IMAU hourly AWS data
%   icemodel.forcing.buildKtransectData                                     - Build native K-transect 30-minute AWS Data
%   icemodel.forcing.buildKtransectMet                                      - Build native met from K-transect 30-minute AWS data
%   icemodel.forcing.buildMarData                                           - Build a Data timetable from MAR v3.11 yearly NetCDF files
%   icemodel.forcing.buildMarMet                                            - Build an icemodel met timetable from MAR v3.11 data
%   icemodel.forcing.buildMerraData                                         - Build a Data timetable from MERRA-2 daily NetCDF files
%   icemodel.forcing.buildMerraMet                                          - Build an icemodel met timetable from MERRA-2 data
%   icemodel.forcing.buildPromiceData                                       - Build PROMICE data for model evaluation
%   icemodel.forcing.buildPromiceMet                                        - Build an icemodel met timetable from PROMICE AWS data
%   icemodel.forcing.buildRacmoData                                         - Build a Data timetable from RACMO2.3 NetCDF files
%   icemodel.forcing.buildSamimiDye2Data                                    - Build RetMIP Dye-2 2016 native Data from Samimi AWS
%   icemodel.forcing.buildSamimiDye2Met                                     - Build RetMIP Dye-2 2016 native met from Samimi AWS
%   icemodel.forcing.data2met                                               - Convert a Data timetable to an icemodel met timetable
%   icemodel.forcing.destepSurface                                          - Detect and optionally correct step-shifts in a surface
%   icemodel.forcing.fillPromiceAlbedo                                      - Fill PROMICE albedo gaps with the winter-fill policy
%   icemodel.forcing.marGridInfo                                            - Read the MAR grid coordinates and static fields
%   icemodel.forcing.modisToMetCadence                                      - Attach staged daily MODIS albedo onto a met time axis
%   icemodel.forcing.promiceAlbedoSourceValid                               - Identify finite native PROMICE albedo samples
%   icemodel.forcing.readGeusModis                                          - Read the GEUS MODIS daily albedo at points or a polygon
%   icemodel.forcing.readMar3p11                                            - Read one MAR v3.11 variable in standard units
%   README.md
%   icemodel.forcing.readMerra2                                             - Read one MERRA-2 variable in standard units
%   icemodel.forcing.readPromiceAws                                         - Read a pypromice L3 AWS NetCDF into icemodel channels
%   icemodel.forcing.readRacmo2p3                                           - Read one RACMO 2.3 (FGRN11) variable in standard units
%   icemodel.forcing.stageModisAlbedo                                       - Stage per-site GEUS MODIS daily albedo userdata artifacts
%
%   +ICEMODEL/+FORCING/+HELPERS
%   icemodel.forcing.helpers.alignMarDailyMetadata                          - Align MAR per-day provenance to a retained time axis
%   icemodel.forcing.helpers.applyMarDailyQualityControl                    - Constrain MAR hourly mass data by daily totals
%   icemodel.forcing.helpers.applyMarSnowDepthQualityControl                - Mask source-discontinuous SHSN2 years
%   icemodel.forcing.helpers.applyMerraTimeSupport                          - Apply the MERRA interval-start and support rules
%   icemodel.forcing.helpers.applyRacmoPrecipitationQualityControl          - Enforce nonnegative RACMO ppt
%   icemodel.forcing.helpers.artifactCadenceMatches                         - Check a saved timetable for the requested cadence
%   icemodel.forcing.helpers.artifactIdentityMatches                        - Reject reuse across concrete provenance conflicts
%   icemodel.forcing.helpers.artifactMetadata                               - Build a source-light top-level artifact metadata record
%   icemodel.forcing.helpers.artifactScalarIdentityMatches                  - Compare concrete scalar provenance facts
%   icemodel.forcing.helpers.attachLocationMetadata                         - Add location CustomProperties to a Data timetable
%   icemodel.forcing.helpers.columnizeMetadata                              - Store metadata vectors as columns for inspection
%   icemodel.forcing.helpers.completeMetVariables                           - Add missing met-contract variables as NaN placeholders
%   icemodel.forcing.helpers.dailyAlbedoAnomalyFlags                        - Flag transient reflected-shortwave collapses
%   icemodel.forcing.helpers.dailyToHourly                                  - Interpolate daily data onto an hourly (or finer) time axis
%   icemodel.forcing.helpers.data2metCollection                             - Convert one Data timetable or a cell collection to met
%   icemodel.forcing.helpers.findEnclosingWindowFile                        - Name of a staged window file bracketing a query
%   icemodel.forcing.helpers.gcnetHourlyAxis                                - Use the documented hourly row-index time convention
%   icemodel.forcing.helpers.gcnetTime                                      - Convert Vandecrux/GC-Net numeric time to UTC datetimes
%   icemodel.forcing.helpers.gcnetVandecruxCatalog                          - The nine Vandecrux/GC-Net stations, defined once
%   icemodel.forcing.helpers.gcnetVandecruxInputs                           - Resolve shared Vandecrux/GC-Net roots and aliases
%   icemodel.forcing.helpers.gcnetVandecruxStation                          - Normalize Vandecrux/GC-Net station aliases
%   icemodel.forcing.helpers.gcnetVandecruxStationMetadata                  - Return station aliases and location
%   icemodel.forcing.helpers.geusModisCoverageMetadata                      - Build GEUS MODIS coverage provenance
%   icemodel.forcing.helpers.geusModisProjection                            - Native polar-stereographic projection of the GEUS grid
%   icemodel.forcing.helpers.gridLocation                                   - Map a point or polygon onto a grid hyperslab + collapse rule
%   icemodel.forcing.helpers.hasCanonicalMerraTimeSupport                   - True for the complete MERRA time contract
%   icemodel.forcing.helpers.hasConstantMerraTavg3Support                   - True when glacier channels hold each UTC block
%   icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid                  - True for an exact native tavg3 inventory
%   icemodel.forcing.helpers.interpRcm                                      - Interpolate RCM grid-cell data to query points
%   icemodel.forcing.helpers.intervalMaximumSolarElevation                  - Maximum sun angle over each time support
%   icemodel.forcing.helpers.listSourceFiles                                - Return flat and nested files under a source/cache root
%   icemodel.forcing.helpers.locateImauHourlyFile                           - Find a flat or hourly-subfolder PANGAEA tab file
%   icemodel.forcing.helpers.locateKtransectFiles                           - Find every cached K-transect annual file for a station
%   icemodel.forcing.helpers.locateNormalizedFile                           - Find a file by normalized path/name token
%   icemodel.forcing.helpers.locateSamimiDye2Workbook                       - Find the Dye-2 summer 2016 AWS workbook
%   icemodel.forcing.helpers.marDiagnosticMetadata                          - Describe optional MAR mass-balance diagnostics
%   icemodel.forcing.helpers.marDynamicProfileQa                            - Diagnose MAR dynamic snow/firn layer consistency
%   icemodel.forcing.helpers.marRefreezeMetadata                            - Describe signed native MAR RZ values in one artifact
%   icemodel.forcing.helpers.metchecks                                      - Gap-fill and clamp met variables to physically valid ranges
%   icemodel.forcing.helpers.metfilename                                    - Build a standard icemodel met-file name
%   icemodel.forcing.helpers.metTimestepSuffix                              - Return the file tag for a model-met cadence
%   icemodel.forcing.helpers.metvariables                                   - Met-file variable names for the forcing builders
%   icemodel.forcing.helpers.modisAlbedoChannel                             - GEUS MODIS daily albedo on a time axis, per location
%   icemodel.forcing.helpers.normalizedFileToken                            - Compare filenames case-insensitively across separators
%   icemodel.forcing.helpers.normalizeGeusModisAlbedo                       - Mask undocumented GEUS albedo sentinels
%   icemodel.forcing.helpers.normalizeLocations                             - Normalize one location or a point list to a row cell
%   icemodel.forcing.helpers.optionalDate                                   - Normalize optional date-like inputs to UTC datetimes
%   icemodel.forcing.helpers.precipitationConsistency                       - Check finite, nonnegative phase mass balance
%   icemodel.forcing.helpers.precipitationValidity                          - Validate complete and partial precipitation splits
%   icemodel.forcing.helpers.precipitationVariables                         - Total and partitioned precipitation names
%   icemodel.forcing.helpers.projectLocation                                - Ensure a location struct carries EPSG:3413 coordinates
%   icemodel.forcing.helpers.promiceShortwave                               - Select physical public PROMICE shortwave channels
%   icemodel.forcing.helpers.pruneSupersededWindowFiles                     - Remove shorter windows contained by a new file
%   icemodel.forcing.helpers.psnProjection                                  - Polar stereographic north projection used by the builders
%   icemodel.forcing.helpers.readGcnetDonor                                 - Load one GC-Net surface file as an observed-only donor
%   icemodel.forcing.helpers.readImauHourlyTable                            - Parse one IMAU hourly PANGAEA AWS table
%   icemodel.forcing.helpers.readKtransectHeights                           - Parse the K-transect sensor-height workbook
%   icemodel.forcing.helpers.readKtransectTable                             - Parse one K-transect annual PANGAEA AWS table
%   icemodel.forcing.helpers.readMarDensitySnapshots                        - Read requested MAR RO1 density profiles
%   icemodel.forcing.helpers.readMerra2Time                                 - Decode only a MERRA-2 file's native UTC time coordinate
%   icemodel.forcing.helpers.readNetcdfAttribute                            - Return one NetCDF attribute or empty string when absent
%   icemodel.forcing.helpers.readPangaeaTab                                 - Read one PANGAEA tab-delimited dataset export
%   icemodel.forcing.helpers.regexpOnce                                     - Return one stripped regexp token or an empty string
%   icemodel.forcing.helpers.remapPolygon                                   - Conservative area-weighted remap of a grid block to a polygon
%   icemodel.forcing.helpers.resampleMetTimestep                            - Resample model met without crossing source outages
%   icemodel.forcing.helpers.slabMean                                       - Mean of a static grid field over hyperslab target cells
%   icemodel.forcing.helpers.solarElevation                                 - Approximate NOAA geometric solar elevation in degrees
%   icemodel.forcing.helpers.sourceAlbedo                                   - Derive broadband albedo where shortwave input is valid
%   icemodel.forcing.helpers.sourceSearchDirs                               - Candidate dirs for a staged file: per-source subfolder first
%   icemodel.forcing.helpers.sourceVariableAttributes                       - Preserve NetCDF variable attributes by variable
%   icemodel.forcing.helpers.stampMetadata                                  - Embed CF-ish metadata in a timetable's properties
%   icemodel.forcing.helpers.stationTransitionTimes                         - Within-record AWS handover times for a PROMICE site
%   icemodel.forcing.helpers.surfaceFlags                                   - Per-sample quality flags for a PROMICE surface-height series
%   icemodel.forcing.helpers.timeWindowMask                                 - Select an optional datetime window from a source time axis
%   icemodel.forcing.helpers.uniformCadenceSeconds                          - Derive one timetable's exact regular cadence
%   icemodel.forcing.helpers.validatemet                                    - Assert that MET satisfies the icemodel met-file contract
%   icemodel.forcing.helpers.variableUnits                                  - Unit string for each forcing-builder channel
%   icemodel.forcing.helpers.verificationSourceDir                          - Resolve repo-local verification source data roots
%   icemodel.forcing.helpers.windFromComponents                             - Wind speed and direction from zonal/meridional components
%   icemodel.forcing.helpers.writemet                                       - Validate and save an icemodel met file
%   icemodel.forcing.helpers.writeuserdata                                  - Save a Data timetable as met-swap userdata files
%
%   +ICEMODEL/+FORCING/+RECONSTRUCT
%   icemodel.forcing.reconstruct.acceptanceWindow                           - Per-site forcing-ready policy window from staged proxies
%   icemodel.forcing.reconstruct.admissionGate                              - Apply the approved per-variable admission thresholds
%   icemodel.forcing.reconstruct.applyDonorTransfer                         - Apply a fitted donor transfer to donor samples
%   icemodel.forcing.reconstruct.applyProxyCalibration                      - Apply a fitted proxy calibration to model samples
%   icemodel.forcing.reconstruct.assertNotEvaluationDestination             - Refuse reconstruction writes under eval data
%   icemodel.forcing.reconstruct.assertPromiceFilledArtifact                - Prove product and station provenance
%   icemodel.forcing.reconstruct.auditSegments                              - Build one audit row per contiguous selected segment
%   icemodel.forcing.reconstruct.blendFallbackSeams                         - Apply the policy seam taper to fallback fills
%   icemodel.forcing.reconstruct.blendSeams                                 - Taper excess anchored boundary mismatch across one run
%   icemodel.forcing.reconstruct.bucketEdges                                - Return the gap-duration bucket edges in hours
%   icemodel.forcing.reconstruct.clearSkyIndex                              - Normalize shortwave flux by station-specific TOA irradiance
%   icemodel.forcing.reconstruct.climatologyFill                            - Day-of-year climatology estimate of one channel
%   icemodel.forcing.reconstruct.commonSupportSkill                         - Compare candidate and baseline on identical samples
%   icemodel.forcing.reconstruct.deriveUpwardShortwave                      - Fill missing swu from final albedo and swd
%   icemodel.forcing.reconstruct.elevationAdjust                            - Adjust a donor channel across an elevation difference
%   icemodel.forcing.reconstruct.fillPromiceStation                         - Produce the gap-filled met product for one station
%   icemodel.forcing.reconstruct.fillShortGaps                              - Tier-1 bounded interior interpolation of one channel
%   icemodel.forcing.reconstruct.fillTwilightClimatology                    - Fill one-posting SWD gaps beside known night
%   icemodel.forcing.reconstruct.fitDonorTransfer                           - Fit an overlap-calibrated donor-to-target transfer
%   icemodel.forcing.reconstruct.fitProxyCalibration                        - Calibrate a model proxy on its observed overlap
%   icemodel.forcing.reconstruct.flatRunScreen                              - Flag multi-day buried/rime-encased sensor runs in met data
%   icemodel.forcing.reconstruct.gapCensus                                  - Census contiguous missing runs in one role-contract series
%   icemodel.forcing.reconstruct.gapDurationBucket                          - Assign positive durations to right-closed policy bins
%   icemodel.forcing.reconstruct.icemodelRequiredChannels                   - The POLICY A5 seven-channel icemodel set
%   icemodel.forcing.reconstruct.interpolationCapHours                      - Approved per-channel interpolation ceilings
%   icemodel.forcing.reconstruct.lastResortProxies                          - Adopt aligned proxy values for residual gaps
%   icemodel.forcing.reconstruct.loadWidestTimetable                        - Load the staged timetable with the widest time axis
%   icemodel.forcing.reconstruct.lwdEstimator                               - Empirical downward-longwave candidate from temperature and RH
%   icemodel.forcing.reconstruct.mustBeCapHours                             - Require a gap cap within any approved channel ceiling
%   icemodel.forcing.reconstruct.mustBeStationToken                         - Require canonical lowercase alphanumeric station IDs
%   icemodel.forcing.reconstruct.partitionPrecipitation                     - Split total precipitation by air temperature
%   icemodel.forcing.reconstruct.persistenceEstimate                        - Hold pre-gap values without held-out-data leakage
%   icemodel.forcing.reconstruct.physicalBounds                             - Return the approved physical bounds for one channel
%   icemodel.forcing.reconstruct.physicalValidity                           - Enforce scalar and relational reconstruction bounds
%   POLICY.md
%   icemodel.forcing.reconstruct.policySha256                               - Return the SHA-256 fingerprint of reconstruction POLICY.md
%   icemodel.forcing.reconstruct.promiceFilledVerificationMatches           - Match a prevalidated runtime identity
%   icemodel.forcing.reconstruct.provenanceCodes                            - Return the per-sample reconstruction provenance registry
%   icemodel.forcing.reconstruct.proxyArtifactIdentity                      - Verify one staged proxy's target and producer
%   README.md
%   icemodel.forcing.reconstruct.reconstructSeries                          - Compose admitted fill methods into one target series
%   icemodel.forcing.reconstruct.scalarValidity                             - Check finite samples against the A15 scalar registry
%   icemodel.forcing.reconstruct.seasonOf                                   - Meteorological season label of each timestamp
%   icemodel.forcing.reconstruct.selectedDataRoot                           - Resolve one selected met path to its data and met roots
%   icemodel.forcing.reconstruct.setopts                                    - Central options for the reconstruction pipeline
%   icemodel.forcing.reconstruct.smoothShortwaveSeams                       - Repair empirical outlier boundaries in filled SWD
%   icemodel.forcing.reconstruct.solarElevationBands                        - Solar-elevation thresholds for the swd science path
%   icemodel.forcing.reconstruct.stampGapfillIdentity                       - Stamp the gapfill_* identity fields on a filled met
%   icemodel.forcing.reconstruct.stationMethodPlan                          - Select and fit admitted fill methods for one station
%   icemodel.forcing.reconstruct.stepScale                                  - Per-season median absolute step of one observed channel
%   icemodel.forcing.reconstruct.syntheticMissingness                       - Draw blocked synthetic gaps into observed segments
%   icemodel.forcing.reconstruct.toaIrradiance                              - Top-of-atmosphere irradiance on a horizontal surface
%   icemodel.forcing.reconstruct.validationMetrics                          - Grade reconstructed samples against withheld truth
%   icemodel.forcing.reconstruct.validationSplit                            - Partition station years into selection and evaluation sets
%   icemodel.forcing.reconstruct.verifyPromiceFilledReadiness               - Gate derived PROMICE forcing by coverage
%
%   +ICEMODEL/+HELPERS
%   icemodel.helpers.absolutePath                                           - Return a path anchored at the current MATLAB folder
%   icemodel.helpers.canonicalPath                                          - Return the canonical absolute form of a file path
%   icemodel.helpers.copyFields                                             - Copy SOURCE struct fields onto TARGET
%   icemodel.helpers.ensureDirExists                                        - Create one directory when it does not already exist
%   README.md
%   icemodel.helpers.rmttleapinds
%
%   +ICEMODEL/+INTERNAL
%   icemodel.internal.contact                                               - Set or get the icemodel contact
%   icemodel.internal.fullpath                                              - Build full path to toolbox folder or file
%   icemodel.internal.installRequiredFiles                                  - Install required files from Github
%   icemodel.internal.isTestRun                                             - True when the MATLAB test framework is on the call stack
%   icemodel.internal.makecontents                                          - Write a Contents.m file in each namespace folder
%   icemodel.internal.readCffVersion                                        - Read the top-level software version from a CFF file
%   README.md
%   icemodel.internal.reference                                             - Set or get the icemodel reference
%   icemodel.internal.releaseMetadata                                       - Prepare, observe, or finalize release metadata
%   icemodel.internal.version                                               - Set or get the IceModel version number
%
%   +ICEMODEL/+KERNELS
%   icemodel.kernels.air_kinematic_viscosity                                - Approximate air kinematic viscosity
%   icemodel.kernels.buckVaporModel                                         - Buck (1981) vapor model, a self-contained reference archive
%   icemodel.kernels.latentEnthalpyWater                                    - Canonical Romps/Ambaum latent enthalpy reference
%   icemodel.kernels.potential_surface_vapor_demand                         - Convert latent heat to surface demand
%   README.md
%   icemodel.kernels.saturationVaporPressure                                - Canonical Romps/Ambaum saturation vapor pressure
%   icemodel.kernels.thermal_conductivity_air                               - Dry-air thermal conductivity and derivative
%   icemodel.kernels.thermal_conductivity_firn                              - Firn conductivity from temperature and density
%   icemodel.kernels.thermal_conductivity_ice                               - Thermal conductivity of ice from temperature
%   icemodel.kernels.thermal_conductivity_snow                              - Archived multi-option snow conductivity helper
%   icemodel.kernels.thermal_conductivity_water                             - Liquid-water thermal conductivity and derivative
%
%   +ICEMODEL/+NAMELISTS
%   icemodel.namelists.benchmark                                            - Return the supported formal benchmark file names
%   icemodel.namelists.benchmarksamplingprofile                             - Return the supported benchmark runner profiles
%   icemodel.namelists.budgetoutputs                                        - Return the mass-budget output channels
%   icemodel.namelists.completions                                          - Return the supported icemodel.completions selector names
%   icemodel.namelists.config                                               - Return the supported icemodel.config selector names
%   icemodel.namelists.cumulativeoutputs                                    - Return cumulative column diagnostic channels
%   icemodel.namelists.cvconvert                                            - Return the supported cvconvert dimension names
%   icemodel.namelists.forcings                                             - Return the supported forcing-source names
%   icemodel.namelists.getpath                                              - Return the supported icemodel.getpath path kinds
%   icemodel.namelists.physicalconstant                                     - Return the supported physical constant names
%   README.md
%   icemodel.namelists.rollingbaseline                                      - Return the supported mutable baseline selectors
%   icemodel.namelists.sitename                                             - Return the supported core run.point site names
%   icemodel.namelists.smbmodel                                             - Return supported smbmodel names by group
%   icemodel.namelists.solver                                               - Return the supported icemodel solver ids
%   icemodel.namelists.surfaceoutputs                                       - Return surface (ice1) output channel names
%   icemodel.namelists.testsmbmodel                                         - Return the supported formal test smbmodel selectors
%   icemodel.namelists.testtier                                             - Return the supported formal test suite tiers
%   icemodel.namelists.testverbosity                                        - Return the supported unit-runner verbosity names
%   icemodel.namelists.unittest                                             - Return the supported unit test file names
%   icemodel.namelists.userdata                                             - Return the supported core userdata source names
%   icemodel.namelists.uservars                                             - Return the supported core userdata variable names
%
%   +ICEMODEL/+NETCDF
%   icemodel.netcdf.config                                                  - Configure icemodel.netcdf API preferences
%   icemodel.netcdf.create                                                  - Create a new NetCDF file with the given properties and global
%   icemodel.netcdf.defdatavars                                             - Define the icemodel data variables and attributes
%   icemodel.netcdf.defdimid                                                - Define icemodel netcdf file dimensions
%   icemodel.netcdf.defdimvars                                              - Define icemodel netcdf grid and time dims and attributes
%   icemodel.netcdf.getchunksize
%   icemodel.netcdf.getdefaults
%   icemodel.netcdf.getdimdata                                              - Get dimensions of icemodel simulation data
%   icemodel.netcdf.getdimsize                                              - Return the size of each named dimension
%   icemodel.netcdf.getvardata                                              - Read icemodel data into memory and fill the arrays
%   icemodel.netcdf.getvarinfo                                              - Load one data file to get the shape and variable names
%   icemodel.netcdf.makencfile
%   icemodel.netcdf.maxcells                                                - Calculate the maximum number of gridcells for icemodel nc file
%   icemodel.netcdf.ncread                                                  - Read icemodel nc file into memory
%   icemodel.netcdf.nctype2mat                                              - Map NetCDF data types to MATLAB data types
%   README.md
%   icemodel.netcdf.redefatt
%   icemodel.netcdf.setfilename
%   icemodel.netcdf.writedims                                               - Write dimensions to icemodel nc file
%   icemodel.netcdf.writeice1                                               - Write ice1 data to an icemodel nc file
%   icemodel.netcdf.writeice2                                               - Write ice2 data to icemodel nc file
%
%   +ICEMODEL/+NETCDF/+DEFAULTS
%   icemodel.netcdf.defaults.axes                                           - Define the grid and time dimension units. Use empty char '' for variables
%   icemodel.netcdf.defaults.cfStandardNames                                - Load the official CF Standard Name Table as a set
%   icemodel.netcdf.defaults.longnames                                      - Define the grid and time dimension names
%   icemodel.netcdf.defaults.standardnames                                  - Not all values have standard names, so I constructed some
%   icemodel.netcdf.defaults.units                                          - Define the grid and time dimension units
%   icemodel.netcdf.defaults.variable                                       - Canonical {standard_name, long_name, unit, is_cf} for a channel
%   icemodel.netcdf.defaults.variables                                      - Canonical variable-metadata map for every icemodel channel
%   icemodel.netcdf.defaults.varnames                                       - Define the grid and time dimension names
%
%   +ICEMODEL/+NUMERICS
%   icemodel.numerics.aitkenscalar                                          - Apply scalar Aitken Delta-squared acceleration with safeguards
%   icemodel.numerics.complexstep                                           - Find a scalar root with a complex-step Newton iteration
%   icemodel.numerics.complexstep_derivative                                - Estimate a scalar derivative with a complex step
%   icemodel.numerics.fsearchzero                                           - Search for a scalar nonlinear root using bounded Brent fallback
%   README.md
%   icemodel.numerics.secantscalar                                          - Apply a safeguarded scalar secant step
%   icemodel.numerics.sign_or_one                                           - Return the sign of X, treating zero as positive one
%   icemodel.numerics.trisolve                                              - Solve a tridiagonal matrix equation Ax = b
%
%   +ICEMODEL/+PLOT
%   icemodel.plot.bulkextcoefs                                              - Icemodel.plot.bulk_extinction_coefficients
%   icemodel.plot.canonicalTimeDimension                                    - Enforce the icemodel timetable row-time name
%   icemodel.plot.compareTimeseries                                         - Overlay one variable from multiple timetables
%   icemodel.plot.enbal                                                     - Plot the surface energy balance
%   icemodel.plot.filterDateWindow                                          - Apply a shared inclusive date window to plot inputs
%   icemodel.plot.forcing                                                   - Plot one or more icemodel met-file forcings
%   icemodel.plot.formatDuration                                            - Render a duration in hours as a human-readable label
%   icemodel.plot.markTimeSpan                                              - Mark a time span on an axes without touching the legend
%   icemodel.plot.mesh                                                      - Plot the nodes and edges of an exponential grid
%   icemodel.plot.newFigure                                                 - Create a white, export-ready figure with a stable pixel size
%   icemodel.plot.parseDate                                                 - Normalize optional date-like input to a datetime
%   icemodel.plot.profile                                                   - Plot one or more depth-profile tables
%   README.md
%   icemodel.plot.scatterplot                                               - Plot paired data with 1:1 and fitted-line overlays
%   icemodel.plot.sourceColor                                               - Return stable verification colors keyed by source identity
%   icemodel.plot.temperature                                               - ICEMODEL.PLOT.TEMPERATURE
%   icemodel.plot.thermal_conductivity_air                                  - Plot dry-air thermal conductivity diagnostics
%   icemodel.plot.thermal_conductivity_firn                                 - Plot firn conductivity references
%   icemodel.plot.thermal_conductivity_ice                                  - Plot ice conductivity reference curves
%   icemodel.plot.thermal_conductivity_snow                                 - Plot archived snow conductivity options
%   icemodel.plot.thermal_conductivity_water                                - Plot liquid-water conductivity diagnostics
%   icemodel.plot.timeseries                                                - Plot one or more datetime-indexed series on one axes
%   icemodel.plot.vaporModel                                                - Three-formulation diagnostic comparison plot
%   icemodel.plot.variableUnit                                              - Return source table units with canonical metadata fallback
%
%   +ICEMODEL/+RADIATION
%   icemodel.radiation.bulk_extinction_coefficients                         - Compute bulk extinction coefficients
%   icemodel.radiation.get_scattering_coefficients                          - Extract spectral scattering coefficients from
%   icemodel.radiation.get_solar_spectrum                                   - Interpolate the reference solar spectrum to the model
%   icemodel.radiation.initialize_spectral_model                            - Initialize spectral geometry and coefficients
%   icemodel.radiation.load_spectral_tables                                 - Load raw spectral input tables from disk
%   icemodel.radiation.make_bulk_extinction_lookup                          - Precompute bulk-extinction coefficients
%   README.md
%   icemodel.radiation.rescale_spectral_extinction_coefficients             - Scale extinction by absorption
%   icemodel.radiation.smoothtwostream                                      - Smooth the up/down flux profiles from the two-stream solve
%   icemodel.radiation.solvetwostream                                       - Solve Schlatter's two-stream radiative transfer system
%   icemodel.radiation.spectral_extinction_coefficients                     - Compute spectral extinction coefficients
%   icemodel.radiation.update_extinction_coefficients                       - Update spectral extinction coefficients for
%
%   +ICEMODEL/+RUN
%   icemodel.run.point                                                      - Run one point-scale simulation and post-process its output
%   README.md
%
%   +ICEMODEL/+SURFACE
%   icemodel.surface.advective_heat_flux                                    - Compute heat advected to the surface by rainfall
%   icemodel.surface.apply_surface_vapor_exchange                           - Apply surface vapor energy demand
%   icemodel.surface.atmospheric_pressure_from_elevation                    - Estimate pressure from elevation
%   icemodel.surface.atmospheric_vapor_pressure                             - Relative humidity to atmospheric vapor pressure
%   icemodel.surface.conductive_heat_flux                                   - Conductive heat flux into the surface and derivative
%   icemodel.surface.diagnose_melt_freeze_energy                            - Diagnose surplus/deficit energy relative to Tf
%   icemodel.surface.diagnose_surface_ablation                              - Diagnose cumulative surface ablation terms
%   icemodel.surface.diagnose_surface_energy_balance                        - Diagnose the full surface energy budget
%   icemodel.surface.diagnose_surface_runoff                                - Diagnose cumulative runoff from surface fluxes
%   icemodel.surface.diagnose_turbulent_heat_fluxes                         - Run the configured THF scheme
%   icemodel.surface.dump_turbulent_heat_flux_debug_state                   - Save THF/SEB failure diagnostics
%   icemodel.surface.empirical_incoming_longwave_radiation                  - Estimate downwelling longwave
%   icemodel.surface.evaluate_surface_energy_balance                        - Evaluate the SEB from known flux terms
%   icemodel.surface.incoming_shortwave_radiation                           - Estimate downwelling shortwave radiation
%   icemodel.surface.initialize_surface_forcings                            - Load the meteorological forcing vectors
%   icemodel.surface.initialize_surface_state                               - Precompute forcing-derived surface state vectors
%   icemodel.surface.net_longwave_radiation                                 - Net surface longwave radiation and T_sfc derivative
%   icemodel.surface.net_shortwave_radiation                                - Compute net absorbed shortwave radiation
%   icemodel.surface.numerical_surface_flux                                 - Evaluate the SEB residual and derivative numerically
%   icemodel.surface.outgoing_longwave_radiation                            - Outgoing longwave radiation and T_sfc derivative
%   icemodel.surface.physical_surface_temperature                           - Cap surface temperature at the melting point
%   icemodel.surface.potential_surface_vapor_demand                         - Diagnose top-cell vapor energy demand
%   icemodel.surface.potential_surface_vapor_exchange                       - Partition surface vapor demand
%   README.md
%   icemodel.surface.resolve_forcing_snow_depth                             - Resolve scalar snow-depth for the THF scheme
%   icemodel.surface.solve_surface_energy_balance                           - Solve the nonlinear surface energy balance
%   icemodel.surface.solve_surface_temperature                              - Solve the explicit bulk-Richardson SEB for T_sfc
%   icemodel.surface.step_observation_heights                               - Select scalar observation heights for one step
%   icemodel.surface.surface_bulk_density                                   - Compute the bulk density of the top model layer
%   icemodel.surface.surface_energy_balance_residual                        - Return the SEB residual at T_sfc
%   icemodel.surface.surface_energy_balance_terms                           - Evaluate the SEB term set at T_sfc
%   icemodel.surface.surface_flux_linearization                             - Linearize the non-conductive surface flux
%   icemodel.surface.surface_roughness_length                               - Select the momentum roughness length z0m
%   icemodel.surface.surface_vapor_mass_flux                                - Convert a surface vapor fraction to a mass flux
%   icemodel.surface.terrain_adjusted_shortwave_radiation                   - Estimate terrain-adjusted shortwave
%   icemodel.surface.update_surface_state                                   - Update the surface state at substep entry
%
%   +ICEMODEL/+SURFACE/+TURBULENCE
%   README.md
%
%   +ICEMODEL/+SURFACE/+TURBULENCE/+BULK_RICHARDSON
%   icemodel.surface.turbulence.bulk_richardson.bulk_richardson_diagnostics - Assemble the full bulk-Richardson diagnostics
%   icemodel.surface.turbulence.bulk_richardson.exchange_coefficients       - Compute bulk-Richardson exchange coefficients
%   icemodel.surface.turbulence.bulk_richardson.latent_heat_flux            - Compute the turbulent latent heat flux
%   icemodel.surface.turbulence.bulk_richardson.richardson_number           - Compute the bulk Richardson number
%   icemodel.surface.turbulence.bulk_richardson.sensible_heat_flux          - Compute the turbulent sensible heat flux
%   icemodel.surface.turbulence.bulk_richardson.stability_factor            - Compute the stability function and derivative wrt T_sfc
%   icemodel.surface.turbulence.bulk_richardson.surface_flux_linearization  - Linearize the surface energy balance equation
%   icemodel.surface.turbulence.bulk_richardson.turbulent_heat_flux         - Evaluate the bulk-Richardson THF scheme
%
%   +ICEMODEL/+SURFACE/+TURBULENCE/+MONIN_OBUKHOV
%   icemodel.surface.turbulence.monin_obukhov.monin_obukhov_length          - Return the Monin-Obukhov stability length
%   icemodel.surface.turbulence.monin_obukhov.psi_h_paulson                 - Dyer/Paulson unstable scalar profile correction
%   icemodel.surface.turbulence.monin_obukhov.psi_holtslag                  - Holtslag and de Bruin stable profile correction
%   icemodel.surface.turbulence.monin_obukhov.psi_m_paulson                 - Paulson unstable momentum profile correction
%   icemodel.surface.turbulence.monin_obukhov.scalar_roughness_lengths      - Return scalar roughness lengths for bulk-MO
%   icemodel.surface.turbulence.monin_obukhov.stability_corrections         - Return Monin-Obukhov profile corrections
%   icemodel.surface.turbulence.monin_obukhov.surface_flux_linearization    - Linearize the bulk-MO surface flux
%   icemodel.surface.turbulence.monin_obukhov.turbulent_heat_flux           - Evaluate the Monin-Obukhov THF scheme
%
%   +ICEMODEL/+TEST
%   README.md
%
%   +ICEMODEL/+TEST/+FIXTURES
%   icemodel.test.fixtures.cleanupSyntheticWorkspace                        - Restore env vars and remove a temp workspace
%   icemodel.test.fixtures.makeReconstructSeries                            - One year of hourly synthetic met, smooth channels
%   icemodel.test.fixtures.makeSyntheticColumnState                         - Build a resolved synthetic column kernel state
%   icemodel.test.fixtures.makeSyntheticMetFile                             - Build a simple synthetic forcing timetable for tests
%   icemodel.test.fixtures.makeSyntheticWorkspace                           - Create an isolated icemodel test workspace
%   icemodel.test.fixtures.writeSyntheticMetFile                            - Write a synthetic met file for tests
%   icemodel.test.fixtures.writeSyntheticUserdataFile                       - Write a synthetic yearly userdata timetable
%
%   +ICEMODEL/+TEST/+HELPERS
%   icemodel.test.helpers.ambientAnchorVerdict                              - Decide whether ambient conditions held for a run
%   icemodel.test.helpers.archiveManagedBaseline                            - Archive a rolling baseline before overwrite
%   icemodel.test.helpers.artifactFilePath                                  - Return the canonical artifact file path
%   icemodel.test.helpers.assertAmbientBaselineAcceptance                   - Validate the final baseline anchor
%   icemodel.test.helpers.assertCleanPerfSession                            - Refuse an in-session formal run in a dirty session
%   icemodel.test.helpers.assertCleanSnapshotWorktree                       - Require a clean worktree for a release snapshot
%   icemodel.test.helpers.assertFormalBaselineCandidate                     - Reject incomplete state before publication
%   icemodel.test.helpers.assertFormalBaselineForcing                       - Verify a baseline's registered forcing identity
%   icemodel.test.helpers.assertFormalBenchmarkCandidate                    - Reject invalid component timing evidence
%   icemodel.test.helpers.assertNewReleaseBaselineTarget                    - Reject an existing immutable release file
%   icemodel.test.helpers.assertPerfBuildQuality                            - Require measurement quality before a managed build
%   icemodel.test.helpers.assertReleasePerfBaselineSource                   - Verify a rolling perf source for release
%   icemodel.test.helpers.baselineFilePath                                  - Return the canonical baseline file path
%   icemodel.test.helpers.baselineProfilerDir                               - Return the profiler-artifact folder for a baseline file
%   icemodel.test.helpers.benchmarkSuiteSignature                           - Hash the managed core benchmark suite
%   icemodel.test.helpers.bootstrapTestEnvironment                          - Add test paths and install one scoped data config
%   icemodel.test.helpers.buildSyntheticOpts                                - Build resolved OPTS for synthetic unit-test runs
%   icemodel.test.helpers.buildThfValidationCases                           - Build focused real-case THF validation cases
%   icemodel.test.helpers.captureBaselineProfile                            - Save a build-time profiler report alongside a baseline
%   icemodel.test.helpers.captureExpectedWarning                            - Run fcn once, verify its warning, and capture output
%   icemodel.test.helpers.commitBaselineProfilePublication                  - Remove a retained sidecar backup
%   icemodel.test.helpers.displayPerfResults                                - Display compact performance results from run_perf_suite
%   icemodel.test.helpers.displayPerfSummary                                - Display a compact perf case summary and benchmark table
%   icemodel.test.helpers.displayRegressionResults                          - Display compact results from run_regression_suite
%   icemodel.test.helpers.displayRegressionSummary                          - Display compact regression compare summaries
%   icemodel.test.helpers.findCaseRow                                       - Find the first baseline/report row matching CASE_ID
%   icemodel.test.helpers.findRunoffReferenceRow                            - Resolve runoff reference row for one formal case
%   icemodel.test.helpers.formalBaselinePolicy                              - Return forcing and default-root policy for a baseline
%   icemodel.test.helpers.formalPerformanceVerdict                          - Evaluate one formal timing comparison row
%   icemodel.test.helpers.formalRegressionMetricEvidence                    - Check one saved metric is comparable
%   icemodel.test.helpers.getFormalForcing                                  - Return the forcing identity for one formal baseline
%   icemodel.test.helpers.getFormalTestSuiteCases                           - Return the canonical formal test-suite cases
%   icemodel.test.helpers.getPerfCaseMatrix                                 - Return the canonical formal performance case matrix
%   icemodel.test.helpers.getRegressionCaseMatrix                           - Return the canonical formal regression case matrix
%   icemodel.test.helpers.getRunoffSite                                     - Map formal station cases to runoff-validation catchments
%   icemodel.test.helpers.loadArtifact                                      - Load a saved test artifact
%   icemodel.test.helpers.loadBaseline                                      - Load a rolling or release baseline table
%   icemodel.test.helpers.loadProcessedMetForOutputYears                    - Load processed met limited to output years
%   icemodel.test.helpers.loadReference                                     - Load a test reference table
%   icemodel.test.helpers.loadSavedTable                                    - Load a saved table-like object from a MAT file
%   icemodel.test.helpers.machineHostname                                   - Return a normalized, network-independent machine identity
%   icemodel.test.helpers.makeFormalCaseId                                  - Return canonical formal-suite identifier for one model run
%   icemodel.test.helpers.managedBaselineSiblings                           - Return the baseline files one build writes
%   icemodel.test.helpers.markTestSessionDirty                              - Record that this MATLAB session ran a test suite
%   icemodel.test.helpers.measurePerfCase                                   - Measure one case under the selected isolation protocol
%   icemodel.test.helpers.normalizeFormalCaseId                             - Normalize legacy formal-suite identifiers
%   icemodel.test.helpers.normalizeMachineIdentity                          - Fold a machine name to a stable identity
%   icemodel.test.helpers.perfBaselineCompatibility                         - Decide whether wall-time comparison is fair
%   icemodel.test.helpers.perfMeasurementPolicy                             - Formal timing gate thresholds
%   icemodel.test.helpers.perfMeasurementQuality                            - Judge a perf run on its measurement quality
%   icemodel.test.helpers.performanceGate                                   - Compare one runtime to a two-sided accepted band
%   icemodel.test.helpers.perfSampleValidity                                - Decide whether one case's timing samples are usable
%   icemodel.test.helpers.prepareBaselineBuild                              - Resolve shared setup for perf/regression baseline builds
%   icemodel.test.helpers.printFilePath                                     - Print a file path truncated to the test/ directory
%   icemodel.test.helpers.publishBaselineBundleSet                          - Publish a complete model baseline set
%   icemodel.test.helpers.publishBaselineProfile                            - Replace one managed profiler sidecar
%   icemodel.test.helpers.referenceFilePath                                 - Return the canonical reference file path
%   icemodel.test.helpers.regressionCaseGates                               - Evaluate every regression gate for one formal case
%   icemodel.test.helpers.regressionFailures                                - List the failed cases and gates of one regression run
%   icemodel.test.helpers.removeBaselineProfileStage                        - Remove one owned profiler staging directory
%   icemodel.test.helpers.removeReleaseSnapshotArtifacts                    - Remove a snapshot MAT file and sidecar
%   icemodel.test.helpers.resolveBaselineBuild                              - Resolve baseline type/tag and default output file
%   icemodel.test.helpers.resolveBaselineSelector                           - Parse rolling vs release baseline selectors
%   icemodel.test.helpers.resolveBootstrapRelease                           - Load or create registered release baselines
%   icemodel.test.helpers.resolveReleaseDataRoots                           - Resolve model and fixture roots for a baseline
%   icemodel.test.helpers.resolveRequestedSmbmodels                         - Expand one requested formal smbmodel selector
%   icemodel.test.helpers.resolveRunStamp                                   - Resolve shared batch run identifiers for test artifacts
%   icemodel.test.helpers.retryInvalidMeasurement                           - Measure once; re-measure once if invalid
%   icemodel.test.helpers.rollbackBaselineProfilePublication                - Restore the prior profiler sidecar
%   icemodel.test.helpers.runBenchmarkDiagnostics                           - Run and compare managed component benchmarks
%   icemodel.test.helpers.runModelCase                                      - Resolve, execute, and postprocess one supported model case
%   icemodel.test.helpers.runPerfCase                                       - Run one formal performance case and normalize the result
%   icemodel.test.helpers.runPerfCaseSubprocess                             - Run one formal perf case in this fresh session
%   icemodel.test.helpers.runSmbModel                                       - Dispatch to the requested core SMB model kernel
%   icemodel.test.helpers.runThfValidationCase                              - Run one real-case THF validation scenario
%   icemodel.test.helpers.sampleMachineState                                - Sample the machine conditions that disturb a timing
%   icemodel.test.helpers.sanitizeTag                                       - Replace punctuation and whitespace for filename-safe tags
%   icemodel.test.helpers.setModelOptsForCase                               - Build resolved model OPTS for one case
%   icemodel.test.helpers.smbmodelTag                                       - Return canonical smbmodel tag for filenames and identifiers
%   icemodel.test.helpers.snapshotBaseline                                  - Save a release snapshot from the rolling test baseline
%   icemodel.test.helpers.summarizeIce1Metrics                              - Extract formal regression metrics from output and refs
%   icemodel.test.helpers.summarizeMachineState                             - Reduce machine-state samples to one attestation
%   icemodel.test.helpers.testSessionActivity                               - Return the suites this MATLAB session has run
%   icemodel.test.helpers.transactionalSnapshotSet                          - Create and validate an aggregate release set
%   icemodel.test.helpers.worktreeRevision                                  - Return the Git description of the source tree
%
%   +ICEMODEL/+TEST/+VERIFY
%   icemodel.test.verify.verifyEqualNested                                  - Recursively compare nested structs, timetables, and arrays
%   icemodel.test.verify.verifyProcessedOutputBounds                        - Verify physical bounds of the processed output
%
%   +ICEMODEL/+TIMESTEPPING
%   icemodel.timestepping.acceptsubstep                                     - Accept the substep: checkpoint state and time bookkeeping
%   icemodel.timestepping.checksubstep                                      - Accept the substep or decide what the next attempt looks like
%   icemodel.timestepping.getforcings                                       - Load the forcing data for this timestep
%   icemodel.timestepping.getsubstepforcings                                - Return scalar forcing values for the current substep
%   icemodel.timestepping.initialize_timesteps                              - Initialize the model timestep counters
%   icemodel.timestepping.newtimestep                                       - Initialize forcing-step accumulators and diagnostics
%   icemodel.timestepping.nexttimestep                                      - Advance the forcing index and adapt the next substep size
%   README.md
%   icemodel.timestepping.resetsubstep                                      - Restore the accepted state and shorten the retry timestep
%
%   +ICEMODEL/+VALIDATORS
%   icemodel.validators.mustBeBenchmarkSamplingProfileName                  - Validate one benchmark profile name
%   icemodel.validators.mustBeForcingName                                   - Validate that input is a valid forcing-source name
%   icemodel.validators.mustBeFormalSmbmodelName                            - Validate one concrete formal-suite smbmodel
%   icemodel.validators.mustBeRollingBaselineName                           - Validate the mutable build-baseline selector
%   icemodel.validators.mustBeSiteName                                      - Validate that input is a valid core point-run site name
%   icemodel.validators.mustBeSmbmodelName                                  - Validate that input is a valid core smbmodel name
%   icemodel.validators.mustBeSolverFilter                                  - Validate an optional solver-filter vector
%   icemodel.validators.mustBeTestSmbmodelSelector                          - Validate one formal suite smbmodel selector
%   icemodel.validators.mustBeTestTierName                                  - Validate one formal suite tier selector
%   icemodel.validators.mustBeTestVerbosityName                             - Validate one unit-runner verbosity selector
%   icemodel.validators.mustBeUserdataName                                  - Validate that input is a valid userdata source name
%   icemodel.validators.mustBeUservarName                                   - Validate that input is a valid userdata variable name
%   README.md
%
%   +ICEMODEL/+VAPOR
%   icemodel.vapor.dew_point_temperature                                    - Dew point temperature from air temperature and
%   icemodel.vapor.initialize_vapor_model                                   - Initialize Ambaum (2020) Rankine-Kirchhoff vapor
%   icemodel.vapor.latent_enthalpy_switch                                   - Return Ls or Lv based on surface phase state
%   icemodel.vapor.moist_air_density                                        - Return moist-air density from partial pressures
%   README.md
%   icemodel.vapor.relative_humidity_from_specific_humidity                 - Specific humidity to RH [%]
%   icemodel.vapor.relative_humidity_from_vapor_pressure                    - Relative humidity from vapor pressure
%   icemodel.vapor.saturation_vapor_density                                 - Compute saturation vapor density in porous ice
%   icemodel.vapor.saturation_vapor_pressure                                - Saturation vapor pressure over liquid or ice
%   icemodel.vapor.specific_humidity_from_vapor_pressure                    - Convert vapor pressure to q
%   icemodel.vapor.vapor_diffusivity                                        - Effective water-vapor diffusion coefficient in porous ice
%   icemodel.vapor.vapor_pressure_from_specific_humidity                    - Convert specific humidity to vapor
%   icemodel.vapor.vapor_thermal_conductivity                               - Effective thermal conductivity from vapor
%   icemodel.vapor.wet_bulb_temperature                                     - Wet-bulb temperature from air temperature and RH
%
%   +ICEMODEL/+VERIFICATION
%   icemodel.verification.ablationPerformanceMetrics                        - Score every modeled ablation diagnostic
%   icemodel.verification.auditArtifacts                                    - Read-only QA/QC for manifest-referenced artifacts
%   icemodel.verification.candidateFromIcemodelOutput                       - Convert icemodel outputs for verification
%   icemodel.verification.compareAblation                                   - Compare PROMICE lowering with modeled solid-ice loss
%   icemodel.verification.comparecase                                       - Compare one staged verification target against a candidate
%   icemodel.verification.comparisonCompatibility                           - Derive staged verification comparison pairs
%   icemodel.verification.listcases                                         - Enumerate staged verification cases from family manifests
%   icemodel.verification.loadmanifest                                      - Return one resolved verification case manifest
%   MAR_DENSITY_PROFILES.md
%   icemodel.verification.matchObservations                                 - Match interval SMB and dated subsurface profiles
%   icemodel.verification.observationRateOutliers                           - Flag site-years whose observed ablation rate is far
%   icemodel.verification.plotcase                                          - Plot staged verification data without requiring model output
%   icemodel.verification.plotFirnArtifacts                                 - Plot staged firn-family artifacts for visual QA
%   icemodel.verification.plotscatter                                       - Plot target-versus-candidate scatter panels for site cases
%   icemodel.verification.plotVerificationArtifacts                         - Visualize staged verification artifacts
%   README.md
%   README_FIXTURES.md
%   icemodel.verification.runIcemodelCandidate                              - Run icemodel and return a verification candidate
%   icemodel.verification.syntheticSnowModelRun                             - Return snow-model-like icemodel outputs
%
%   +ICEMODEL/+VERIFICATION/+COLBECK
%   icemodel.verification.colbeck.analyticalSolution                        - Compute the analytical Colbeck infiltration solution
%   icemodel.verification.colbeck.caseDefinition                            - Return the canonical Colbeck 1976 verification case
%   icemodel.verification.colbeck.compareSolutions                          - Compare cached + computed Colbeck solutions side-by-side
%   icemodel.verification.colbeck.runCase                                   - Build a Colbeck candidate bundle for the verification suite
%
%   +ICEMODEL/+VERIFICATION/+HELPERS
%   icemodel.verification.helpers.ablationLedgerIncrements                  - Per-interval ablation terms from the mass ledger
%   icemodel.verification.helpers.alignObservationSeries                    - Align one observation/model series by its support
%   icemodel.verification.helpers.assertArtifactSha256                      - Require current bytes to match a pinned SHA-256
%   icemodel.verification.helpers.assertRootRelativeArtifactSha256          - Verify one root-scoped artifact identity
%   icemodel.verification.helpers.classifyObservationSupport                - Apply the flag-support rules to observation rows
%   icemodel.verification.helpers.classifySnowDepth                         - Classify snow support without accepting negative depth
%   icemodel.verification.helpers.esmRuntimeMetFiles                        - Resolve an atomic ESM case's standard runtime met paths
%   icemodel.verification.helpers.esmSnowmipWaterYear                       - Return one snow water year for an ESM-SnowMIP site
%   icemodel.verification.helpers.evaluationDataRoot                        - Resolve the base evaluation-data root
%   icemodel.verification.helpers.evaluationSeason                          - Summertime display and evaluation bounds for one year
%   icemodel.verification.helpers.familyManifestFiles                       - List staged verification family manifest files
%   icemodel.verification.helpers.fieldOr                                   - Return a struct field or default value
%   icemodel.verification.helpers.inputDataRoot                             - Resolve the base icemodel input-data root
%   icemodel.verification.helpers.isPhysicsFingerprint                      - Return true for one well-formed physics stamp
%   icemodel.verification.helpers.loadArtifact                              - Load one named staged verification artifact from a MAT file
%   icemodel.verification.helpers.loadColocatedData                         - Assemble a timeseries bundle from staged per-source files
%   icemodel.verification.helpers.metricRowSchema                           - Return canonical comparison-metric field names/defaults
%   icemodel.verification.helpers.observationSupportFields                  - Columns the observation support rules read
%   icemodel.verification.helpers.physicsFingerprint                        - Fingerprint the model's default physics configuration
%   icemodel.verification.helpers.profileGroups                             - Split profile rows by stable source identity and UTC date
%   icemodel.verification.helpers.readFamilyManifest                        - Read one verification family manifest JSON file
%   icemodel.verification.helpers.residualMetrics                           - Bias, MAE, RMSE, max error, and NSE for one paired series
%   icemodel.verification.helpers.resolveCandidateBundle                    - Resolve the comparison bundle for one case
%   icemodel.verification.helpers.sampleQuantile                            - Linear-interpolated sample quantile, no extra toolboxes
%   icemodel.verification.helpers.sumupColocation                           - Flag a SUMup point as co-located with a mixed anchor
%   icemodel.verification.helpers.validateAblationModelSchema               - Check a saved cohort still supports the report
%   icemodel.verification.helpers.writeRunReport                            - Write a concise markdown report for a verification run
%
%   +ICEMODEL/+VERIFICATION/+NAMELISTS
%   icemodel.verification.namelists.ablationReportChannels                  - Model channels read by the ablation report
%   icemodel.verification.namelists.caseid                                  - Return supported runnable snow-verification case ids
%   icemodel.verification.namelists.casetype                                - Return the supported snow-verification case types
%   icemodel.verification.namelists.completions                             - Return the supported verification namelist selector names
%   icemodel.verification.namelists.datasetfamily                           - Return the supported snow-verification dataset families
%   icemodel.verification.namelists.evaltarget                              - Return the supported case-manifest eval-target descriptors
%   icemodel.verification.namelists.firndatasetfamily                       - Return verification families used by firn staging previews
%   icemodel.verification.namelists.laughtests                              - Canonical Laugh-Tests case-id namelist
%   icemodel.verification.namelists.permafrostzone                          - Return the supported case-manifest permafrost-zone values
%   icemodel.verification.namelists.promiceAblationPolicy                   - Return the fixed PROMICE ablation comparison policy
%   icemodel.verification.namelists.promiceAblationReadiness                - Return the PROMICE ablation admission policy
%   icemodel.verification.namelists.promicesite                             - Auto-discovered PROMICE station-id namelist
%   icemodel.verification.namelists.rcmMetSources                           - Verification RCM labels that currently write met files
%   icemodel.verification.namelists.rcmProductIds                           - Map RCM runtime/storage labels to explicit product ids
%   icemodel.verification.namelists.rcmsources                              - Verification RCM source labels in canonical staging order
%   icemodel.verification.namelists.snowmipsite                             - Canonical ESM-SnowMIP site-name namelist
%   icemodel.verification.namelists.surfacezone                             - Return the supported case-manifest surface-zone values
%
%   +ICEMODEL/+VERIFICATION/+REPORT
%   build_snow_artifact_qa.py
%   icemodel.verification.report.buildAblationEvaluationReport              - Render saved PROMICE ablation results
%   icemodel.verification.report.buildGapFillReport                         - Build the gap-fill before/after Quarto report inputs
%   icemodel.verification.report.buildTestSuiteReport                       - Render numerical or performance suite results
%   check_snow_artifact_qa.py
%   icemodel.verification.report.configureCategoryAxis                      - Label horizontal evidence rows without clipping
%   icemodel.verification.report.escapeMarkdownText                         - Show saved text literally, without markup or raw HTML
%   icemodel.verification.report.exportAndClose                             - Export one report figure and release graphics state
%   icemodel.verification.report.formatReportAxes                           - Isolate exported graphics from interactive theme defaults
%   icemodel.verification.report.formatValue                                - Format one scalar table value for Markdown
%   icemodel.verification.report.gapfillFigureStyle                         - Define the gap-fill report figure colors
%   icemodel.verification.report.generateFinalFirnPreview                   - Build canonical firn QA and figure products
%   icemodel.verification.report.generateFinalSnowPreview                   - Build canonical seasonal QA, figures, and readiness
%   icemodel.verification.report.markdownCode                               - Wrap saved metadata in a code span that renders literally
%   icemodel.verification.report.markdownTable                              - Convert a compact table to inert Markdown
%   icemodel.verification.report.methodFillLayers                           - Split one channel into observed/own-fill/other-fill layers
%   README.md
%   icemodel.verification.report.safeLabel                                  - Collapse control characters for MATLAB graphics text
%   icemodel.verification.report.sanitizeText                               - Collapse controls and neutralize raw HTML delimiters
%   snow-artifact-qa.qmd
%
%   +ICEMODEL/+VERIFICATION/+SETUP
%   icemodel.verification.setup.anchorColocation                            - Flag a point as co-located with the nearest anchor
%   icemodel.verification.setup.buildDatasetFamilyManifest                  - Build and merge-write a family manifest
%   icemodel.verification.setup.buildEsmSnowmipForcing                      - Convert ESM-SnowMIP NetCDF to native forcing
%   icemodel.verification.setup.buildEsmSnowmipObservations                 - Convert ESM-SnowMIP obs NetCDF to verification targets
%   icemodel.verification.setup.buildFetchProductStatus                     - Build ordered status rows from a product registry
%   icemodel.verification.setup.buildLaughTestsArtifacts                    - Build one Laugh-Tests evaluation/reference bundle
%   icemodel.verification.setup.buildSumupObservations                      - Convert SUMup firn records to verification targets
%   icemodel.verification.setup.bytesSha256                                 - Return the lowercase hex SHA-256 of a byte vector
%   icemodel.verification.setup.caseManifestFieldNames                      - Return canonical case-manifest fields
%   icemodel.verification.setup.colocationSourceLists                       - Derive manifest source lists from colocation legs
%   icemodel.verification.setup.datasetFamilyStagingPaths                   - Build the shared dataset-family output paths
%   icemodel.verification.setup.deduplicateSumupRecords                     - Keep one row per SUMup scientific identity
%   icemodel.verification.setup.emptyFetchProductStatusRow                  - Return the literal shared fetch-row prototype
%   icemodel.verification.setup.ensureUtc                                   - Coerce a date/datetime input to a UTC-tagged datetime
%   icemodel.verification.setup.esmSnowmipSiteCatalog                       - Return the ESM-SnowMIP source-site catalog
%   icemodel.verification.setup.familyManifestFieldNames                    - Return canonical family-manifest fields
%   icemodel.verification.setup.fetchEsmSnowmip                             - Locate or verify the ESM-SnowMIP source NetCDF files
%   icemodel.verification.setup.fetchFixtures                               - Transactionally provision or verify release data
%   icemodel.verification.setup.fetchGcnet                                  - Locate or verify local Vandecrux/GC-Net source caches
%   icemodel.verification.setup.fetchImau                                   - Locate or verify local IMAU PANGAEA source caches
%   icemodel.verification.setup.fetchKtransect                              - Locate or verify the local K-transect PANGAEA source cache
%   icemodel.verification.setup.fetchLaughTests                             - Locate or verify the Laugh-Tests source checkout
%   icemodel.verification.setup.fetchMissingStatus                          - Convert missing source patterns to fetch status rows
%   icemodel.verification.setup.fetchProductFiles                           - Return files matching product-cache patterns
%   icemodel.verification.setup.fetchProductNames                           - Return ordered product selectors from a fetch registry
%   icemodel.verification.setup.fetchProductStatusRow                       - Build the standard fetch product status record
%   icemodel.verification.setup.fetchPromice                                - Locate or verify local PROMICE pypromice L3 source caches
%   icemodel.verification.setup.fetchRetmip                                 - Locate or verify local RetMIP source caches
%   icemodel.verification.setup.fetchSumup                                  - Locate or verify the SUMup firn source files
%   icemodel.verification.setup.fileSha256                                  - Return the lowercase hex SHA-256 of a file's bytes
%   icemodel.verification.setup.finishFetchStatus                           - Apply the shared fetch status/strict/silent contract
%   icemodel.verification.setup.firnCaseManifestFieldNames                  - Return canonical firn case-manifest fields
%   icemodel.verification.setup.fixtureCallerSymlink                        - Find the first caller-controlled link in a path
%   icemodel.verification.setup.fixtureDataRoot                             - Return the data root registered for a release
%   icemodel.verification.setup.fixtureFetchCommand                         - Build a copy-paste release-data repair command
%   icemodel.verification.setup.fixtureFileList                             - Return manifest paths for selected release capabilities
%   icemodel.verification.setup.fixtureRelativePosix                        - Return one root-relative path with POSIX separators
%   icemodel.verification.setup.formatManifestTime                          - Serialize manifest timestamps with explicit clock time
%   icemodel.verification.setup.gcnetInventory                              - Index Vandecrux/GC-Net products without loading arrays
%   icemodel.verification.setup.gcnetProductNames                           - Return canonical Vandecrux/GC-Net product selectors
%   icemodel.verification.setup.gcnetProductSpec                            - Return Vandecrux/GC-Net DOI metadata and file rules
%   icemodel.verification.setup.hasObservationRecords                       - True when any observation sub-bundle carries rows
%   icemodel.verification.setup.imauSiteCatalog                             - Return the IMAU hourly-AWS source-site catalog
%   icemodel.verification.setup.importEsmSnowmip                            - Stage ESM-SnowMIP site fixtures for the verification suite
%   icemodel.verification.setup.importImau                                  - Stage the IMAU hourly AWS verification family
%   icemodel.verification.setup.importKtransect                             - Stage the K-transect annual AWS verification family
%   icemodel.verification.setup.importLaughTests                            - Stage selected Laugh-Tests synthetic snow benchmarks
%   icemodel.verification.setup.importPromiceSites                          - Stage PROMICE-anchored firn-evaluation cases
%   icemodel.verification.setup.importResearchSites                         - Stage generic research-site firn targets
%   icemodel.verification.setup.importRetmip                                - Stage the RetMIP protocol verification family
%   icemodel.verification.setup.importSumup                                 - Stage co-located SUMup firn evaluation cases
%   icemodel.verification.setup.ktransectAliasCrosswalk                     - Return the K-transect station alias hypothesis table
%   icemodel.verification.setup.ktransectSiteCatalog                        - Return the K-transect AWS source-site catalog
%   icemodel.verification.setup.loadPriorDatasetFamilyCases                 - Read cases needed by an additive native refresh
%   icemodel.verification.setup.makeCaseManifestEntry                       - Build one case manifest entry from canonical fields
%   icemodel.verification.setup.makeFamilyManifest                          - Build one verification family manifest struct
%   icemodel.verification.setup.makeFirnCaseManifestEntry                   - Build one firn case manifest entry
%   icemodel.verification.setup.manifestWindow                              - Serialize one start/end pair for a JSON manifest
%   icemodel.verification.setup.mergeColocation                             - Copy every field from ADD onto a colocation struct
%   icemodel.verification.setup.metadataStruct                              - Build a metadata struct from a 2-column cell array
%   icemodel.verification.setup.metArtifactReadiness                        - Diagnose the exact saved scalar-window met artifact
%   icemodel.verification.setup.metForcingReady                             - Test unfilled readiness and inventory complete windows
%   icemodel.verification.setup.mixedAnchorCatalog                          - Read staged firn/research anchors from manifests
%   icemodel.verification.setup.mustBeKnownFetchProducts                    - Reject selectors absent from a fetch registry
%   icemodel.verification.setup.normalizeForcingSources                     - Normalize one public forcing-source selection
%   icemodel.verification.setup.packFixtures                                - Pack selected release-data capabilities into separate archives
%   icemodel.verification.setup.periodBounds                                - Parse a manifest period to UTC datetimes
%   icemodel.verification.setup.preferPrimary                               - Coalesce two same-shape series, keeping primary where finite
%   icemodel.verification.setup.prepareCaseRoot                             - Create one case folder and select artifacts needing writes
%   icemodel.verification.setup.prepareReplacementCaseEntry                 - Retain observations but clear prior runtime legs
%   icemodel.verification.setup.preservePriorNativeLeg                      - Merge prior native artifacts into a fresh case leg
%   icemodel.verification.setup.preserveRcmLegs                             - Keep compatible staged RCM legs after a failed refresh
%   icemodel.verification.setup.previewFirnStaging                          - Stage short build_forcing=true firn QA previews
%   icemodel.verification.setup.printFetchProductBanner                     - Print shared manual-cache retrieval instructions
%   icemodel.verification.setup.priorCaseById                               - Return one prior manifest case by canonical case id
%   icemodel.verification.setup.promiceSiteCatalog                          - Return the PROMICE source-site catalog
%   promote_snow_verification_artifacts.py
%   icemodel.verification.setup.rcmArtifactOutputDirs                       - Resolve shared default RCM artifact output roots
%   icemodel.verification.setup.rcmSourceCoverage                           - Probe on-disk year coverage of each forcing source
%   icemodel.verification.setup.rcmStorageAlias                             - Return the collision-safe RCM artifact identity for a case
%   icemodel.verification.setup.readBestSnowDepth                           - Site-aware snow-depth selector for ESM-SnowMIP obs
%   icemodel.verification.setup.readNetcdfTime                              - Read a NetCDF time coordinate as UTC datetime
%   icemodel.verification.setup.readNetcdfVariable                          - Read one ESM-SnowMIP-style NetCDF variable with NaN-fill
%   icemodel.verification.setup.readRcmArtifactMetadata                     - Read saved RCM provenance without payload arrays
%   icemodel.verification.setup.readRetmipProfileTable                      - Read a RetMIP initial profile table
%   icemodel.verification.setup.readRetmipProtocolTable                     - Read a RetMIP tab-delimited protocol time series
%   icemodel.verification.setup.refreshManifestSourceLists                  - Recompute source lists without rebuilding data
%   icemodel.verification.setup.refreshPromiceMetIdentities                 - Pin existing native met bytes in the manifest
%   icemodel.verification.setup.regexpOnce                                  - Return one stripped regexp token or an empty string
%   icemodel.verification.setup.releaseManifestFile                         - Return the release-data manifest path for a version
%   icemodel.verification.setup.relpaths                                    - Reduce absolute staged paths to base-relative names for JSON
%   icemodel.verification.setup.repairMetTimeSupport                        - Repair legacy linear 15-minute met artifacts
%   icemodel.verification.setup.repairRcmArtifactMetadata                   - Classify and repair current-token RCM artifacts
%   icemodel.verification.setup.reportPromiceCoverage                       - Print requested-vs-actual source coverage
%   icemodel.verification.setup.researchSiteCatalog                         - Return the catchall research source-site catalog
%   icemodel.verification.setup.resolveFetchCacheDir                        - Apply the shared fetch-cache defaulting rule
%   icemodel.verification.setup.resolveLegWindows                           - Decouple each gridded RCM leg's window from the met window
%   icemodel.verification.setup.resolveStagingRoots                         - Resolve paired eval/input staging roots for importers
%   icemodel.verification.setup.retmipCaseCatalog                           - Return the RetMIP protocol-case source catalog
%   icemodel.verification.setup.retmipOutputInventory                       - Return variables in a RetMIP model-output NetCDF
%   icemodel.verification.setup.reuseDatasetFamilyCases                     - Load staged cases for forcing-only attachment
%   icemodel.verification.setup.runDatasetFamilyDryRun                      - Build a dry-run manifest through shared orchestration
%   icemodel.verification.setup.runDatasetFamilyImport                      - Persist native state and optional requested RCMs
%   icemodel.verification.setup.selectSiteCatalogEntries                    - Select known ids from a source-site catalog
%   icemodel.verification.setup.stageDatasetFamilyCases                     - Stage requested cases with shared skip handling
%   icemodel.verification.setup.stageDatasetRcmForcing                      - Delegate one RCM source at a time to stageRcmForcing
%   icemodel.verification.setup.stageMarDensityProfiles                     - Add optional MAR RO1 profiles to SUMup cases
%   icemodel.verification.setup.stageRcmForcing                             - Stage RCM forcing + Data for a list of points
%   icemodel.verification.setup.stampArtifactMetadata                       - Add variable metadata to staged artifact tables
%   icemodel.verification.setup.stateCaseEntry                              - Refresh one staged state record into a manifest case entry
%   icemodel.verification.setup.sumupCacheDir                               - Resolve the canonical SUMup verification source cache
%   icemodel.verification.setup.sumupComparisonVariables                    - List nonempty SUMup observation groups
%   icemodel.verification.setup.textSha256                                  - Return the lowercase hex SHA-256 of a text value's UTF-8 bytes
%   icemodel.verification.setup.validateEvalTarget                          - Validate a case-manifest eval_target value
%   icemodel.verification.setup.validatePermafrostZone                      - Validate a case-manifest permafrost_zone value
%   icemodel.verification.setup.validateSurfaceZone                         - Validate a case-manifest surface_zone value
%   icemodel.verification.setup.writeFamilyManifestMerge                    - Merge new case entries into a family manifest
%   icemodel.verification.setup.writeJson                                   - Write pretty-printed JSON as UTF-8 with one trailing newline
%   icemodel.verification.setup.writeManifest                               - Write one verification family manifest JSON file
%   icemodel.verification.setup.writePromiceAblationReadiness               - Write the PROMICE ablation readiness ledger
%   icemodel.verification.setup.writePromiceSnowModelReadyYears             - Write the annual PROMICE snow handoff
%
%   +ICEMODEL/+VERIFICATION/+VALIDATORS
%   icemodel.verification.validators.mustBeCaseIdSubset                     - Validate a subset of canonical snow-verification case ids
%   icemodel.verification.validators.mustBeDatasetFamilyFilter              - Validate one optional dataset-family selector
%   icemodel.verification.validators.mustBeDatasetFamilySelection           - Validate dataset-family selectors plus "all"
%   icemodel.verification.validators.mustBeFirnDatasetFamilySelection       - Validate firn-family selectors plus "all"
%   icemodel.verification.validators.mustBeGcnetProductSelection            - Validate Vandecrux/GC-Net product selectors
%   icemodel.verification.validators.mustBeLaughTestCase                    - Validate cases against the Laugh-Tests namelist
%   icemodel.verification.validators.mustBeRcmSourceSelection               - Validate verification forcing-source selectors
%   icemodel.verification.validators.mustBeSnowmipSite                      - Validate sitename against the canonical ESM-SnowMIP
%
%   updatecontents.m generated this file on 18 Sep 2026 at 01:12:17.
