% +HELPERS
%
%   Contents file for +HELPERS and its subfolders.
%
%   +HELPERS
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
%   updatecontents.m generated this file on 14 Sep 2026 at 20:13:58.
