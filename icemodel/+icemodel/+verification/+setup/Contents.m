% +SETUP
%
%   Contents file for +SETUP and its subfolders.
%
%   +SETUP
%   icemodel.verification.setup.anchorColocation                - Flag a point as co-located with the nearest anchor
%   icemodel.verification.setup.buildDatasetFamilyManifest      - Build and merge-write a family manifest
%   icemodel.verification.setup.buildEsmSnowmipForcing          - Convert ESM-SnowMIP NetCDF to native forcing
%   icemodel.verification.setup.buildEsmSnowmipObservations     - Convert ESM-SnowMIP obs NetCDF to verification targets
%   icemodel.verification.setup.buildFetchProductStatus         - Build ordered status rows from a product registry
%   icemodel.verification.setup.buildLaughTestsArtifacts        - Build one Laugh-Tests evaluation/reference bundle
%   icemodel.verification.setup.buildSumupObservations          - Convert SUMup firn records to verification targets
%   icemodel.verification.setup.bytesSha256                     - Return the lowercase hex SHA-256 of a byte vector
%   icemodel.verification.setup.caseManifestFieldNames          - Return canonical case-manifest fields
%   icemodel.verification.setup.colocationSourceLists           - Derive manifest source lists from colocation legs
%   icemodel.verification.setup.datasetFamilyStagingPaths       - Build the shared dataset-family output paths
%   icemodel.verification.setup.deduplicateSumupRecords         - Keep one row per SUMup scientific identity
%   icemodel.verification.setup.emptyFetchProductStatusRow      - Return the literal shared fetch-row prototype
%   icemodel.verification.setup.ensureUtc                       - Coerce a date/datetime input to a UTC-tagged datetime
%   icemodel.verification.setup.esmSnowmipSiteCatalog           - Return the ESM-SnowMIP source-site catalog
%   icemodel.verification.setup.familyManifestFieldNames        - Return canonical family-manifest fields
%   icemodel.verification.setup.fetchEsmSnowmip                 - Locate or verify the ESM-SnowMIP source NetCDF files
%   icemodel.verification.setup.fetchFixtures                   - Transactionally provision or verify release data
%   icemodel.verification.setup.fetchGcnet                      - Locate or verify local Vandecrux/GC-Net source caches
%   icemodel.verification.setup.fetchImau                       - Locate or verify local IMAU PANGAEA source caches
%   icemodel.verification.setup.fetchKtransect                  - Locate or verify the local K-transect PANGAEA source cache
%   icemodel.verification.setup.fetchLaughTests                 - Locate or verify the Laugh-Tests source checkout
%   icemodel.verification.setup.fetchMissingStatus              - Convert missing source patterns to fetch status rows
%   icemodel.verification.setup.fetchProductFiles               - Return files matching product-cache patterns
%   icemodel.verification.setup.fetchProductNames               - Return ordered product selectors from a fetch registry
%   icemodel.verification.setup.fetchProductStatusRow           - Build the standard fetch product status record
%   icemodel.verification.setup.fetchPromice                    - Locate or verify local PROMICE pypromice L3 source caches
%   icemodel.verification.setup.fetchRetmip                     - Locate or verify local RetMIP source caches
%   icemodel.verification.setup.fetchSumup                      - Locate or verify the SUMup firn source files
%   icemodel.verification.setup.fileSha256                      - Return the lowercase hex SHA-256 of a file's bytes
%   icemodel.verification.setup.finishFetchStatus               - Apply the shared fetch status/strict/silent contract
%   icemodel.verification.setup.firnCaseManifestFieldNames      - Return canonical firn case-manifest fields
%   icemodel.verification.setup.fixtureCallerSymlink            - Find the first caller-controlled link in a path
%   icemodel.verification.setup.fixtureDataRoot                 - Return the data root registered for a release
%   icemodel.verification.setup.fixtureFetchCommand             - Build a copy-paste release-data repair command
%   icemodel.verification.setup.fixtureFileList                 - Return manifest paths for selected release capabilities
%   icemodel.verification.setup.fixtureRelativePosix            - Return one root-relative path with POSIX separators
%   icemodel.verification.setup.formatManifestTime              - Serialize manifest timestamps with explicit clock time
%   icemodel.verification.setup.gcnetInventory                  - Index Vandecrux/GC-Net products without loading arrays
%   icemodel.verification.setup.gcnetProductNames               - Return canonical Vandecrux/GC-Net product selectors
%   icemodel.verification.setup.gcnetProductSpec                - Return Vandecrux/GC-Net DOI metadata and file rules
%   icemodel.verification.setup.hasObservationRecords           - True when any observation sub-bundle carries rows
%   icemodel.verification.setup.imauSiteCatalog                 - Return the IMAU hourly-AWS source-site catalog
%   icemodel.verification.setup.importEsmSnowmip                - Stage ESM-SnowMIP site fixtures for the verification suite
%   icemodel.verification.setup.importImau                      - Stage the IMAU hourly AWS verification family
%   icemodel.verification.setup.importKtransect                 - Stage the K-transect annual AWS verification family
%   icemodel.verification.setup.importLaughTests                - Stage selected Laugh-Tests synthetic snow benchmarks
%   icemodel.verification.setup.importPromiceSites              - Stage PROMICE-anchored firn-evaluation cases
%   icemodel.verification.setup.importResearchSites             - Stage generic research-site firn targets
%   icemodel.verification.setup.importRetmip                    - Stage the RetMIP protocol verification family
%   icemodel.verification.setup.importSumup                     - Stage co-located SUMup firn evaluation cases
%   icemodel.verification.setup.ktransectAliasCrosswalk         - Return the K-transect station alias hypothesis table
%   icemodel.verification.setup.ktransectSiteCatalog            - Return the K-transect AWS source-site catalog
%   icemodel.verification.setup.loadPriorDatasetFamilyCases     - Read cases needed by an additive native refresh
%   icemodel.verification.setup.makeCaseManifestEntry           - Build one case manifest entry from canonical fields
%   icemodel.verification.setup.makeFamilyManifest              - Build one verification family manifest struct
%   icemodel.verification.setup.makeFirnCaseManifestEntry       - Build one firn case manifest entry
%   icemodel.verification.setup.manifestWindow                  - Serialize one start/end pair for a JSON manifest
%   icemodel.verification.setup.mergeColocation                 - Copy every field from ADD onto a colocation struct
%   icemodel.verification.setup.metadataStruct                  - Build a metadata struct from a 2-column cell array
%   icemodel.verification.setup.metArtifactReadiness            - Diagnose the exact saved scalar-window met artifact
%   icemodel.verification.setup.metForcingReady                 - Test unfilled readiness and inventory complete windows
%   icemodel.verification.setup.mixedAnchorCatalog              - Read staged firn/research anchors from manifests
%   icemodel.verification.setup.mustBeKnownFetchProducts        - Reject selectors absent from a fetch registry
%   icemodel.verification.setup.normalizeForcingSources         - Normalize one public forcing-source selection
%   icemodel.verification.setup.packFixtures                    - Pack selected release-data capabilities into separate archives
%   icemodel.verification.setup.periodBounds                    - Parse a manifest period to UTC datetimes
%   icemodel.verification.setup.preferPrimary                   - Coalesce two same-shape series, keeping primary where finite
%   icemodel.verification.setup.prepareCaseRoot                 - Create one case folder and select artifacts needing writes
%   icemodel.verification.setup.prepareReplacementCaseEntry     - Retain observations but clear prior runtime legs
%   icemodel.verification.setup.preservePriorNativeLeg          - Merge prior native artifacts into a fresh case leg
%   icemodel.verification.setup.preserveRcmLegs                 - Keep compatible staged RCM legs after a failed refresh
%   icemodel.verification.setup.previewFirnStaging              - Stage short build_forcing=true firn QA previews
%   icemodel.verification.setup.printFetchProductBanner         - Print shared manual-cache retrieval instructions
%   icemodel.verification.setup.priorCaseById                   - Return one prior manifest case by canonical case id
%   icemodel.verification.setup.promiceSiteCatalog              - Return the PROMICE source-site catalog
%   promote_snow_verification_artifacts.py
%   icemodel.verification.setup.rcmArtifactOutputDirs           - Resolve shared default RCM artifact output roots
%   icemodel.verification.setup.rcmSourceCoverage               - Probe on-disk year coverage of each forcing source
%   icemodel.verification.setup.rcmStorageAlias                 - Return the collision-safe RCM artifact identity for a case
%   icemodel.verification.setup.readBestSnowDepth               - Site-aware snow-depth selector for ESM-SnowMIP obs
%   icemodel.verification.setup.readNetcdfTime                  - Read a NetCDF time coordinate as UTC datetime
%   icemodel.verification.setup.readNetcdfVariable              - Read one ESM-SnowMIP-style NetCDF variable with NaN-fill
%   icemodel.verification.setup.readRcmArtifactMetadata         - Read saved RCM provenance without payload arrays
%   icemodel.verification.setup.readRetmipProfileTable          - Read a RetMIP initial profile table
%   icemodel.verification.setup.readRetmipProtocolTable         - Read a RetMIP tab-delimited protocol time series
%   icemodel.verification.setup.refreshManifestSourceLists      - Recompute source lists without rebuilding data
%   icemodel.verification.setup.refreshPromiceMetIdentities     - Pin existing native met bytes in the manifest
%   icemodel.verification.setup.regexpOnce                      - Return one stripped regexp token or an empty string
%   icemodel.verification.setup.releaseManifestFile             - Return the release-data manifest path for a version
%   icemodel.verification.setup.relpaths                        - Reduce absolute staged paths to base-relative names for JSON
%   icemodel.verification.setup.repairMetTimeSupport            - Repair legacy linear 15-minute met artifacts
%   icemodel.verification.setup.repairRcmArtifactMetadata       - Classify and repair current-token RCM artifacts
%   icemodel.verification.setup.reportPromiceCoverage           - Print requested-vs-actual source coverage
%   icemodel.verification.setup.researchSiteCatalog             - Return the catchall research source-site catalog
%   icemodel.verification.setup.resolveFetchCacheDir            - Apply the shared fetch-cache defaulting rule
%   icemodel.verification.setup.resolveLegWindows               - Decouple each gridded RCM leg's window from the met window
%   icemodel.verification.setup.resolveStagingRoots             - Resolve paired eval/input staging roots for importers
%   icemodel.verification.setup.retmipCaseCatalog               - Return the RetMIP protocol-case source catalog
%   icemodel.verification.setup.retmipOutputInventory           - Return variables in a RetMIP model-output NetCDF
%   icemodel.verification.setup.reuseDatasetFamilyCases         - Load staged cases for forcing-only attachment
%   icemodel.verification.setup.runDatasetFamilyDryRun          - Build a dry-run manifest through shared orchestration
%   icemodel.verification.setup.runDatasetFamilyImport          - Persist native state and optional requested RCMs
%   icemodel.verification.setup.selectSiteCatalogEntries        - Select known ids from a source-site catalog
%   icemodel.verification.setup.stageDatasetFamilyCases         - Stage requested cases with shared skip handling
%   icemodel.verification.setup.stageDatasetRcmForcing          - Delegate one RCM source at a time to stageRcmForcing
%   icemodel.verification.setup.stageMarDensityProfiles         - Add optional MAR RO1 profiles to SUMup cases
%   icemodel.verification.setup.stageRcmForcing                 - Stage RCM forcing + Data for a list of points
%   icemodel.verification.setup.stampArtifactMetadata           - Add variable metadata to staged artifact tables
%   icemodel.verification.setup.stateCaseEntry                  - Refresh one staged state record into a manifest case entry
%   icemodel.verification.setup.sumupCacheDir                   - Resolve the canonical SUMup verification source cache
%   icemodel.verification.setup.sumupComparisonVariables        - List nonempty SUMup observation groups
%   icemodel.verification.setup.textSha256                      - Return the lowercase hex SHA-256 of a text value's UTF-8 bytes
%   icemodel.verification.setup.validateEvalTarget              - Validate a case-manifest eval_target value
%   icemodel.verification.setup.validatePermafrostZone          - Validate a case-manifest permafrost_zone value
%   icemodel.verification.setup.validateSurfaceZone             - Validate a case-manifest surface_zone value
%   icemodel.verification.setup.writeFamilyManifestMerge        - Merge new case entries into a family manifest
%   icemodel.verification.setup.writeJson                       - Write pretty-printed JSON as UTF-8 with one trailing newline
%   icemodel.verification.setup.writeManifest                   - Write one verification family manifest JSON file
%   icemodel.verification.setup.writePromiceAblationReadiness   - Write the PROMICE ablation readiness ledger
%   icemodel.verification.setup.writePromiceSnowModelReadyYears - Write the annual PROMICE snow handoff
%
%   updatecontents.m generated this file on 14 Sep 2026 at 20:14:02.
