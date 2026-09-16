% +VERIFICATION
%
%   Contents file for +VERIFICATION and its subfolders.
%
%   +VERIFICATION
%   icemodel.verification.ablationPerformanceMetrics                  - Score every modeled ablation diagnostic
%   icemodel.verification.auditArtifacts                              - Read-only QA/QC for manifest-referenced artifacts
%   icemodel.verification.candidateFromIcemodelOutput                 - Convert icemodel outputs for verification
%   icemodel.verification.compareAblation                             - Compare PROMICE lowering with modeled solid-ice loss
%   icemodel.verification.comparecase                                 - Compare one staged verification target against a candidate
%   icemodel.verification.comparisonCompatibility                     - Derive staged verification comparison pairs
%   icemodel.verification.listcases                                   - Enumerate staged verification cases from family manifests
%   icemodel.verification.loadmanifest                                - Return one resolved verification case manifest
%   MAR_DENSITY_PROFILES.md
%   icemodel.verification.matchObservations                           - Match interval SMB and dated subsurface profiles
%   icemodel.verification.observationRateOutliers                     - Flag site-years whose observed ablation rate is far
%   icemodel.verification.plotcase                                    - Plot staged verification data without requiring model output
%   icemodel.verification.plotFirnArtifacts                           - Plot staged firn-family artifacts for visual QA
%   icemodel.verification.plotscatter                                 - Plot target-versus-candidate scatter panels for site cases
%   icemodel.verification.plotVerificationArtifacts                   - Visualize staged verification artifacts
%   README.md
%   README_FIXTURES.md
%   icemodel.verification.runIcemodelCandidate                        - Run icemodel and return a verification candidate
%   icemodel.verification.syntheticSnowModelRun                       - Return snow-model-like icemodel outputs
%
%   +VERIFICATION/+COLBECK
%   icemodel.verification.colbeck.analyticalSolution                  - Compute the analytical Colbeck infiltration solution
%   icemodel.verification.colbeck.caseDefinition                      - Return the canonical Colbeck 1976 verification case
%   icemodel.verification.colbeck.compareSolutions                    - Compare cached + computed Colbeck solutions side-by-side
%   icemodel.verification.colbeck.runCase                             - Build a Colbeck candidate bundle for the verification suite
%
%   +VERIFICATION/+HELPERS
%   icemodel.verification.helpers.ablationLedgerIncrements            - Per-interval ablation terms from the mass ledger
%   icemodel.verification.helpers.alignObservationSeries              - Align one observation/model series by its support
%   icemodel.verification.helpers.assertArtifactSha256                - Require current bytes to match a pinned SHA-256
%   icemodel.verification.helpers.assertRootRelativeArtifactSha256    - Verify one root-scoped artifact identity
%   icemodel.verification.helpers.classifyObservationSupport          - Apply the flag-support rules to observation rows
%   icemodel.verification.helpers.classifySnowDepth                   - Classify snow support without accepting negative depth
%   icemodel.verification.helpers.esmRuntimeMetFiles                  - Resolve an atomic ESM case's standard runtime met paths
%   icemodel.verification.helpers.esmSnowmipWaterYear                 - Return one snow water year for an ESM-SnowMIP site
%   icemodel.verification.helpers.evaluationDataRoot                  - Resolve the base evaluation-data root
%   icemodel.verification.helpers.evaluationSeason                    - Summertime display and evaluation bounds for one year
%   icemodel.verification.helpers.familyManifestFiles                 - List staged verification family manifest files
%   icemodel.verification.helpers.fieldOr                             - Return a struct field or default value
%   icemodel.verification.helpers.inputDataRoot                       - Resolve the base icemodel input-data root
%   icemodel.verification.helpers.isPhysicsFingerprint                - Return true for one well-formed physics stamp
%   icemodel.verification.helpers.loadArtifact                        - Load one named staged verification artifact from a MAT file
%   icemodel.verification.helpers.loadColocatedData                   - Assemble a timeseries bundle from staged per-source files
%   icemodel.verification.helpers.metricRowSchema                     - Return canonical comparison-metric field names/defaults
%   icemodel.verification.helpers.observationSupportFields            - Columns the observation support rules read
%   icemodel.verification.helpers.physicsFingerprint                  - Fingerprint the model's default physics configuration
%   icemodel.verification.helpers.profileGroups                       - Split profile rows by stable source identity and UTC date
%   icemodel.verification.helpers.readFamilyManifest                  - Read one verification family manifest JSON file
%   icemodel.verification.helpers.residualMetrics                     - Bias, MAE, RMSE, max error, and NSE for one paired series
%   icemodel.verification.helpers.resolveCandidateBundle              - Resolve the comparison bundle for one case
%   icemodel.verification.helpers.sampleQuantile                      - Linear-interpolated sample quantile, no extra toolboxes
%   icemodel.verification.helpers.sumupColocation                     - Flag a SUMup point as co-located with a mixed anchor
%   icemodel.verification.helpers.validateAblationModelSchema         - Check a saved cohort still supports the report
%   icemodel.verification.helpers.writeRunReport                      - Write a concise markdown report for a verification run
%
%   +VERIFICATION/+NAMELISTS
%   icemodel.verification.namelists.ablationReportChannels            - Model channels read by the ablation report
%   icemodel.verification.namelists.caseid                            - Return supported runnable snow-verification case ids
%   icemodel.verification.namelists.casetype                          - Return the supported snow-verification case types
%   icemodel.verification.namelists.completions                       - Return the supported verification namelist selector names
%   icemodel.verification.namelists.datasetfamily                     - Return the supported snow-verification dataset families
%   icemodel.verification.namelists.evaltarget                        - Return the supported case-manifest eval-target descriptors
%   icemodel.verification.namelists.firndatasetfamily                 - Return verification families used by firn staging previews
%   icemodel.verification.namelists.laughtests                        - Canonical Laugh-Tests case-id namelist
%   icemodel.verification.namelists.permafrostzone                    - Return the supported case-manifest permafrost-zone values
%   icemodel.verification.namelists.promiceAblationPolicy             - Return the fixed PROMICE ablation comparison policy
%   icemodel.verification.namelists.promiceAblationReadiness          - Return the PROMICE ablation admission policy
%   icemodel.verification.namelists.promicesite                       - Auto-discovered PROMICE station-id namelist
%   icemodel.verification.namelists.rcmMetSources                     - Verification RCM labels that currently write met files
%   icemodel.verification.namelists.rcmProductIds                     - Map RCM runtime/storage labels to explicit product ids
%   icemodel.verification.namelists.rcmsources                        - Verification RCM source labels in canonical staging order
%   icemodel.verification.namelists.snowmipsite                       - Canonical ESM-SnowMIP site-name namelist
%   icemodel.verification.namelists.surfacezone                       - Return the supported case-manifest surface-zone values
%
%   +VERIFICATION/+REPORT
%   build_snow_artifact_qa.py
%   icemodel.verification.report.buildAblationEvaluationReport        - Render saved PROMICE ablation results
%   icemodel.verification.report.buildGapFillReport                   - Build the gap-fill before/after Quarto report inputs
%   icemodel.verification.report.buildTestSuiteReport                 - Render numerical or performance suite results
%   check_snow_artifact_qa.py
%   icemodel.verification.report.configureCategoryAxis                - Label horizontal evidence rows without clipping
%   icemodel.verification.report.escapeMarkdownText                   - Show saved text literally, without markup or raw HTML
%   icemodel.verification.report.exportAndClose                       - Export one report figure and release graphics state
%   icemodel.verification.report.formatReportAxes                     - Isolate exported graphics from interactive theme defaults
%   icemodel.verification.report.formatValue                          - Format one scalar table value for Markdown
%   icemodel.verification.report.gapfillFigureStyle                   - Define the gap-fill report figure colors
%   icemodel.verification.report.generateFinalFirnPreview             - Build canonical firn QA and figure products
%   icemodel.verification.report.generateFinalSnowPreview             - Build canonical seasonal QA, figures, and readiness
%   icemodel.verification.report.markdownCode                         - Wrap saved metadata in a code span that renders literally
%   icemodel.verification.report.markdownTable                        - Convert a compact table to inert Markdown
%   icemodel.verification.report.methodFillLayers                     - Split one channel into observed/own-fill/other-fill layers
%   README.md
%   icemodel.verification.report.safeLabel                            - Collapse control characters for MATLAB graphics text
%   icemodel.verification.report.sanitizeText                         - Collapse controls and neutralize raw HTML delimiters
%   snow-artifact-qa.qmd
%
%   +VERIFICATION/+SETUP
%   icemodel.verification.setup.anchorColocation                      - Flag a point as co-located with the nearest anchor
%   icemodel.verification.setup.buildDatasetFamilyManifest            - Build and merge-write a family manifest
%   icemodel.verification.setup.buildEsmSnowmipForcing                - Convert ESM-SnowMIP NetCDF to native forcing
%   icemodel.verification.setup.buildEsmSnowmipObservations           - Convert ESM-SnowMIP obs NetCDF to verification targets
%   icemodel.verification.setup.buildFetchProductStatus               - Build ordered status rows from a product registry
%   icemodel.verification.setup.buildLaughTestsArtifacts              - Build one Laugh-Tests evaluation/reference bundle
%   icemodel.verification.setup.buildSumupObservations                - Convert SUMup firn records to verification targets
%   icemodel.verification.setup.bytesSha256                           - Return the lowercase hex SHA-256 of a byte vector
%   icemodel.verification.setup.caseManifestFieldNames                - Return canonical case-manifest fields
%   icemodel.verification.setup.colocationSourceLists                 - Derive manifest source lists from colocation legs
%   icemodel.verification.setup.datasetFamilyStagingPaths             - Build the shared dataset-family output paths
%   icemodel.verification.setup.deduplicateSumupRecords               - Keep one row per SUMup scientific identity
%   icemodel.verification.setup.emptyFetchProductStatusRow            - Return the literal shared fetch-row prototype
%   icemodel.verification.setup.ensureUtc                             - Coerce a date/datetime input to a UTC-tagged datetime
%   icemodel.verification.setup.esmSnowmipSiteCatalog                 - Return the ESM-SnowMIP source-site catalog
%   icemodel.verification.setup.familyManifestFieldNames              - Return canonical family-manifest fields
%   icemodel.verification.setup.fetchEsmSnowmip                       - Locate or verify the ESM-SnowMIP source NetCDF files
%   icemodel.verification.setup.fetchFixtures                         - Transactionally provision or verify release data
%   icemodel.verification.setup.fetchGcnet                            - Locate or verify local Vandecrux/GC-Net source caches
%   icemodel.verification.setup.fetchImau                             - Locate or verify local IMAU PANGAEA source caches
%   icemodel.verification.setup.fetchKtransect                        - Locate or verify the local K-transect PANGAEA source cache
%   icemodel.verification.setup.fetchLaughTests                       - Locate or verify the Laugh-Tests source checkout
%   icemodel.verification.setup.fetchMissingStatus                    - Convert missing source patterns to fetch status rows
%   icemodel.verification.setup.fetchProductFiles                     - Return files matching product-cache patterns
%   icemodel.verification.setup.fetchProductNames                     - Return ordered product selectors from a fetch registry
%   icemodel.verification.setup.fetchProductStatusRow                 - Build the standard fetch product status record
%   icemodel.verification.setup.fetchPromice                          - Locate or verify local PROMICE pypromice L3 source caches
%   icemodel.verification.setup.fetchRetmip                           - Locate or verify local RetMIP source caches
%   icemodel.verification.setup.fetchSumup                            - Locate or verify the SUMup firn source files
%   icemodel.verification.setup.fileSha256                            - Return the lowercase hex SHA-256 of a file's bytes
%   icemodel.verification.setup.finishFetchStatus                     - Apply the shared fetch status/strict/silent contract
%   icemodel.verification.setup.firnCaseManifestFieldNames            - Return canonical firn case-manifest fields
%   icemodel.verification.setup.fixtureCallerSymlink                  - Find the first caller-controlled link in a path
%   icemodel.verification.setup.fixtureDataRoot                       - Return the data root registered for a release
%   icemodel.verification.setup.fixtureFetchCommand                   - Build a copy-paste release-data repair command
%   icemodel.verification.setup.fixtureFileList                       - Return manifest paths for selected release capabilities
%   icemodel.verification.setup.fixtureRelativePosix                  - Return one root-relative path with POSIX separators
%   icemodel.verification.setup.formatManifestTime                    - Serialize manifest timestamps with explicit clock time
%   icemodel.verification.setup.gcnetInventory                        - Index Vandecrux/GC-Net products without loading arrays
%   icemodel.verification.setup.gcnetProductNames                     - Return canonical Vandecrux/GC-Net product selectors
%   icemodel.verification.setup.gcnetProductSpec                      - Return Vandecrux/GC-Net DOI metadata and file rules
%   icemodel.verification.setup.hasObservationRecords                 - True when any observation sub-bundle carries rows
%   icemodel.verification.setup.imauSiteCatalog                       - Return the IMAU hourly-AWS source-site catalog
%   icemodel.verification.setup.importEsmSnowmip                      - Stage ESM-SnowMIP site fixtures for the verification suite
%   icemodel.verification.setup.importImau                            - Stage the IMAU hourly AWS verification family
%   icemodel.verification.setup.importKtransect                       - Stage the K-transect annual AWS verification family
%   icemodel.verification.setup.importLaughTests                      - Stage selected Laugh-Tests synthetic snow benchmarks
%   icemodel.verification.setup.importPromiceSites                    - Stage PROMICE-anchored firn-evaluation cases
%   icemodel.verification.setup.importResearchSites                   - Stage generic research-site firn targets
%   icemodel.verification.setup.importRetmip                          - Stage the RetMIP protocol verification family
%   icemodel.verification.setup.importSumup                           - Stage co-located SUMup firn evaluation cases
%   icemodel.verification.setup.ktransectAliasCrosswalk               - Return the K-transect station alias hypothesis table
%   icemodel.verification.setup.ktransectSiteCatalog                  - Return the K-transect AWS source-site catalog
%   icemodel.verification.setup.loadPriorDatasetFamilyCases           - Read cases needed by an additive native refresh
%   icemodel.verification.setup.makeCaseManifestEntry                 - Build one case manifest entry from canonical fields
%   icemodel.verification.setup.makeFamilyManifest                    - Build one verification family manifest struct
%   icemodel.verification.setup.makeFirnCaseManifestEntry             - Build one firn case manifest entry
%   icemodel.verification.setup.manifestWindow                        - Serialize one start/end pair for a JSON manifest
%   icemodel.verification.setup.mergeColocation                       - Copy every field from ADD onto a colocation struct
%   icemodel.verification.setup.metadataStruct                        - Build a metadata struct from a 2-column cell array
%   icemodel.verification.setup.metArtifactReadiness                  - Diagnose the exact saved scalar-window met artifact
%   icemodel.verification.setup.metForcingReady                       - Test unfilled readiness and inventory complete windows
%   icemodel.verification.setup.mixedAnchorCatalog                    - Read staged firn/research anchors from manifests
%   icemodel.verification.setup.mustBeKnownFetchProducts              - Reject selectors absent from a fetch registry
%   icemodel.verification.setup.normalizeForcingSources               - Normalize one public forcing-source selection
%   icemodel.verification.setup.packFixtures                          - Pack selected release-data capabilities into separate archives
%   icemodel.verification.setup.periodBounds                          - Parse a manifest period to UTC datetimes
%   icemodel.verification.setup.preferPrimary                         - Coalesce two same-shape series, keeping primary where finite
%   icemodel.verification.setup.prepareCaseRoot                       - Create one case folder and select artifacts needing writes
%   icemodel.verification.setup.prepareReplacementCaseEntry           - Retain observations but clear prior runtime legs
%   icemodel.verification.setup.preservePriorNativeLeg                - Merge prior native artifacts into a fresh case leg
%   icemodel.verification.setup.preserveRcmLegs                       - Keep compatible staged RCM legs after a failed refresh
%   icemodel.verification.setup.previewFirnStaging                    - Stage short build_forcing=true firn QA previews
%   icemodel.verification.setup.printFetchProductBanner               - Print shared manual-cache retrieval instructions
%   icemodel.verification.setup.priorCaseById                         - Return one prior manifest case by canonical case id
%   icemodel.verification.setup.promiceSiteCatalog                    - Return the PROMICE source-site catalog
%   promote_snow_verification_artifacts.py
%   icemodel.verification.setup.rcmArtifactOutputDirs                 - Resolve shared default RCM artifact output roots
%   icemodel.verification.setup.rcmSourceCoverage                     - Probe on-disk year coverage of each forcing source
%   icemodel.verification.setup.rcmStorageAlias                       - Return the collision-safe RCM artifact identity for a case
%   icemodel.verification.setup.readBestSnowDepth                     - Site-aware snow-depth selector for ESM-SnowMIP obs
%   icemodel.verification.setup.readNetcdfTime                        - Read a NetCDF time coordinate as UTC datetime
%   icemodel.verification.setup.readNetcdfVariable                    - Read one ESM-SnowMIP-style NetCDF variable with NaN-fill
%   icemodel.verification.setup.readRcmArtifactMetadata               - Read saved RCM provenance without payload arrays
%   icemodel.verification.setup.readRetmipProfileTable                - Read a RetMIP initial profile table
%   icemodel.verification.setup.readRetmipProtocolTable               - Read a RetMIP tab-delimited protocol time series
%   icemodel.verification.setup.refreshManifestSourceLists            - Recompute source lists without rebuilding data
%   icemodel.verification.setup.refreshPromiceMetIdentities           - Pin existing native met bytes in the manifest
%   icemodel.verification.setup.regexpOnce                            - Return one stripped regexp token or an empty string
%   icemodel.verification.setup.releaseManifestFile                   - Return the release-data manifest path for a version
%   icemodel.verification.setup.relpaths                              - Reduce absolute staged paths to base-relative names for JSON
%   icemodel.verification.setup.repairMetTimeSupport                  - Repair legacy linear 15-minute met artifacts
%   icemodel.verification.setup.repairRcmArtifactMetadata             - Classify and repair current-token RCM artifacts
%   icemodel.verification.setup.reportPromiceCoverage                 - Print requested-vs-actual source coverage
%   icemodel.verification.setup.researchSiteCatalog                   - Return the catchall research source-site catalog
%   icemodel.verification.setup.resolveFetchCacheDir                  - Apply the shared fetch-cache defaulting rule
%   icemodel.verification.setup.resolveLegWindows                     - Decouple each gridded RCM leg's window from the met window
%   icemodel.verification.setup.resolveStagingRoots                   - Resolve paired eval/input staging roots for importers
%   icemodel.verification.setup.retmipCaseCatalog                     - Return the RetMIP protocol-case source catalog
%   icemodel.verification.setup.retmipOutputInventory                 - Return variables in a RetMIP model-output NetCDF
%   icemodel.verification.setup.reuseDatasetFamilyCases               - Load staged cases for forcing-only attachment
%   icemodel.verification.setup.runDatasetFamilyDryRun                - Build a dry-run manifest through shared orchestration
%   icemodel.verification.setup.runDatasetFamilyImport                - Persist native state and optional requested RCMs
%   icemodel.verification.setup.selectSiteCatalogEntries              - Select known ids from a source-site catalog
%   icemodel.verification.setup.stageDatasetFamilyCases               - Stage requested cases with shared skip handling
%   icemodel.verification.setup.stageDatasetRcmForcing                - Delegate one RCM source at a time to stageRcmForcing
%   icemodel.verification.setup.stageMarDensityProfiles               - Add optional MAR RO1 profiles to SUMup cases
%   icemodel.verification.setup.stageRcmForcing                       - Stage RCM forcing + Data for a list of points
%   icemodel.verification.setup.stampArtifactMetadata                 - Add variable metadata to staged artifact tables
%   icemodel.verification.setup.stateCaseEntry                        - Refresh one staged state record into a manifest case entry
%   icemodel.verification.setup.sumupCacheDir                         - Resolve the canonical SUMup verification source cache
%   icemodel.verification.setup.sumupComparisonVariables              - List nonempty SUMup observation groups
%   icemodel.verification.setup.textSha256                            - Return the lowercase hex SHA-256 of a text value's UTF-8 bytes
%   icemodel.verification.setup.validateEvalTarget                    - Validate a case-manifest eval_target value
%   icemodel.verification.setup.validatePermafrostZone                - Validate a case-manifest permafrost_zone value
%   icemodel.verification.setup.validateSurfaceZone                   - Validate a case-manifest surface_zone value
%   icemodel.verification.setup.writeFamilyManifestMerge              - Merge new case entries into a family manifest
%   icemodel.verification.setup.writeJson                             - Write pretty-printed JSON as UTF-8 with one trailing newline
%   icemodel.verification.setup.writeManifest                         - Write one verification family manifest JSON file
%   icemodel.verification.setup.writePromiceAblationReadiness         - Write the PROMICE ablation readiness ledger
%   icemodel.verification.setup.writePromiceSnowModelReadyYears       - Write the annual PROMICE snow handoff
%
%   +VERIFICATION/+VALIDATORS
%   icemodel.verification.validators.mustBeCaseIdSubset               - Validate a subset of canonical snow-verification case ids
%   icemodel.verification.validators.mustBeDatasetFamilyFilter        - Validate one optional dataset-family selector
%   icemodel.verification.validators.mustBeDatasetFamilySelection     - Validate dataset-family selectors plus "all"
%   icemodel.verification.validators.mustBeFirnDatasetFamilySelection - Validate firn-family selectors plus "all"
%   icemodel.verification.validators.mustBeGcnetProductSelection      - Validate Vandecrux/GC-Net product selectors
%   icemodel.verification.validators.mustBeLaughTestCase              - Validate cases against the Laugh-Tests namelist
%   icemodel.verification.validators.mustBeRcmSourceSelection         - Validate verification forcing-source selectors
%   icemodel.verification.validators.mustBeSnowmipSite                - Validate sitename against the canonical ESM-SnowMIP
%
%   updatecontents.m generated this file on 14 Sep 2026 at 20:14:01.
