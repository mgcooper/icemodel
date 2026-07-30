# PROMICE Gap-Fill Reconstruction Policy

Owning DesignSpec: `.agents/plans/design-specs/2026-07-23-promice-gap-filling-and-ktransect.md`
Owning ExecPlan: `.agents/plans/execution/2026-07-23-icemodel-g1n-promice-gap-filling-and-ktransect.md`
Status: HARDENED CONTRACT (adversarial hardening pass merged 2026-07-26
with user sign-off; supersedes the Gate-A draft-with-inline-revisions
form). Sections A and B are the review-loop oracle; Section C parameters
are never review-blocking; Section D is the authoritative decision log.
Runtime defaults live in `icemodel.forcing.reconstruct.setopts` and must
match Section C by mechanical diff.

Authority tags: `[ratified]` Gate A or explicit later approval;
`[user-directed]` explicit user order; `[agent-proposed]` autonomy-era
rule ratified via D-11; `[codex-amended]` review-cycle contract
correction; `[parameter]` tunable via `setopts`.

---

## A. Normative requirements

A1 `[ratified]` **Native immutability.** Native pypromice L3 staged values
are never overwritten; the filled product is the separately named
`met_<site>_promice_filled` artifact with per-sample uint8 provenance.
Byte-pinned identity: a staged native target or PROMICE donor must match
the importer manifest's size+SHA-256 before it trains reconstruction.
Pinning is a fingerprint ledger, never a staging burden: when artifacts
predate a pinning rule, the approved path is refreshing the manifest —
re-hashing the files already on disk, a seconds-scale metadata operation
— and verification blocks only a genuine on-disk mismatch, exactly the
case worth stopping for. Identity checks are never bypassed and never
require re-importing sources (D-22).

A2 `[user-directed]` **Canonical runnable product.** `promice_filled` is
the canonical runnable PROMICE forcing; native `promice` is retained
unmodified as the provenance/QC source.

A3 `[user-directed]` **Runnability over optional-quality gates.** No
optional-quality metadata may block a run, a donor, or readiness. Boom
height resolves measured → interpolated (including across
station-composition handovers) → nominal 2.6 m constant, each step
recorded with provenance and a warning. `icemodel-1ps.16`'s
maintenance-visit registry later improves the interpolation step; it is
not a gate. K-transect height provenance annotates/downgrades donor
records; it never blocks staging or donor use. The `boom_height >
z0_bulk` check validates measured values; failures fall down the chain
rather than erroring.

A4 `[user-directed]` **Runtime readiness gate = requested-window
coverage.** `loadmet` gates `promice_filled` on complete required-channel
coverage of every requested timestep between startdate and enddate (water
years and arbitrary windows are first-class); calendar-year ledger
verdicts are bookkeeping and reporting, never the runtime gate. Identity
checks retained: producer-manifest hash, provenance registry + typed
channels, product/site identity, readiness snapshot binding
`[codex-amended]`.

A5 `[user-directed]` **Two readiness verdicts.** The ledger records
`ready_icemodel` (requires tair, rh, wspd, psfc, swd, lwd, albedo) and
`ready_snowmodel` (additionally requires total ppt OR snowf) per
station-year. `swu` is derived, never required; `rainf` is never required
(rain inclusion is a model option). Both verdicts require ZERO residual
missing in their set within the product span (A6) — the model consumes
every sample; an absent required channel counts as wholly missing. Stored
verdicts are absolute; the A6 policy view never overrides them.

A6 `[user-directed]` **The proxy window BOUNDS the product.** The
per-site acceptance window — the union span of staged MAR/MERRA proxy met
actually loaded (`acceptanceWindow`; validated filename↔timetable
agreement, continuous support, no disjoint unions `[codex-amended]`) — is
the filled product's time span: native samples outside it are not
reconstructed and are not part of the filled artifact (they remain
untouched in the native product). No station is ever "not ready" because
of an era no proxy covers; donor availability never extends the span —
the MAR/MERRA period defines the gap-filling period. The window derives
from staged files at run time (nothing stores it), so staging more proxy
met automatically widens the product on regeneration. Report claims pin
to the producer-recorded window and options; the live view serves
interactive queries. If no validated 15-minute MAR/MERRA proxy is staged,
the acceptance-window union is empty, no fillable product span exists, and
production refuses with `noProxyWindow`; an empty union never authorizes
full-native-record reconstruction from donors or climatology
`[codex-amended]`.

A7 `[ratified]` **Provenance honesty.** Every filled sample carries its
method code; every contiguous filled segment carries an exact audit row
joined by `context_id` to one persisted plan record; provisional
`unfilled` rows are reconciled after later tiers `[codex-amended]`.
Builder-constructed values never masquerade as observations:
winter-albedo stamps and annual interpolation fills are masked on input
and refilled honestly; builder darkness zeros (swd, swu) are replayed
against the raw source and owned by the darkness/derivation methods;
estimated-LWD-flagged artifacts are masked before planning; raw-fallback
and clamped shortwave carry codes 13/14 `[codex-amended]`. A finite
albedo sample whose raw source link is unverifiable is treated as
not-observed (fail closed) — by A3 this affects provenance labels only;
masked samples remain fillable by the tiers.

A8 `[user-directed]` **Donor eligibility.** Donors are other stations'
MEASUREMENT products: reconstructed, gap-filled, or estimated values are
never donors (GC-Net donor samples require per-channel origin-observed
flags; unflagged channels are ineligible). Values derived from a
station's own coincident measurements are measurements — K-transect
measured-ratio albedo (swu/swd) is donor-eligible. The target is never
its own donor, through every alias spelling `[codex-amended]`.

A9 `[ratified]` **Methods compete; nothing is grandfathered.** Admission
is exclusively via the held-out harness (B7/B8). In-sample overlap
statistics are never acceptance evidence. A record with no disjoint whole
evaluation year cannot admit a statistical method — but whole-record
proxy replacement (A11) is guaranteed available, so single-year stations
remain fillable `[user-directed]`.

A10 `[user-directed]` **Precipitation.** Native `rainf` (corrected
tipping-bucket observations) is preserved wherever finite. Total `ppt`
adopts from MAR, else MERRA-2; where MAR provides its own
`rainf`/`snowf`, those ship in the product as MAR's split. NO partitioning
is applied at reconstruction time. At runtime a model option selects the
phase source: `source` (the product's rainf/snowf, e.g. MAR's
energy-balance split; default) or `threshold` (partition the product's
ppt with the transition-temperature function). The solver currently
consumes zero rain until the advective-rain physics is finished (D-0b);
the product still carries real rainf. Components must be nonnegative and
sum to the total; stations are never total-precip sources; precipitation
never enters the statistically planned method set `[codex-amended]`.

A11 `[user-directed]` **Last resort (distinct from competition).**
COMPETITION is held-out admission of planned methods per stratum; the
LAST RESORT is a fixed, parameterized per-channel preference ORDER
applied to samples no admitted method covered (albedo: MODIS → MAR →
MERRA-2; every other channel: MAR → MERRA-2). One source per contiguous
outage (coupled-channel consistency); overlap-calibrated where overlap
exists; bounds-checked; proxy provenance codes; audited as
`proxy:<source>:last_resort`. Zero-overlap sources are excluded from
competition but ARE admitted as last resort with a low-confidence flag
where nothing else exists. `swu` is excluded (derived, B10). Samples with
no in-bounds fallback value proceed to the post-final short-gap pass (B3);
only samples that still have no admissible bridge remain missing and
ledgered.

A12 `[ratified]` **Proxy scope.** MAR (default) and MERRA-2 (alternative
and post-MAR-era precip source) are the only proxy met families; order is
one `setopts` argument; RACMO stays subsurface-only; a future source is a
deliberate code change at one marked extension point. MODIS is an
albedo-only observational source (A11, B12), not a met proxy family.
CARRA noted, out of scope.

A13 `[user-directed]` **Native-support advisory.** A station-year with
< 30% (parameter) in-record native core-channel (tair, rh, psfc)
coverage carries a low-confidence advisory and its coverage in the
ledger. The advisory never replaces or demotes the forcing-readiness
verdict: every fillable in-window required sample must still be valid,
and a complete product is still ready. Coverage is judged within the A6
product span.

A14 `[user-directed]` **Report contract.** Fixed structure (Executive
Summary, Background, Methods, Results, Summary, Appendices); shared
plotting layer only; human-readable durations; inputs verified by
size+SHA-256 and site/product identity; report gate numbers and windows
come from producer-pinned records, never live defaults. Per-site appendix
layout (D-23): one OVERVIEW figure per station — a single column of 8
tiles (one per required channel) spanning the full period covering that
station's filled gaps — plus up to ~6 (parameter) representative
gap-detail figures in the windowed style (context pad per side = filled
period, floored and capped by parameters), CURATED FOR METHOD DIVERSITY
so the detail set doubles as a method comparison (the per-method sketch's
purpose). Never one figure per filled segment. The main report embeds
the scientific Results figures and station overviews only; every
method/gap appendix figure and its per-station diagnostic tables live in
a companion detail report. Both reports render PDF and HTML from their
saved QMDs. Ledger reconciliation covers exactly the rendered set.
Results must explain every remaining missing forcing run and every
non-ready station-year in plain language from producer-recorded evidence.
The report must verify those run rows against the shipped product, state
explicitly when none remain, and independently prove why any staged
station has no product; file absence alone may not invent a refusal
reason (D-38).

A15 `[ratified]` **Physical validity.** Shared per-channel scalar bounds
are hard limits at every FILL tier including last resort; readiness
verdicts grade completeness plus the scalar bounds, while RELATIONAL
exceedances (swd vs TOA, swu vs swd) are reported as ledger DIAGNOSTIC
columns and never flip a verdict (D-28: twilight diffuse light is real
incident energy the TOA model cannot represent, and native sensor noise
must not brand usable years not-ready). Relational rules remain hard
gates for fill CANDIDATES except the D-32 post-final SWD flux-linear
bridge: its interval-start anchors can legitimately exceed an
instantaneous horizontal-TOA ceiling, so that bridge enforces the scalar
bound and retains the solar relation as a diagnostic rather than clipping
a continuous local interpolation into an artificial low seam. The RULE is
normative; the bound VALUES are
Section-C parameters (D-11, D-25). Current values: tair [193, 300] K;
rh [5, 100] % with calibrated candidates clamped into bounds (D-27);
wspd [0, 60] m s⁻¹; psfc [60000, 108000] Pa; swd
[0, max(1.05×TOA, 5)] W m⁻² for candidates (the 5 W m⁻² term is a
minimum CEILING in darkness admitting thermal-offset noise, not a floor
on data — D-28); swu [0, swd]; albedo [0.05, 0.98]; lwd [40, 470]
W m⁻² (floor lowered from 100 and ceiling raised from 400 on 2026-07-27
evidence — MAR produces genuine 69–90 W m⁻² cold clear-sky values, a
blackbody at the tair floor emits ~79 W m⁻², and warm-fjord stations
observe genuine 406–451 W m⁻² under mild overcast; D-25/D-26); ppt ≥ 0;
measured boom height strictly above `z0_bulk` (else fall back per A3).

A16 `[user-directed]` **Channel naming (D-24).** Source-schema names die
at the read boundary: inside every icemodel-owned artifact (met,
userdata, filled products, provenance, reports) the icemodel canonical
names apply — upward shortwave is `swu`. `usr` exists only while reading
pypromice sources (ingest alias). This is the general rule, not a
product-specific exception; `swu` itself ships in the filled product
purely as the derived QC/report channel (B10) and is never a required
forcing input (A5).

A17 `[user-directed]` **Corroborated native sensor failure.** A native
run that satisfies the multi-channel burial/rime screen (B14) remains
byte-identical in the staged artifact but is treated as missing in the
reconstruction working copy. Its implicated channels cannot train a
method, define a seam scale, or serve as donor truth; they remain
fillable by the normal tier ladder. The exact finding is persisted with
the station plan and represented by a fixed scientific-interpretation
category in the report. A single flat channel is insufficient evidence.

---

## B. Method contract

B1 **Tier order** `[ratified]`: darkness pre-pass (swd/swu) → tier-1
bounded interpolation → donor transfer → proxy calibration → empirical
estimators → climatology → last resort (A11) → post-final bounded
interpolation (B3) → one-posting twilight climatology (D-44) →
winter-albedo bridge (B13). *Plain language: try the
most information-rich source first, then reconnect any short seam slivers
that the completed ladder leaves; every fill is validated before it may
ship.*

B2 **Darkness rule** `[ratified]`: when the sun is below the
meaningful-sun threshold, missing shortwave is a KNOWN zero, filled first
with its own provenance code; builder-made night zeros are replayed
against the raw source so this method owns them (A7). *Night sunlight is
knowledge, not estimation — and it splits multi-day outages into daylight
pieces the validation actually graded.*

B3 **Tier-1 interpolation** `[ratified]`: straight-line bridging of
interior gaps up to the PER-CHANNEL cap (6 h except SWD 9 h under D-39,
RH 9 h under D-42/D-50, and albedo 30 h under D-41), state variables plus
short albedo slivers;
shortwave interpolates in clear-sky-index space (the
fraction of possible sunlight), never raw across daylight; a bridge never
spans two seasons or crosses darkness, except an SWD gap within its existing
9 h cap may borrow an adjacent finite anchor across, or itself straddle, the
calendar boundary under D-43/D-49; RH may cross the calendar boundary only
within its existing 9 h cap under D-48/D-50 `[codex-amended]`. Loosening any
channel's cap requires held-out evidence at the target duration against
the persistence baseline; `psfc` (smooth, synoptic) is the flagged
candidate for a 12–24 h extension. The same bounded operation runs once
more after last resort, using the final finite neighbors as anchors, to
close residual interior seam slivers that were not bridgeable before the
other tiers ran. For SWD this pass tries CSI first, then uses a
flux-linear bridge only for still-missing slivers; KANL heldouts over
137,708 observed samples show flux-linear RMSE 45.9 W m⁻² versus
persistence 103.9, and shoulder-hour RMSE 28.9 versus 71.6 where CSI has
no coverage. It uses the same caps, season/darkness boundaries,
provenance code, and exact audit; scalar bounds remain hard and the solar
relation follows the A15 diagnostic exception `[user-directed]`. *Short
gaps get connected smoothly; the sun's daily arc is preserved by CSI
where defined and by the observed local ramp at its dawn/dusk boundary.*
For guarded hourly or half-hourly sources, the same residual closer runs
once more after mean-preserving restoration to the delivered 15-minute
axis (D-47); any newly finite SWD or albedo operand is followed immediately
by the ordered SWU derivation and final audit reconciliation.

B4 **Donor transfer** `[ratified]`: per-season regression from a nearby
station on overlapping data, with the monotone spline adjustment as a
validated hyperparameter, bounded lag search, elevation adjustment (lapse
rate for tair, barometric for psfc) fitted on overlap with recorded
fallbacks, extrapolation refused beyond the fitted range plus margin,
geometry gate, elevation-aware bracketing (e.g. KAN_M from S6+S9) as a
candidate. SWD donors transfer in CSI space `[codex-amended]`. K-transect
geometry authenticates against the byte-pinned artifact. *A neighbor's
weather is translated using the periods both stations measured, and never
trusted outside what that translation saw.*

B5 **Proxy calibration** `[ratified]`: MAR/MERRA values bias-corrected on
the station overlap (seasonal where supported; additive for state
channels; multiplicative for shortwave and nonnegative wind speed, with
meaningful-sun screening and a positive proxy denominator for shortwave
`[codex-amended]`); identity corrections are admissible only through the
gates or as flagged last resort (A11).

B6 **Seam blending** `[ratified]`: where a fill meets real data, a step
larger than the seasonal jump limit is smoothed away over the blend
window — excess only, skilled seams untouched; the per-run limit is
floored at the linear bridge's implied step (physically real jumps — dawn
beside darkness, fronts — are never grounds for refusal; a single-sample
run with both anchors binding takes the bridge midpoint); disjoint
fragments taper independently with zero endpoints inside the gap; the
scale is frozen from native observations before any filling
`[codex-amended]`. After full composition, SWD transitions are screened
against observed steps in the same season and solar-elevation band, with
the native-anchor linear-bridge step as a floor. An outlier may mask and
flux-linearly reconnect only the reconstructed posting beside the
boundary; the parameterized pass count is bounded, native samples are
immutable, and no observed reference forces a review verdict rather than
a silent pass (D-33). Refusal is reserved for bounds violations and
missing estimates. *Fills are stitched in without fake jumps, and the
stitching threshold cannot be corrupted by the fills themselves.*

B7 **Held-out validation** `[ratified]`: synthetic gaps drawn from each
station's real gap-length census per exact stratum (season ×
duration bucket; strata and buckets from `bucketEdges`), inserted into
observed spans (daylight-bounded for shortwave); whole-year
selection/evaluation split (~70/30, persisted, schema-verified, disjoint,
joint-core-supported years only `[codex-amended]`); selection evidence is
union-masked from evaluation everywhere — fits, climatology pools,
persistence anchors, seam scales `[codex-amended]`. *Methods audition on
artificial gaps where the truth is known, and the exam questions are
never in the study set.*

B8 **Admission gates** `[ratified]`: per stratum — (a) |bias| within
instrument-class caps; (b) RMSE improves on the baseline (persistence for
the ≤6 h bucket, station climatology otherwise) by the required fraction
on IDENTICAL common support meeting a minimum-support floor
`[codex-amended]`; (c) zero bound violations; (d) coverage above the
usefulness floor; (e) finite variability ratio within the band, with
unmeasurable within-gap spread counting as denial, not exemption
`[codex-amended]`; (f) composition never applies a method to a real gap
longer than its longest validated gap in that bucket. Admitted methods
rank by common-support fractional improvement, RMSE tie-break. Metrics
follow the multi-component practice of Knoben et al. (2019),
doi:10.5194/hess-23-4323-2019, without hiding components in one score.
*A method must beat the dumb baseline on the same questions, must not
smooth reality away, and never works beyond where it was tested.*

B9 **Climatology** `[ratified]`: day-of-year × exact posting-time lookup
on a 365-day calendar (Feb 29 pools into Feb 28, parameter note);
admission baseline and competing method.

B10 **swu derivation** `[ratified]`: swu = albedo × swd, derived after
all operand tiers, provenance code 12, audited as the algebraic relation;
native finite swu preserved; never independently admitted, adopted, or
required.

B11 **Cross-season outages** `[codex-amended → ratified]`: production
outages split at season boundaries into independently admitted segments;
completed neighbors never become synthetic seam anchors.

B12 **MODIS albedo (activated, D-15).** Full GEUS C6 daily albedo
2000–2019 is available at `/Volumes/S03/DATA/greenland/geus/albedo/gris`
(earlier "2012-only" notes described the repo-local fixture cache).
Pipeline: inventory currently staged coverage (legacy met `modis`
columns, modis userdata), stage per-site MODIS userdata at native daily
cadence, attach through ONE daily→met-cadence interpolation helper at
staging/reconstruction time (linear between daily values; runtime
attachment rejected — artifacts stay self-contained), provenance code 7.
Competes under the gates AND holds first position in albedo's last-resort
order (A11). *Daily satellite albedo is adequate sub-daily forcing here
because albedo varies with snowfall and synoptic events, not the hour of
day.*

B13 **Winter-albedo bridge (D-15a).** The legacy 0.8 constant was a
run-enabling heuristic; for snow/firn modeling the fallback becomes a
bounded seasonal bridge: linear interpolation between the gap-edge
observations, floored at the dry-snow value (`parameterLookup`), applied
only to masked winter albedo that nothing better fills (methods, MODIS,
RCM last resort all rank ahead). Stamped with its own audit trail.
*Winter sunlight is near zero so the forcing impact is small; the bridge
fixes the unphysical step edges where the sun returns.*

B14 **Native flat-run screen (D-37).** Before target training or donor
assembly, screen daily native observations for a sustained near-flat
temperature run corroborated by near-zero shortwave despite substantial
TOA irradiance and longwave near the blackbody temperature. Mask only the
channels named by the corroborating evidence in the working copy; never
rewrite the staged artifact. Persist start, end, duration, implicated
channels, and supporting diagnostics. *A buried or rimed station is
source-QC evidence, not trustworthy weather and not permission to erase
the source record.*

---

## C. Parameters (change via `reconstruct.setopts`; never review-blocking)

Per-channel tier-1 caps (all 6 h default; psfc flagged extension
candidate, B3); jump factor 3× seasonal median step; blend window 6 h
(seam lands at half the limit); meaningful-sun TOA threshold 10 W m⁻²;
donor geometry 60 km / 600 m; elevation-adjust threshold 100 m; lapse
rate −0.0060 K m⁻¹; pressure scale temperature 255 K; min overlap 8760
elapsed donor-support hours `[codex-amended]`; lag search ±18 h, min lag
gain 0.02; extrapolation fraction 0.10; RMSE improvement 0.10; coverage
floor 0.10 (revised from 0.9 on KAN evidence, D-log); variability band
[0.67, 1.50]; seasonal-calibration support 300; climatology window ±7 d,
min support 5; planner draws 12/stratum; synthetic draws 25, context
24 h; knot candidates [0 6]; max donors 3/stratum; selection fraction
0.7; triage floor 0.30; instrument-class bias caps (tair 0.5 K; rh 5%;
wspd 1 m s⁻¹ or 10% via stratum typical speed; psfc 100 Pa; swd/swu 20
W m⁻² daily-mean; lwd 15 W m⁻²; albedo 0.05); winter-albedo bridge floor
0.80 + months (`parameterLookup`); rain/snow transition temperature = Tf;
per-channel last-resort order (albedo [modis, mar, merra]; others [mar,
merra]); report context floor 2 d / cap 366 d, detail-figure budget ~6,
overview layout 8×1. Merge check: mechanical diff against `setopts`.
The dedicated native flat-run screen owns its thresholds: minimum 3 d,
daily tair range ≤0.5 K, swd ≤5 W m⁻² while TOA ≥100 W m⁻², and lwd
within 15 W m⁻² of the temperature blackbody.

---

## D. Decision log (authoritative rulings; history detail in the ExecPlan)

- D-G (2026-07-24, Gate A): global rules, variable table, validation
  design and thresholds, triage, deviations ledger, alias crosswalk;
  MODIS then-dormant; standing autonomy granted (REVOKED 2026-07-26, see
  approval record).
- D-0a (2026-07-26): boom-height fail-closed REVERSED → A3 chain.
- D-0b (2026-07-26): solver keeps zero-rain until advective-rain physics
  is finished; product keeps real rainf.
- D-0c (2026-07-26): strictness dispositions — wind-beats-MAR and window
  honesty kept; zero-overlap last resort softened (A11); MAR/MERRA-only
  accepted (A12).
- D-1…D-7 (2026-07-25/26): runnability principle; window gating;
  single-year proxy guarantee; measured-ratio donor eligibility; native
  rainf preservation; figure pooling (superseded by D-19); producer
  pinning.
- D-8/D-9 (2026-07-26): two verdicts; product-level RCM fallback for
  partial rain, solver zero-rain only with no proxy anywhere.
- D-10 (2026-07-26): builder darkness zeros = known values with darkness
  provenance.
- D-11 (2026-07-26): agent-proposed mechanisms ratified as Method
  Contract; all numerics demoted to Section C.
- D-12 (2026-07-26): Feb-29 pooling ratified (parameter note).
- D-13 (2026-07-26): fail-closed albedo provenance ratified with the
  A3/A7 runnability guarantee.
- D-14 (2026-07-26): competition vs last-resort fallback formally
  distinguished; per-channel last-resort order is a parameter.
- D-15 (2026-07-26): MODIS albedo tier ACTIVATED (full 2000–2019 source
  on S03); staging-time attachment via one SSOT helper. D-15a: winter
  constant replaced by the floored seasonal bridge.
- D-16 (2026-07-26): verdict channel sets confirmed; swu never required.
- D-17 (2026-07-26): the proxy window BOUNDS the filled product; donors
  never extend the span; automatic pickup on future staging.
- D-18 (2026-07-26): no reconstruction-time partitioning; runtime option
  `source` (default) | `threshold`.
- D-19 (2026-07-26): appendix = per-station 8-tile overview + ~6 detail
  figures.
- D-20 (2026-07-26): canonical upward-shortwave name `swu`; `usr` read
  alias.
- D-21 (2026-07-26): per-channel tier-1 caps; loosening requires held-out
  evidence; psfc flagged.
- D-22 (2026-07-27): byte-pinning clarified — manifest refresh (re-hash
  on-disk files) is the approved fast path for WIP regeneration; bypass
  never; re-import never required (A1).
- D-23 (2026-07-27): appendix detail figures are curated for method
  diversity — the layout is the 8-tile overview PLUS the per-method
  comparison role folded into the ~6 detail figures (A14).
- D-24 (2026-07-27): naming rule generalized — source-schema names die at
  the ingest boundary; icemodel names throughout icemodel-owned
  artifacts; no product-specific exceptions (A16).
- D-25 (2026-07-27): A11 one-source-per-outage KEPT STRICT on evidence
  (8-site probe): MAR is never missing in-window; residuals were caused
  by the lwd floor rejecting genuine cold-sky values (fixed — A15 floor
  100→40, values confirmed as Section-C parameters) and by the
  zero-overlap last-resort refusal already reversed by D-0c (bead .28).
  The mixed-source sliver case is negligible once those land.
- Net accumulation (2026-07-26): snowf is never derived from
  surface-height/accumulation records (noise coupling, consistent with
  the deviations ledger); an accumulation-driven model mode is separate
  future work.
- D-26 (2026-07-27, user ruling): lwd ceiling 400→470 W m⁻². Native
  warm-fjord observations of 406–451 W m⁻² are valid physics; observations
  are never clamped or censored — the ceiling exists to reject garbage,
  so it moves above real physics.
- D-27 (2026-07-27, user ruling): calibrated rh candidates clamp into
  [5, 100]. Values above 100 after an overlap correction are calibration
  arithmetic, not physics; refusing them left megasample gaps (SWC 4.5%).
  Clamps are audit-noted.
- D-28 (2026-07-27, user ruling): swd relational validity. The candidate
  ceiling is max(1.05×TOA, 5) W m⁻² — the absolute term is a minimum
  ceiling during darkness (thermal-offset noise up to ~1.4 W m⁻² and
  twilight diffuse light of 8–28 W m⁻² are real incident energy; the
  ceiling additionally carries a solar-elevation twilight allowance where
  implemented). The twilight band is classified from the maximum solar
  elevation over the complete posting interval, matching PROMICE staging
  and reconstruction darkness; a posting that begins in deep night but
  enters twilight is not judged from its start instant. Relational
  exceedances (swd vs TOA, swu vs swd) become
  ledger DIAGNOSTIC columns and never flip verdicts: verdicts grade
  completeness + scalar bounds so real usable years read ready.
- D-29 (2026-07-27, user ruling): for swd, calibrated RCM proxies rank
  ahead of climatology in the candidate preference (climatology's
  day-of-year median structurally inserts clearer-than-context days into
  cloudy weeks: cen fills +0.09 median CSI vs observed neighbors).
  Climatology remains the B8 admission baseline and the fallback of last
  resort in RCM-void eras.
- D-30 (2026-07-27, user ruling): the product's hourly→15-minute
  expansion uses mean-preserving disaggregation for FILLED samples
  (smooth curve whose hourly integral equals the source posting);
  OBSERVED samples remain exact zero-order-hold copies of their native
  postings (honest copies, A7). Staged artifacts are unchanged — no
  restaging; the change is surgical at the product expansion.
- D-31 (2026-07-27, user rulings batch): MODIS is preferred wherever it
  has retrievals for albedo (already first in the last-resort order);
  MERRA-2 precip adoptions derive rainf as the exact complement
  ppt − snowf from the staged MERRA snowf; IMAU s21/s22/s23 join the
  gap-filled cohort; per-method report panels must render only their own
  method in accent color; overview and detail figures live in separate
  folders and report sections.
- D-32 (2026-07-28, user ruling): forcing readiness, not residual
  characterization, is the completion objective. Run bounded interpolation
  again after last-resort adoption so newly bounded short seam slivers use
  the final finite neighbors as anchors. Preserve the B3 caps and
  SWD-aware geometry, validate every candidate, and audit its provenance.
  Where CSI remains undefined or blocked, use the heldout-supported
  flux-linear fallback for residual SWD rather than clipping to the
  instantaneous horizontal-TOA ceiling; keep that relation diagnostic and
  the scalar bound hard. Zero missingness alone is insufficient:
  reconstruction QA must also examine continuity, seam jitter, flat runs,
  physical validity, and cadence artifacts before a product is accepted.
- D-33 (2026-07-28, user-approved evidence ruling): use the narrow
  post-final SWD seam repair described by B6. On KANL, the observed
  99.9th-percentile screen found 75 of 5,690 method boundaries above the
  season/elevation-conditioned native-step envelope. One reconstructed
  posting per boundary and at most two iterative passes changed 106
  hourly postings and reduced the count to 1, below the expected
  six-boundary statistical tail; widening to two postings changed more
  data without improving the result. The shipped 15-minute product has
  zero missing/scalar-invalid required values and zero native mutations.
  All other required KANL channels had zero empirical boundary outliers,
  so no general smoother is justified by this evidence. KANM subsequently
  showed that the empirical pool must remain inside one season and fine
  solar-elevation bin while mathematical local continuity may cross only a
  fine-bin edge, never a coarse dark/twilight/daylight regime or season.
  The one-posting neighborhood of a coarse solar-regime transition is a
  physical dawn/dusk ramp, not a method seam; exact-zero provenance
  transitions require no empirical magnitude estimate. Native, raw,
  clamped, policy-known darkness, and D-44 twilight-climatology values are
  immutable in this pass. `[codex-amended]`
- D-34 (2026-07-28, user ruling): native-core coverage is a confidence
  advisory, not a readiness state. A low-coverage station-year must still
  be reconstructed throughout the A6 fillable window; completeness and
  scalar validity alone determine forcing readiness.
- D-35 (2026-07-28, user ruling): render the main and companion
  method-detail reports to both PDF and HTML, with PDF listed first. The
  main appendix embeds station overviews only; method/gap figures and
  their per-station candidate tables live in the companion so the main
  HTML remains responsive without discarding evidence.
- D-36 (2026-07-28, user-authorized evidence ruling): calibrate wind
  speed multiplicatively so a correction cannot create negative wind.
  On FRE's held-out 2018 observations, additive MAR calibration produced
  9.94% negative estimates and RMSE 3.50 m s^-1; the multiplicative
  correction produced zero invalid estimates and RMSE 1.50 m s^-1.
  This closes finite-proxy wind holes without a permissive bounds clamp.
- D-37 (2026-07-28, user-authorized evidence ruling): activate the B14
  native flat-run screen in the production driver for both targets and
  donors. FRE has one 59-day corroborated run, 2018-02-16 through
  2018-04-15, implicating tair, swd, and lwd; all implicated samples are
  reconstructed with zero observed provenance, the staged file remains
  unchanged, and the wind-calibration correction in D-36 makes all four
  FRE product years forcing-ready. The same screen rejects the suspected
  NUKU late-2011 case because its daily temperature range is not flat.
- D-38 (2026-07-28, user-authorized evidence ruling): the final report
  must answer "why is anything still unfilled?" exhaustively, not through
  aggregate counts or a representative example. Publish every reconciled
  unfilled run (site, year, channel, span, duration, final-tier reason)
  and every non-ready IceModel or snow-model station-year, and fail the
  report if those run rows do not exactly match the shipped missing
  samples. For a staged native station with no product, validate and show
  its native/proxy spans; overlapping spans with no product are an
  unexplained production failure, never an inferred
  `windowRecordDisjoint`. When no in-window gap remains, say so directly.
- D-39 (2026-07-28, extended 2026-07-29 from production evidence): extend
  only the SWD interpolation cap from 6 h to 9 h. The final DY2 pilot exposed six
  otherwise unfillable seven-hour daylight seam slivers after the strict
  whole-outage proxy pass. Observed-only KANL holdouts at the target
  durations show flux-linear RMSE 76.6 versus persistence 193.5 W m-2 at
  7 h (87,570 samples) and 91.0 versus 206.4 W m-2 at 8 h (96,720
  samples); shoulder-hour RMSE is 53.7 versus 131.6 and 68.2 versus
  140.9 W m-2, respectively. Final NAU production then exposed three
  otherwise-unfillable nine-hour SWD slivers (108 delivered samples).
  Across 15,329 disjoint nine-hour NAU honest-native holdouts, the exact
  production interpolation admitted 69,246 samples with RMSE 66.8 W m-2
  and MAE 50.4 W m-2 versus persistence RMSE 260.4 W m-2; a read-only
  nine-hour NAU replay left zero missing SWD, zero blocked years, and
  passed seam QA. The extension uses the existing post-final
  CSI-first/flux-linear path, season and darkness boundaries, scalar
  bounds, provenance, and audit. Every other channel remains capped at
  its separately evidenced limit.
- D-40 (2026-07-28, user-authorized evidence ruling): admit linear
  interpolation for albedo gaps no longer than 2 h. The final SWC pilot
  exposed four 1--2 h seam remnants rejected after an otherwise complete
  MAR adoption. Observed-only SWC holdouts show linear RMSE 0.0170 versus
  persistence 0.0237 at 1 h (6,289 samples) and 0.0167 versus 0.0273 at
  2 h (8,734 samples). The same scalar bounds, season boundary, native
  immutability, provenance, and audit apply; `swu` remains derived from
  final albedo and SWD rather than interpolated independently.
- D-41 (2026-07-28, user-authorized evidence ruling): extend the albedo
  interpolation cap from 2 h to 30 h. The first full-cohort run exposed
  eleven 3--26 h CP1 albedo remnants where the selected whole-outage MAR
  source itself lacked valid albedo. CP1 observed-only holdouts scored
  only finite observations inside each synthetic gap: at 18 h, linear
  RMSE was 0.0463 versus persistence 0.0602 (4,759 samples); at 24 h,
  0.0490 versus 0.0531 (14,638); at 26 h, 0.0487 versus 0.0520
  (13,631); and at 30 h, 0.0464 versus 0.0574 (9,907). The bridge still
  requires finite bracketing anchors within one season, enforces scalar
  bounds, preserves native values, and records bounded-interpolation
  provenance. This closes source-internal day-scale holes without
  relaxing any other channel.
- D-42 (2026-07-29, user-authorized evidence ruling): extend only the RH
  interpolation cap from 6 h to 8 h. The final JAR production run exposed
  two eight-hour residual RH seam slivers after otherwise available proxy
  candidates failed post-blend scalar validity. Observed-only JAR holdouts
  at the target duration show linear RMSE 3.341 versus persistence 5.604
  percentage points over 155,432 samples, with 2.212-point MAE and
  0.011-point bias. The existing post-final bridge still requires finite
  bracketing anchors in one season, preserves native samples, enforces the
  RH scalar bounds, and records bounded-interpolation provenance. Every
  other non-SWD/non-albedo channel remains capped at 6 h.
- D-43 (2026-07-29, user-authorized evidence ruling; extended by Codex
  under the user's scientific-autonomy direction): permit an SWD gap
  no longer than 3 h to use one adjacent finite anchor across a
  calendar-season boundary or to straddle that boundary itself. The strict
  all-anchor season check exposed five otherwise-unfillable NAU postings
  on August 31 (1998, 2013, and 2014); their MAR/MERRA candidates are
  invalid or absent, so the September 1 finite native value is the only
  local physical anchor. Across 211 honest-native NAU samples withheld at
  season boundaries, the exact local bridge has RMSE 11.45 W m-2 at 1 h
  and 19.64 W m-2 at 2 h, versus persistence RMSE 52.80 and 90.09 W m-2.
  At September boundaries specifically, bridge RMSE is 22.80 versus 76.14
  W m-2 and MAE is 17.30 versus 64.94 W m-2. Two final NSE gaps then
  exposed the three-hour May 31 edge class. Across 45 honest-native NSE
  three-hour boundary cases (135 samples), bridge RMSE is 23.61 versus
  persistence 125.97 W m-2; at June boundaries it is 31.41 versus
  242.40 W m-2. NUKL then exposed one two-hour May/June gap whose missing
  samples occupy both labels even though the sunlight curve is continuous.
  Across 105, 225, and 341 honest-native NUKL boundary-straddling cases,
  the exact bridge RMSE is 6.31, 9.54, and 17.21 W m-2 at 1, 2, and 3 h
  versus persistence RMSE 12.22, 28.21, and 51.84 W m-2. The exception is
  SWD-only, uses the central 3 h ceiling, never crosses a darkness boundary,
  and retains scalar validity, provenance, and audit.
- D-44 (2026-07-29, user-authorized evidence ruling): after post-final
  interpolation, a single still-missing SWD posting inside civil twilight
  with exactly one adjacent all-interval darkness posting may use the
  existing station day-of-year/posting median climatology. Support comes
  only from untouched native SWD, must contain at least 30 samples, and
  the candidate must stay within scalar bounds and the 50 W m-2 twilight
  ceiling. Correcting reconstruction darkness from posting-start to
  full-interval geometry exposed 47 such NSE edge runs. Leave-one-year-out
  NSE validation over 9,516 honest-native targets gives climatology RMSE
  3.14 W m-2 and MAE 1.30 versus a zero estimate's RMSE 3.72; dawn is
  effectively zero while dusk improves from 5.31 to 4.48 W m-2 RMSE.
  The tier uses dedicated append-only twilight-climatology provenance,
  preserves the 50 W m-2 ceiling through mean-preserving 15-minute
  expansion, and emits exact audit. It never fills a deep-dark posting,
  never changes native values, and declines longer, non-edge, or
  under-supported gaps.
- D-45 (2026-07-29, Codex evidence-backed ruling under the user's
  authorized scientific-autonomy direction): after post-final
  interpolation, a 2--9-posting SWD residual whose missing postings all
  remain above civil twilight may use one finite daylight anchor's
  clear-sky index when the opposite adjacent posting is a finite,
  policy-known darkness zero. The candidate applies that CSI to each
  target posting's TOA irradiance, uses the configured CSI-darkness
  threshold, and must pass the shared scalar and TOA-relative SWD validity
  checks. It never fills a deep-dark posting; one-posting edges remain
  owned by D-44. On honest-native DY2 holdouts, this one-sided CSI method
  beats linear bridging and zero at every tested 2--9-posting duration:
  sunset RMSE is 10.78 versus 29.59 and 28.63 W m-2 at two postings and
  68.87 versus 116.72 and 292.47 at nine; dawn RMSE is 3.89 versus 12.29
  and 3.89 at two and remains lower through nine. The exact DY2 residual
  class closes all 32 delivered 15-minute holes without changing native
  observations.
- D-46 (2026-07-29, Codex evidence-backed QA ruling under the user's
  authorized scientific-autonomy direction): absence of the configured
  minimum native reference count does not itself fail seam QA. A record with
  no method boundary needs no empirical boundary reference, and a boundary
  proven continuous by the rules below remains supported even when the
  record-wide native-reference count is zero. A nonzero SWD method boundary
  is supported when its jump is no larger than
  either an adjacent same-direction native-to-native step or the
  unavoidable linear slope between the two finite native anchors enclosing
  one contiguous synthetic run. The same exemption applies to an exact
  bounded-interpolation run between any two finite anchors; every
  interpolated value must equal the line between those anchors. This proof
  is boundary-local: the enclosing line may begin in another coarse solar
  regime, but the diagnosed boundary itself must still be same-season,
  same-fine-bin, same-regime, and outside the one-step physical-regime-edge
  neighborhood. An adjacent synthetic step alone is never independent
  support, and a nonlinear run remains unsupported. Empirical limits use
  the fine season/elevation bin
  first; only when that bin misses the configured reference minimum may
  they fall back to native steps in the same season and the same coarse
  dark/twilight/daylight regime, still meeting the identical minimum. This
  resolves DY2's 2010-02-28 monotonic
  14.2 -> 9.1 -> 4.0 W m-2 sequence, whose two 5.1 W m-2 steps were
  physically continuous despite only 30 references in a stratum whose
  configured minimum is 100, plus 18 EGP MAR-to-bounded-interpolation
  label transitions whose flux values lie exactly on their enclosing MAR
  lines. At FRE, six exact mixed-anchor bridges and two direct MAR/native
  steps (7.86 and 12.31 W m-2) are supported by 102 same-winter/daylight
  native steps with a 99th-percentile magnitude of 65.26 W m-2 when their
  fine bins contain only 58 and 27 references. Season and coarse-regime
  guards remain mandatory. At KPCU, two zero-deviation observed-to-MAR
  bridges begin in twilight but end two postings into an eligible daylight
  boundary; their remaining jumps are only 0.35 and 0.72 W m-2.
- D-47 (2026-07-29, Codex evidence-backed completion ruling under the
  user's scientific-autonomy direction): guarded hourly or half-hourly
  reconstruction must rerun the existing D-32 residual closer after
  mean-preserving restoration to the delivered 15-minute axis. This is not
  a new method: it retains the B3 cap, bounded-interpolation provenance,
  finite-anchor and extraterrestrial-anchor guards, scalar bounds, and
  exact audit, while native held copies remain immutable. Re-derive missing
  SWU immediately afterward and reconcile provisional unfilled rows on the
  delivered axis. NUKL exposed eleven one-hour low-sun source postings that
  the coarse posting could not certify; their 44 delivered quarters form
  smooth local bridges. Together with D-43's two-hour May/June correction,
  the fresh NUKL replay moves from 52 missing SWD quarters to zero, makes
  all 13 IceModel and snow-model years ready, leaves zero required unfilled
  audit rows, and passes seam QA with zero unsupported boundaries and zero
  outliers.
- D-48 (2026-07-29, Codex evidence-backed completion ruling under the
  user's scientific-autonomy direction): RH interpolation may cross or
  straddle a calendar-season boundary within the existing D-42 eight-hour
  cap. The calendar label is not a physical humidity discontinuity; all
  other B3 bounds, finite anchors, native immutability, provenance, and
  audit remain unchanged, and no other state channel inherits this
  exception. SCOL exposed one otherwise-unfillable hour at the August/
  September boundary, bracketed by 99.71% and 89.51% RH. Across 135, 315,
  495, 675, 855, 1035, 1215, and 1395 honest-native SCOL boundary cases
  at 1 through 8 h, respectively, linear interpolation beats persistence
  at every duration: bridge RMSE rises from 1.37 to 4.97 percentage points
  versus persistence from 2.65 to 7.64.
- D-49 (2026-07-29, Codex evidence-backed completion ruling under the
  user's scientific-autonomy direction): apply SWD's existing D-39
  nine-hour interpolation cap at calendar-season boundaries instead of
  maintaining a separate three-hour ceiling. Missing samples must remain
  above civil twilight; finite anchors, the extraterrestrial-anchor guard,
  scalar validity, native immutability, bounded-interpolation provenance,
  exact audit, and final seam QA remain mandatory. SDM exposed one
  otherwise-unfillable five-hour May/June ramp, entirely above civil
  twilight and bracketed by 285.9 and 30.4 W m-2. Across 24 to 256
  honest-native SDM boundary-straddling cases per duration, linear bridging
  beats persistence at every 1--9 h duration: bridge RMSE rises from 23.10
  to 119.17 W m-2 versus persistence from 44.73 to 383.20; at the target
  five hours it is 69.21 versus 249.87. The calendar-edge cap is an alias
  of the central SWD cap so the two cannot drift.
- D-50 (2026-07-29, Codex evidence-backed completion ruling under the
  user's scientific-autonomy direction): extend RH's D-42 interpolation
  cap from eight to nine hours. UPE_L exposed one otherwise-unfillable
  nine-hour RH residual, bracketed by 94.79% and 96.98%; the only final-tier
  denial was the eight-hour ceiling. Across 87,884 honest-native UPE_L
  nine-hour windows (3,163,824 samples), the exact linear bridge has RMSE
  5.20 and MAE 3.56 percentage points versus persistence RMSE 8.28 and MAE
  5.58. Finite bracketing anchors, the calendar-label rule, RH scalar
  bounds, native immutability, bounded-interpolation provenance, exact
  audit, and final seam QA remain mandatory; ten-hour gaps remain outside
  tier 1.

---

## E. Reference material

### Evidence base (2026-07-23 snapshot)

- Gap census v2, 51 stations, in-record, daylight-only swd:
  `data/preview/qa/gapfill/promice_gap_census.csv`. Headline: six-hour
  interpolation fixes 2.3% (tair), 2.3% (rh), 4.9% (wspd), 3.1% (psfc),
  0.4% (lwd), 18.5% (swd) of missing samples; lwd is 51.7% missing with
  ~1,500 runs longer than a day; albedo is 3.6% missing in 38 long runs.
- Proximity survey: AWS5→KAN_L 6.4 km (−168 m); AWS6→KAN_L 23.8 km
  (+332 m); AWS6→KAN_M 24.4 km (−269 m); AWS9→KAN_M 25.4 km (+231 m);
  AWS10→KAN_U 0.2 km (+10 m); KAN_M elevation-bracketed by S6/S9; IMAU
  S21–S23 hundreds of km away.
- Compatibility audit: legacy consumed PROMICE v03 ASCII; current staging
  consumes pypromice L3 hourly NetCDF. Never port: −999 handling, column
  maps, `RHice2water`/NUK_K patch, KAN_U daily-v02 splice, hard-coded
  sensor-height offsets, in-code tilt/bias correction.
- Vandecrux triage (19 repos): no Charalampidis S10→KAN_U implementation
  anywhere; KAN_U only a DONOR in legacy code; fallbacks were
  RACMO/HIRHAM5; MODIS MOD10A1 C6 was the SW-up filler; AWS_processing is
  GPL-3.0; SEB_Firn_model has no license.
- MODIS: full GEUS C6 daily coverage 2000–2019 staged at
  `/Volumes/S03/DATA/greenland/geus/albedo/gris`; the repo-local cache
  holds only the 2012 v1.1 fixture. (Supersedes the earlier "2012-only"
  inventory note; tier activated per D-15/B12.)
- Post-cycle-62 state audit (Codex):
  `.agents/plans/drafts/2026-07-26-promice-gapfill-pre-policy-hardening-audit.md`.

### Variable-by-variable notes (supplementing A/B; canonical L3 mappings)

- **tair** — `t_u` → `tair` (K). Donors with lapse-rate elevation
  adjustment (overlap-fitted; recorded fallback). Proxy: MAR t2m, MERRA-2
  alternative.
- **rh** — `rh_u` (wrt water, upstream-corrected) → `rh` (%). No default
  elevation adjustment. Proxy: MAR rh (derive from q2m/t2m if not
  direct). Clamp [5, 100].
- **wspd** — `wspd_u` → `wspd`. Donors admitted only where held-out RMSE
  beats calibrated MAR wind in the same stratum `[codex-amended]`.
  `wdir`: filled only from the same source as wspd (vector components),
  never independently interpolated beyond the cap; low-speed masking for
  K-transect direction.
- **psfc** — `p_u` ×100 → `psfc` (Pa). Barometric elevation adjustment
  mandatory for donors.
- **swd** — `dsr_cor` (fallback `dsr`) → `swd`. Never raw-linear in
  daylight; CSI space throughout; darkness zero-fill (B2). Screened
  `dsr_cor` with raw `dsr` present is a distinct sub-case: prefer
  re-screened raw over model proxy when tilt flagging, not the sensor,
  caused the gap (per-gap QC signature). Raw fallback code 13; clamped
  negatives code 14.
- **swu / albedo** — `usr_cor`/`albedo` → `swu`, `albedo` (naming per
  A16). swu = albedo × swd after operand tiers (B10). Albedo masking and
  refill per A7; winter bridge per B13; MODIS per B12.
- **lwd** — `dlr` → `lwd`. The critical channel (51.7% missing). Donors
  where instruments exist (all four K-transect stations); proxies MAR and
  MERRA-2 in competition with the `fill_lwd`-style vapor-pressure
  estimator; `lwd_estimated=true` artifacts masked before planning.
- **ppt / rainf / snowf** — per A10. `CalculatePrecipitation`-style
  height-change scaling is NOT adopted (couples precip to surface-height
  noise; recorded deviation).
- **instrument heights** — `z_boom_cor_u` → `boom_height`, time-varying,
  with the A3 fallback chain (measured → interpolated → 2.6 m).
  K-transect donor heights from staged per-visit records (units to be
  confirmed against Smeets et al. 2022 — magnitudes are cm).

### Deviations-from-papers ledger (publication bookkeeping)

- MAR replaces HIRHAM5 as the complete-met fallback; CARRA noted future.
- Spline interval count is a validated hyperparameter, not fixed.
- Longwave bias correction estimated from our own overlaps.
- No use of `KAN_U_babis.csv` (derived, unreproducible).
- `ResampleTable`-style whole-channel extrapolation and `FillLastGaps`
  unrestricted interpolation explicitly rejected.
- Precipitation never from stations nor accumulation-to-height scaling;
  MAR/MERRA-2 totals with the runtime phase option (A10).
- Provenance registry distinguishes every method (no legacy
  conflations).

### GC-Net Vandecrux reference product (approved 2026-07-24)

The cached Vandecrux GC-Net product
(`data/verification/gcnet/<station>_surface.nc`; nine stations: CP1,
DYE_2, NASA_E, NASA_SE, NASA_U, Saddle, SouthDome, Summit, TUNU_N;
Vandecrux et al., J. Glaciol.,
https://www.cambridge.org/core/journals/journal-of-glaciology/article/7ECF88C9F3A8CD32C5CA4C333643079E)
is an already-gap-filled forcing series with per-variable `_origin`
codes. Treatment: (1) donor role only for origin-observed samples, per
channel (A8); (2) reference-reconstruction role — our held-out fills are
compared against Vandecrux's on common gaps and tabulated in the report
(divergences are diagnostics, not failures); (3) origin-observed samples
may extend donor records backward under the same gates; (4) the nine
stations cross-reference the worth-filling triage for accumulation-zone
targets.

### Approval record

- 2026-07-24 (Gate A): original approval of the draft policy (global
  rules, variable table, MODIS-dormant verdict, validation design,
  triage, deviations ledger, crosswalk), with standing autonomy.
- 2026-07-26 (hardening pass): the user reviewed the adversarial
  hardening draft (v1 comments + v2 sign-off), ratifying decisions
  D-0a…D-21 and this restructured contract. The 2026-07-24 grant of
  standing autonomy for policy changes is REVOKED: every policy change
  now requires a Section-D decision-log entry recording the user's
  ruling. Parameters (Section C) remain freely tunable via `setopts`.
- 2026-07-28 (forcing-readiness correction): the user explicitly
  authorized evidence-driven policy changes needed to make every
  in-window product forcing-ready and publication-grade. Each such
  change still requires its own Section-D evidence entry.
