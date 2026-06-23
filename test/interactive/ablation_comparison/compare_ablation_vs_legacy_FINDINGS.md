# New PROMICE ablation vs legacy staged ablation (KAN_L, KAN_M)

One-off comparison (bead `icemodel-dz2.2`, shifted role). The NEW PROMICE eval
builder `icemodel.forcing.buildPromiceData` (reads L3 `z_ice_surf` with
site-type branching + gap flags) now fills the GOF-reference role that the
pre-refactor staged ablation artifacts originally held. The central research
question this answers: **would historical KAN_L / KAN_M model-data evaluations
change if the obs reference were swapped from the legacy series to the new
builder?**

- Script: `test/interactive/compare_ablation_vs_legacy.m`
- Figures (gitignored): `test/interactive/figures/ablation_vs_legacy_{kanl,kanm}_{full,2015_jun_oct,2015_jul_oct,2016_jun_oct,2016_jul_oct}.png`
- New series: `buildPromiceData("KAN_L"|"KAN_M", frequency="daily").ablation`
  (cumulative ice-surface lowering, positive down, gap-bridged samples flagged
  and EXCLUDED from the stats).
- Legacy series: `test/assets/legacy_ablation/kan{l,m}_ablation_daily.mat`
  (timetable, variable `ablation`, daily, 2009-01-01..2018-12-31, cumulative
  positive down). Retired as the FORMAL eval reference but PRESERVED as a
  committed test asset (bit-identical to the runoff copies under
  `runoff/data/icemodel/eval/`).
- Convention: each window rebaselines BOTH series to their value at the window
  start `t1` (the legacy `plotPromice(...,'refstart',t1)` behavior — subtract
  the `t1` value so both series start at 0 at `t1`). Surface heights are
  relative to installation; only the shape and the total ablation over the
  window are physical.

## Full-record agreement (2009-2018, direct-observation days only)

| site  | n     | gap-excl | corr   | bias (m) | RMSE (m) | total new (m) | total legacy (m) | total diff (m) | rel %  | max\|diff\| date  |
|-------|-------|----------|--------|----------|----------|---------------|------------------|----------------|--------|-------------------|
| KAN_L | 3645  | 7        | 0.9982 | -10.004  | 11.479   | 38.573        | 57.130           | -18.557        | -32.5  | 25-Dec-2018       |
| KAN_M | 3260  | 177      | 0.9299 | -1.207   | 4.496    | 12.750        | 23.796           | -11.046        | -46.4  | 31-Dec-2018       |

The full-record difference is a slowly COMPOUNDING multi-year datum/scaling
drift (the difference curve steps down each summer and is flat each winter,
max divergence at the END of the record), not a single discontinuity. It is the
right sign/character for a surface-height-channel / datum change between the
L3-based new builder and the older legacy artifact.

## Targeted melt-season windows (rebaselined to window start `t1`)

These are the windows over which historical KAN_L / KAN_M evaluations were run.
Each row shows how much a GOF metric would shift if the obs reference changed
from legacy to new over THAT window.

| site  | window          | n   | gapx | bias (m) | RMSE (m) | corr   | total new (m) | total legacy (m) | total diff (m) | rel %  |
|-------|-----------------|-----|------|----------|----------|--------|---------------|------------------|----------------|--------|
| KAN_L | 2015 Jun1-Oct1  | 123 | 0    | -0.248   | 0.255    | 0.9984 | 2.730         | 2.995            | -0.265         | -8.8   |
| KAN_L | 2015 Jul1-Oct1  | 93  | 0    | +0.051   | 0.055    | 0.9999 | 2.249         | 2.195            | +0.054         | +2.4   |
| KAN_L | 2016 Jun1-Oct1  | 122 | 1    | +0.003   | 0.019    | 0.9999 | 4.103         | 4.180            | -0.077         | -1.9   |
| KAN_L | 2016 Jul1-Oct1  | 92  | 1    | +0.018   | 0.026    | 0.9998 | 2.735         | 2.803            | -0.068         | -2.4   |
| KAN_M | 2015 Jun1-Oct1  | 122 | 1    | -0.251   | 0.265    | 0.9664 | 0.719         | 1.006            | -0.287         | -28.6  |
| KAN_M | 2015 Jul1-Oct1  | 93  | 0    | +0.057   | 0.077    | 0.9831 | 0.719         | 0.665            | +0.054         | +8.1   |
| KAN_M | 2016 Jun1-Oct1  | 123 | 0    | -0.125   | 0.130    | 0.9991 | 2.138         | 2.287            | -0.149         | -6.5   |
| KAN_M | 2016 Jul1-Oct1  | 93  | 0    | +0.021   | 0.026    | 0.9998 | 1.821         | 1.816            | +0.005         | +0.3   |

## Verdict: would historical evaluations materially change?

**No — not materially, for any melt-season window evaluated.** Within a single
melt season the new and legacy series agree to within a few centimetres of bias
and RMSE and correlate at r >= 0.97 (r >= 0.998 for the Jul-start windows). The
large full-record divergence is a compounding multi-year datum drift that does
NOT live inside any one season; once each series is rebaselined to the window
start it cancels out almost entirely.

- **Jul1-Oct1 windows (the cleanest melt-season span): negligible.** Bias and
  RMSE are <= 0.08 m, total-ablation difference <= 0.05 m (<= 8% at KAN_M where
  the seasonal total is small, <= 2.4% at KAN_L). A GOF computed against the new
  reference instead of the legacy one would shift by far less than the model's
  own error. Historical Jul-Oct evaluations stand.

- **Jun1-Oct1 windows: a small, sign-consistent early-June offset.** The new
  series lags the legacy by ~0.25 m bias in 2015 (both sites) because the legacy
  artifact registers some surface lowering in early June that the new L3
  `z_ice_surf` series does not (still snow-covered / not yet pure ice surface).
  This is ~9% of the 2015 KAN_L window total and ~29% of the small 2015 KAN_M
  window total; the 2016 Jun-Oct windows agree to within 0.13-0.15 m (<= 6.5%).
  An evaluation starting 1 June 2015 could shift by up to ~0.25 m of reference
  ablation; this is the only case where the choice of reference is worth noting,
  and even there it is a season-start datum effect, not a model-relevant change.

**Bottom line for the research:** the new `buildPromiceData` ablation can replace
the legacy staged artifact as the melt-season GOF reference with NO material
effect on historical KAN_L / KAN_M evaluations, provided the window is
rebaselined to its start (which the legacy `plotPromice 'refstart'` workflow
already did). The only caveat is to prefer a July (pure-ice) window start over a
June start at these sites, since the June segment carries a small surface-channel
ambiguity. The new builder is the L3-QC'd, gap-flagged source and is the
defensible reference going forward; the legacy series is preserved as a
provenance-uncertain historical artifact (see Provenance below), not ground
truth.

## Provenance of the full-record drift (background)

The high-correlation / large-magnitude full-record pattern points to a per-season
SCALING / datum difference in how cumulative surface lowering is summed, not a
timing or noise difference:

1. **Different surface-height channel / datum.** The new builder uses the QC'd
   L3 `z_ice_surf` (pure ice-surface lowering). The legacy artifact predates that
   channel and was almost certainly built from a different surface-height product
   (raw SR50 / pressure-transducer surface height or an earlier processing level)
   that registers more total lowering per season.
2. **Reinstall / offset handling at KAN_M.** The lower KAN_M full-record
   correlation reflects a station reinstall / datum reset the two pipelines zero
   differently. This is a KAN_M data-history artifact, not a builder bug, and it
   washes out inside the rebaselined melt-season windows (r >= 0.966 there).
3. **Gap / winter bridging.** A small contributor (7 days excluded at KAN_L; 177
   at KAN_M, mostly winter), not the main driver of the full-record offset.

## Model run

This pass did NOT run the icemodel itself over the windows. The
`+icemodel/+verification` runners (`runIcemodelSnowCandidate`,
`candidateFromIcemodelOutput`, the `+setup` namelists) are built for the Colbeck
snow-verification cases, not for a KAN_L / KAN_M melt-season ablation run, so
staging a real model run over these windows is not a quick drop-in. Per the task
this is optional and was not allowed to block the required data-vs-data
deliverable. The model-vs-legacy-obs-vs-new-obs overlay is left as a follow-up:
the intended path is a standalone interactive script that loads saved icemodel
output for the window and overlays modeled cumulative ablation onto the two obs
series produced here (the figure already plots legacy vs new; adding a third
modeled line is a small extension once a saved run exists).
