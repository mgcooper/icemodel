# New PROMICE ablation vs legacy staged ablation (KAN_L, KAN_M)

One-off comparison (bead `icemodel-dz2.2`, shifted role). The NEW PROMICE eval
builder `icemodel.forcing.buildPromiceData` (reads L3 `z_ice_surf` with
site-type branching + gap flags) now fills the GOF-reference role that the
pre-refactor staged ablation artifacts originally held.

- Script: `test/interactive/compare_ablation_vs_legacy.m`
- Figures (gitignored): `test/interactive/figures/ablation_vs_legacy_{kanl,kanm}.png`
- New series: `buildPromiceData("KAN_L"|"KAN_M", frequency="daily").ablation`
  (cumulative ice-surface lowering, positive down, gap-bridged samples flagged
  and EXCLUDED from the stats).
- Legacy series: `demo/data/eval/kan{l,m}_ablation_daily.mat` (timetable,
  variable `ablation`, daily, 2009-01-01..2018-12-31, cumulative positive down).
- Both rebaselined to their value at the common-window start (surface heights
  are relative to installation; only shape + total are physical).

## Agreement stats (direct-observation days only)

| site  | n     | gap-excl | corr   | bias (m) | RMSE (m) | total new (m) | total legacy (m) | total diff (m) | max\|diff\| (date)   |
|-------|-------|----------|--------|----------|----------|---------------|------------------|----------------|----------------------|
| KAN_L | 3650  | 2        | 0.9982 | -10.00   | 11.48    | 38.57         | 57.13            | -18.56         | 18.56 (25-Dec-2018)  |
| KAN_M | 3422  | 15       | 0.9261 | -1.01    | 4.43     | 12.75         | 23.80            | -11.05         | 11.05 (31-Dec-2018)  |

Common window: 2009-01-01 .. 2018-12-31 (the legacy record span).

## Findings

**KAN_L — same shape, smaller magnitude per melt season (high correlation,
large compounding offset).** The two series track the SAME seasonal staircase
(summer melt steps, winter plateaus), so the correlation is 0.998. But the new
series accumulates noticeably LESS ablation each melt season; the difference
curve steps down every summer and is flat every winter, growing monotonically
to -18.6 m by 2018. This is a per-melt-season systematic magnitude offset
(new is roughly two-thirds of legacy per season), NOT noise and NOT a single
step. The maximum divergence is at the END of the record, confirming
progressive (compounding) divergence.

**KAN_M — datum / reinstall discontinuity plus later divergence (lower
correlation).** The new series rises from 2009 while the legacy stays near zero
until mid-2012, so early on the NEW series is AHEAD of legacy (diff up to
+3.7 m). After ~2015 the legacy overtakes and ends ~11 m higher (diff -11 m).
The crossover and the early-period legacy flatness indicate the two series are
not measuring the same quantity over the full record at KAN_M (different start
datum / surface channel, and almost certainly a station reinstall the two
pipelines handle differently). Correlation is only 0.926 because of the
crossover.

## Likely cause

The high-correlation/large-magnitude pattern at KAN_L points to a per-season
SCALING difference in how surface lowering is summed, not a timing or noise
difference. The most probable causes, in order:

1. **Different surface-height channel / datum.** The new builder uses the QC'd
   L3 `z_ice_surf` (pure ice-surface lowering). The legacy artifact predates
   that channel and was almost certainly built from a different surface-height
   product (raw pressure-transducer / SR50 surface height, or an earlier
   GEUS/PROMICE processing level) that registers more total lowering per season
   (e.g. it folds in snow-surface drop or uses an un-QC'd height that drifts).
2. **Reinstall / offset handling at KAN_M.** The KAN_M crossover and the
   legacy's flat 2009-2012 segment look like a station reinstall or datum reset
   that the new L3-based builder and the legacy pipeline zero differently. This
   is a KAN_M-specific data-history artifact, not a builder bug.
3. **Gap / winter bridging.** The legacy series may continue accumulating
   (slope-bridged) across winter gaps where the new builder holds flat and
   flags the gap; only 2-15 days are flagged here, so this is a SMALL
   contributor, not the main driver.

## Bottom line

The new `buildPromiceData` ablation reproduces the SHAPE and seasonal timing of
the legacy reference very well (especially KAN_L, r=0.998), but the cumulative
MAGNITUDE differs substantially and the difference is the right sign and
character for a surface-height-channel / datum change between the L3-based new
builder and the older legacy artifact. Before adopting either as the GOF
reference, the legacy artifact's provenance (which surface-height channel and
processing level it was built from, and the KAN_M reinstall handling) should be
pinned down. The new builder is the L3-QC'd, gap-flagged source and is the
defensible reference going forward; the legacy series should be treated as a
provenance-uncertain historical artifact, not ground truth.
