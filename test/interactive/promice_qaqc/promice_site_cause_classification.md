# PROMICE evaluation-data per-site cause classification

Companion to `promice_eval_summary.md` (auto-generated metrics),
`promice_step_screening.md` (all-site step screening), and the per-site
`promice_eval_<site>.png` figures. This file is the manual interpretation layer
explaining each per-site behaviour and stating what the builders do about it.

## Flagging philosophy (what is GEUS vs ours)

We preserve the authoritative GEUS series and ATTACH per-sample flags; we never
silently edit data. Concretely:

- **GEUS provides** all QC, manual flagging/fixing, gap-filling (slope-bridging
  the surface height across sensor outages), the multi-sensor `z_surf_combined`,
  the `z_ice_surf` re-derivation (one-year lagging minimum of `z_surf_combined`,
  driven by the pressure transducer / stake sonic ranger), `snow_height`, the
  depth-tagged thermistor string, and the standardized 10 m `t_i_10m`.
- **We provide** only: read the L3 channels; derive per `site_type`
  (`ablation`/`surface_height`); clamp `snow_depth >= 0`; unit conversions; and
  per-sample FLAGS (gap-bridged, station-transition window, step
  detected/correctable + magnitude). **No GEUS data value is modified.**

The de-stepping CORRECTION is an opt-in evaluation-time transform
(`icemodel.forcing.destepSurface`); the staged `.mat` is faithful, and only
UNAMBIGUOUS steps are corrected by default at consumption.

## The ablation discriminator (two authorities that agree)

A site is treated as an ablation site when its L3 file ships `z_ice_surf` (the
OPERATIONAL authority used by `buildPromiceData`). This agrees with the readme
Table 1 "Site type" surfaced by `promicesiteinfo.site_type`. Both are
documented and consistent; the operational `z_ice_surf`-presence test is the one
the builder branches on. (Marginal local glaciers whose L3 file ships no
`z_ice_surf` — e.g. ZAC_A, WEG_B — screen operationally as accumulation even
though Table 1 calls them Ablation; that is the honest operational state of
their staged data.) The accumulation/ablation + season relationship is USED as a
QA/QC evidence line in step detection (a melt-signed step in winter is suspect),
not as a relabeling.

## Guidance for model comparison (gap-bridged data)

Per the GEUS readme: *"surface height change over data gaps should not be
regarded as direct observation, still the surface height trend over the entire
period should be unaffected by the gaps."* Therefore:

1. **Subtract the window-start value** of the surface-height channel before
   comparing (all L3 surface heights are relative to installation; the builder
   already zeroes the channel at the first finite window sample).
2. **CUMULATIVE and visual comparison use the FULL series** — the trend is
   preserved through gaps. Do NOT drop gap-bridged samples from a cumulative or
   visual comparison.
3. **Only RATE-based diagnostics** (per-timestep ablation/accumulation rate)
   should flag/exclude `surface_height_flag == 1` segments, where the
   per-timestep rate is slope-interpolated rather than observed.
4. Use `tice10m` (heavy black line, panel c) as the primary subsurface
   comparison; the grey thermistor string is the depth-tagged diagnostic.
5. On accumulation sites there is no ablation/snow_depth channel by design.

## Per-site explanations

| Site | type (op.) | Behaviour | Explanation |
| --- | --- | --- | --- |
| CEN  | accumulation | drop ~end 2021 | Station merge (CEN2/CEN1/GITS) over a heavily gap-bridged historical GC-Net record. Emitted as `surface_height`; gap segments flagged. The change at the merge is expected. |
| CP1  | accumulation | oscillations 2010/11 + gaps | A real GEUS seasonal accumulation signal through a gap-heavy record: the historical `z_surf_combined` is slope-bridged across long gaps, so the apparent oscillations are mostly interpolation between sparse observed segments. The cumulative trend is still authoritative; RATE diagnostics should exclude the bridged stretches. |
| DY2  | accumulation | surface-height, not ablation | DY2 ships no `z_ice_surf`, so it carries only the net `surface_height` (rises ~20 m over 30 yr). No ablation channel applies. |
| EGP  | accumulation | surface step mid-2023 | EGP yields only `surface_height` (no `z_ice_surf`). The residual surface step is a real event in `z_surf_combined`; gaps are flagged. |
| KAN_B | accumulation (bedrock/off-ice) | snow < 0 late | KAN_B sits on tundra/bedrock OFF the ice: the ice-surface channel barely applies. It ships no `z_ice_surf`, so NO `snow_depth` is exported and the raw `snow_height` noise (775 sub-cm negatives) never reaches the eval axis. The surface series is withheld/flagged as off-ice. This remains authoritative GEUS data. |
| KAN_U | accumulation (percolation) | reference | Net `surface_height` -0.5..+7.6 m. `site_type`=Accumulation, `surface_zone`=percolation (firn-core truth). A few >1 m single-hour glitches in `z_surf_combined` are flagged unambiguous (transient spikes). |
| LYN_L | ablation | rapid rise then flat | A short local-glacier record: an early observed ablation segment (~0.1 m, real), then the record reaches a gap / record-end and flatlines. The flat tail is gap-bridged and flagged; the ablation magnitude is real. Remains authoritative GEUS data. |
| LYN_T | ablation | rapid rise then flat | Same pattern as LYN_L; record ends 2024-04, much of it gap-bridged and flagged. SHORT/SPARSE flags fire. |
| MIT  | ablation | ~12 m step early | Two ~5.9 m single-hour jumps at 2009-08-11/12 (~11.9 m total): grossly implausible, single-station (no merge), so an installation/sensor jump. These are the UNAMBIGUOUS correctable steps: the default `destepSurface` transform levels them at analysis time; the staged series stays faithful. |
| NUK_K | ablation | early bad thermistor; late-2024 step at a gap | Surfaced/early-install thermistor readings are removed by the depth<=0 discard + the [-80,1] clamp; `tice10m` is primary. The late-2024 surface step sits in a gap-flagged segment. Residual warm-share is early-record near-melt, flagged for the eye. |
| NUK_L | ablation | snow spike late 2011 | Ablation to ~98 m (low-elevation margin, plausible). The 2011 snow spike (max ~1.9 m) is within range and is real remaining GEUS noise; `snow_depth` is clamped >= 0, the spike is flagged, not deleted. |
| QAS_L | ablation | 2009 spike | `snow_max` 6.87 m fires SNOW_HIGH; ablation to ~98 m is real (very low elevation). Spike flagged, not deleted. |
| QAS_U | ablation | apparent step late 2023 | The apparent ~2 m offset spans a long gap (~0.03 m/day), so it is NOT a discrete step — the GAP flag covers it, and step detection correctly does NOT flag it. QAS_U is a known merge (QAS_Uv3/QAS_U), but with no staged per-station handover date the transition evidence cannot confirm a step here. Reported honestly, not auto-corrected. |
| SCO_U | ablation | snow step late 2025 | A near-real-time data artefact in the still-updating record; the readme notes recent values are recomputed after the next station visit. `snow_depth` clamped >= 0; flagged for the eye. |

## How to regenerate

```matlab
% from the repo root, on the staged L3 product:
T = plot_promice_eval(sites="all", frequency="daily", save_figs=true);
```

Writes `promice_eval_<site>.png` (3 panels: site-type surface channel; the
surface-height series with gap-bridged samples in red, station-transition
windows in green, detected steps marked unambiguous=magenta-square /
ambiguous=grey-x; `tice10m` heavy-black primary over the depth-tagged thermistor
string) plus `promice_eval_summary.md`. The all-site step screening is in
`promice_step_screening.md`. This file is the manual interpretation layer.
