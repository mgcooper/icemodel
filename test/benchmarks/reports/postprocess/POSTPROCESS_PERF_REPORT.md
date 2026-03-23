# Postprocess Performance Study

## Scope

This study compares the legacy timetable `retime(..., 'hourly', 'mean')`
branch against the fixed-step hourly aggregation helper used when the input
series is aligned quarter-hour output.

The study is separate from the spectral workflow.

## Production Code Shape

1. `/Users/mattcooper/MATLAB/projects/icemodel/icemodel/+icemodel/postprocess.m`
   - chooses the fixed-step fast path when the timetable spacing is aligned
   - falls back to the legacy timetable `retime` branch otherwise
2. `/Users/mattcooper/MATLAB/projects/icemodel/icemodel/+icemodel/+internal/retimeHourlyFixedStep.m`
   - exact fixed-step hourly aggregation for aligned quarter-hour data

## How To Reproduce

```matlab
clear functions
report = summarize_postprocess_perf(simyear=2016);
```

## Interpretation Notes

- The benchmark times only the hourly aggregation branch.
- Accuracy is reported as the maximum absolute difference against the legacy
  timetable `retime` output.
- This tool is intentionally separate from `run_spectral_study_bootstrap`.

## Resolution

The fixed-step hourly aggregation is accepted as the production default. The
legacy `retime` branch remains as a fallback for non-aligned timetable spacing.
