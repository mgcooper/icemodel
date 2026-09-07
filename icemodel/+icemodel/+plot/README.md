# icemodel.plot

Purpose: Shared plotting helpers.

Contents:

- `sourceColor` keys colors to MAR, RACMO, and MERRA-2 (based on the
  `runoff/functions/RunoffPlot.m` palette). PROMICE met, userdata, and
  observations have their own palette.
- `timeseries` preserves NaNs and inserts a NaN midpoint for an unambiguous gap
  in a regularly sampled series.
- `newFigure` creates the shared hidden, white, export-sized frame.
- `markTimeSpan` adds interval boundaries that stay out of the legend.
- `formatDuration` formats hour, day, and year labels for report figures
  and tables.
