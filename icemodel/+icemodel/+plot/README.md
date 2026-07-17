# icemodel.plot

Shared plotting helpers for model forcing, verification time series, profiles,
and comparisons.

`sourceColor` keys colors to canonical source labels, so plot order and active
subsets cannot change source identity. MAR, RACMO, and MERRA-2 use the exact
`runoff/functions/RunoffPlot.m` palette; PROMICE model met, PROMICE native
userdata, and observations use stable distinct role colors.

`timeseries` preserves explicit NaNs and inserts a NaN midpoint when repeated
cadence makes an omitted-time gap unambiguous. It does not infer gaps for
two-point interval observations or irregular sparse series.
