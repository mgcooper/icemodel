# Legacy Test Helpers

This folder holds historical implementations. They are outside the core
production path, but they are still useful as reference code for one-off
studies and benchmark comparisons.

Current contents:

- `SPECTRALSOURCETERM_INLINE.m`
  - historical fully inlined spectral source-term implementation
  - kept so the spectral perf study can compare the organized-functions path
    against the inline organization
