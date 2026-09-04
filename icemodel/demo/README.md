# icemodel demos

Live scripts in the plain-text live code format (`.m` with `%[text]`
markup), so they diff and version like ordinary code.

Contents:

- `demo_vaporDensity.m`: saturation vapor pressure and density
  formulations (Buck, Jordan/SNTHERM, Ambaum/Romps) and their
  derivatives
- `demo_vaporDiffusion.m`: the vapor thermal conductivity `k_vap` and
  the Liston/Jordan formulation difference
- `demo_vaporModel.m`: the vapor model pieces end to end
- `demo_vaporSurfaceBoundary.m`: why the interior vapor transport
  operator closes its top face, what face 1 carries instead, and the
  Patankar boundary treatment

Format requirement: the Live Editor recognizes a plain-text `.m` file as
live code only when the file ends with the appendix trailer

```
%[appendix]{"version":"1.0"}
%---
```

(`matlab.desktop.editor.EditorUtils.isLiveCodeFile` returns false
without it, and the file opens as an ordinary script with the `%[text]`
lines printed literally). The Live Editor appends `%[metadata:view]` and
`%[output:...]` blocks after that trailer when it saves outputs. A new
demo must include the trailer by hand.
