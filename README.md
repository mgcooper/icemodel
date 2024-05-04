# IceModel

**A 1-d radiative-thermodynamic heat transfer model that simulates meltwater runoff from glacier ice.**

- [IceModel](#icemodel)
  - [Background](#background)
  - [Getting Started](#getting-started)
  - [Advanced Use](#advanced-use)
  - [Input Data](#input-data)
  - [Summary](#summary)
  - [Naming Conventions](#naming-conventions)
    - [1. met files](#1-met-files)
    - [2. user data](#2-user-data)
  - [References](#references)

## Background

`IceModel` is based on models of glacier ice and snowpack described in Liston et al (1999) and Jordan (1991) (*SNTHERM*). State variables are ice temperature and liquid water content, from which hydrologic and thermodynamic quantities are computed, including meltwater runoff.

Following Liston et al (1999), `IceModel` uses the two-stream radiative transfer model described in Schlatter (1972), updated with spectral detail following Brandt and Warren (1993). The surface energy balance largely follows Liston et al (1999), including the turbulent flux parameterization based on Monin-Obukhov similarity theory (Monin & Obukhov, 1954). `IceModel` draws on both Liston et al (1999) and *SNTHERM* for subsurface thermodynamics, including the vapor diffusion parameterization, but implements an enthalpy-conserving phase change numerical scheme described in Clark et al. (2021) and Swaminathan and Voller (1993).

A brief technical description is available [here](https://github.com/mgcooper/icemodel/blob/main/docs/IcemodelTechnicalDoc.pdf).

<!-- `icemodel` solves the unsteady one-dimensional heat equation:

$$
\dfrac{\partial H}{\partial t} = \dfrac{\partial}{\partial z}\biggr(F\biggr) + S
$$

where *H* [J m$^{-3}$] is enthalpy, *t* [s] is time, *F* [W m$^{-2}$] is net heat flux along vertical dimension *z* [m], and *S* [W m$^{-3}$] is a source/sink term. -->

## Getting Started

Thanks for your interest. To get started, here's what we recommend.

- The main program is `icemodel/icemodel.m`. Open the function to get a sense for the model structure.
- Open and run Example 1 in `demo/demo.m`. This will run an `IceModel` simulation for the KAN_M weather station, located on the Greenland ice sheet, for year 2016 on a 1-hr timestep.
- Set `saveflag=true` and re-run Example 1. Notice how the `demo/output` directory is created, and the model output is saved there.
- Set `backupflag=true` and re-run Example 1. Notice how the files are backed-up. By default, saveflag and backupflag are both false.

## Advanced Use

To set run-specific model options and parameters, open and edit the function `icemodel/+icemodel/setopts.m`. In your matlab terminal:

- Type `edit icemodel.setopts` and press enter. Set the options and resave the function. Then run `icemodel.run.point` with the new options (see `demo.m` for an example of how to call `icemodel.run.point`).

For extended use, specify the model input and output directories by editing the configuration function `icemodel/+icemodel/config.m`. By default, these directories are set to the top-level folders `input/` and `output/`. Note that the `.gitignore` file that ships with this repo ignores these folders.

<!-- There are options in that script to save the model output and evaluate the results against observations. -->

<!-- , such as which input forcing dataset to use, and whether to run `icemodel` or `skinmodel` -->

## Input Data

Example required inputs are in `demo/inputs`. These include the meteorological forcing data in `inputs/met`, user data in `inputs/userdata`, and the inputs to the two-stream spectral model in `inputs/spectral`.

The `spectral` directory contains values for the absorption coefficient of pure ice from Warren et al. 2008, in-situ absorption coefficients for glacier ice from Cooper et al. 2021, a downwelling solar spectrum for the Arctic atmosphere generated with `ATRAN` (Lord, 1991), and a library of mie-scattering coefficients as described in Cooper et al. 2021. 

The `userdata` directory contains alternative model forcing data that is "swapped out" with the standard model forcing to test hypotheses about processes and model sensitivity to forcings. For example, users can select input forcing data generated from climate model output by setting the `userdata` and `uservars` configuration parameters (see `demo/demo.m`), or swap weather station albedo with MODIS albedo. Setting `userdata="modis"` and `uservars="albedo"` would replace the albedo values in the input meteorological forcing file with modis albedo for the same time and location.

## Summary

1. Install this repo and place it on your Matlab path
2. Edit and run the following files

| Function Name | Description | How to Run |
| --- | --- | --- |
| `icemodel.config` | Set global configuration settings. | Type `edit icemodel.config` then press enter. Edit the environment variables and save the function. |
| `icemodel.setopts` | Set run-specific model configuration. | Type `edit icemodel.setopts` then press enter. Edit the model options and save the function. |
| `icemodel.run.point` | Run the model at a point. | See the example in `demo.m`, and the function arguments in `icemodel.run.point` for additional configuration. |
| `demo.m` | Script to run and evaluate the model output. | Place this repo on your matlab path, edit and run the script. |

## Naming Conventions

### 1. met files

The met (forcing) file naming convention is:
`met_SITENAME_FORCINGS_YYYY_TIMESTEP`

For example:

`met_KANM_KANM_2016_1hr.mat` = met (forcing) data for site KAN-M with KAN-M forcings for year 2016 at a 1 hour timestep.

`met_KANM_KANM_2016_15m.mat` = met (forcing) data for site KAN-M with KAN-M forcings for year 2016 at a 15 minute timestep.

`met_KANM_MERRA_2016_15m.mat` = met (forcing) data for site KAN-M with MERRA-2 forcings for year 2016 at a 15 minute timestep.

### 2. user data

The "userdata" file naming convention is:
`FORCINGS_SITENAME_YYYY`

For example:

`MERRA_KANM_2016.mat` = forcing data from the MERRA-2 climate model for the KAN-M weather station location for year 2016.

<!-- 
The function `METINIT.m` will then swap out the KAN-M albedo data in the met forcing data with the modis albedo.

In practice, however, these options are set programmatically by passing the `sitename`, `forcingdata`, `userdata`, `uservars`, `meltmodel`, `startyear`, and `endyear` variables to the `icemodel_opts.m` function, which builds the `opts.metfname` string. The function `icemodel_run.m` is used to set these variables, pass them to `icemodel_opts.m`, and then pass `opts` to the main program `icemodel.m`. -->

<!-- ## Example model configuration

```matlab
    % set the main configuration options
    sitename    = 'KANM';
    startyear   = 2018;
    endyear     = 2018;
    meltmodel   = 'icemodel';
    forcingdata = 'KANM';
    userdata    = 'modis';
    uservars    = 'albedo';

    % build the 'opts' model configuration structure
    opts = icemodel_opts(   sitename,meltmodel,forcingdata,userdata,uservars,  ...
                            startyear,endyear);

    % run the model
    [ice1,ice2,met,opts] = icemodel(opts);
``` -->

<!-- ## Code reference -->

<!-- will insert zenodo ref here -->
<!-- Cooper, M G (2022). `icemodel` [code](https://doi.org/some-doi-number) -->

<!-- ## Data reference

### Input data

| Data | Source| Download website | Usage |
|-------|---------|-----------------|-----|
| Climate data | MAR3.11 |  | model forcing |
| Climate reanalysis | MERRA2 |  | model forcing | -->

<!-- ### Output data

Cooper, M G (2022). icemodel output [Data set](https://doi.org/some-doi-number)

| Data | Format| Content | Usage |
|-------|---------|-----------------|-----|
| Runoff | netcdf | runoff | model output | -->

<!-- ## Contributing modeling software

| Model | Version | Repository Link | DOI |
|-------|---------|-----------------|-----|
| ATRAN | version | link to code repository | link to DOI dataset | -->

<!-- ## Reference

Cooper, M G, et al. *Greenland Ice Sheet runoff reduced by meltwater refreezing in bare ice*, in revision ([preprint](https://assets.researchsquare.com/files/rs-842710/v1_covered.pdf?c=1631877020)) -->

## References

Anderson E A 1976 A  point energy and mass balance model of a snow cover ed N W S United States [Online](https://repository.library.noaa.gov/view/noaa/6392)

Clark M P, Zolfaghari R, Green K R, Trim S, Knoben W J M, Bennett A, Nijssen B, Ireson A and Spiteri R J 2021 The Numerical Implementation of Land Models: Problem Formulation and Laugh Tests Journal of Hydrometeorology 22 1627–48 [Online](https://journals.ametsoc.org/view/journals/hydr/22/6/JHM-D-20-0175.1.xml)

Cooper M G, Smith L C, Rennermalm Å K, Tedesco M, Muthyala R, Leidman S Z, Moustafa S E and Fayne J V 2021 *Spectral attenuation coefficients from measurements of light transmission in bare ice on the Greenland Ice Sheet* The Cryosphere 15 1931–53, [Online](https://tc.copernicus.org/articles/15/1931/2021/)

Jordan R 1991 *A One-Dimensional Temperature Model For a Snowpack* (Hanover, NH: Cold Regions Research and Engineering Laboratory) [Online](http://hdl.handle.net/11681/11677)

Liston G E, Bruland O, Elvehøy H and Sand K 1999 *Below-surface ice melt on the coastal Antarctic ice sheet* Journal of Glaciology 45 273–85, [Online](https://doi.org/10.3189/002214399793377130)

Lord S D 1992 A new software tool for computing Earth’s atmospheric transmission of near-and far-infrared radiation vol 103957 (Ames Research Center)

Monin A S and Obukhov A M 1954 Basic laws of turbulent mixing in the surface layer of the atmosphere Contrib. Geophys. Inst. Acad. Sci. USSR 151 e187

Schlatter T W 1972 *The Local Surface Energy Balance and Subsurface Temperature Regime in Antarctica* J. Appl. Meteor. 11 1048–62 [Online](http://journals.ametsoc.org/doi/abs/10.1175/1520-0450%281972%29011%3C1048%3ATLSEBA%3E2.0.CO%3B2)

Warren S G and Brandt R E 2008 *Optical constants of ice from the ultraviolet to the microwave: A revised compilation* J. Geophys. Res. 113 D14220, [Online](http://onlinelibrary.wiley.com/doi/10.1029/2007JD009744/abstract)
