# IceModel

**A 1-d radiative-thermodynamic heat transfer model for glacier ice.**

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=mgcooper/icemodel&file=demo/demo.m) [![DOI](https://zenodo.org/badge/534812039.svg)](https://zenodo.org/doi/10.5281/zenodo.11539329) ![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/mgcooper/icemodel?include_prereleases) [![GitHub license](https://img.shields.io/github/license/mgcooper/icemodel)](https://github.com/mgcooper/icemodel/blob/main/LICENSE)

- [IceModel](#icemodel)
  - [Background](#background)
  - [Getting Started](#getting-started)
  - [Advanced Use](#advanced-use)
    - [Global configuration: Specify workspace paths](#global-configuration-specify-workspace-paths)
    - [Runtime configuration: Specify model options](#runtime-configuration-specify-model-options)
  - [Input Data](#input-data)
  - [User Data](#user-data)
  - [Summary](#summary)
  - [Naming Conventions](#naming-conventions)
    - [1. General conventions](#1-general-conventions)
    - [2. Met files](#2-met-files)
    - [3. User data files](#3-user-data-files)
    - [4. Output files](#4-output-files)
  - [References](#references)
  - [System Requirements](#system-requirements)
  - [Installation Guide](#installation-guide)
  - [Contribute](#contribute)
  - [How do I cite this?](#how-do-i-cite-this)

## Background

`IceModel` is based on models of glacier ice and snowpack described in Liston et al (1999) and Jordan (1991) (*SNTHERM*). State variables are ice temperature and liquid water content, from which various hydrologic and thermodynamic quantities are computed, including the flux of meltwater runoff.

Following Liston et al (1999), `IceModel` implements the two-stream radiative transfer model described in Schlatter (1972), updated with spectral detail following Brandt and Warren (1993). The surface energy balance largely follows Liston et al (1999), including the turbulent flux parameterization based on Monin-Obukhov similarity theory (Monin & Obukhov, 1954). `IceModel` follows both Liston et al (1999) and *SNTHERM* for subsurface thermodynamics, including representations of solar radiative heating, conductive heat transfer, and vapor diffusion, but implements an enthalpy-conserving numerical phase change scheme described in Clark et al. (2021) and Swaminathan and Voller (1993).

A technical description is available [here](https://github.com/mgcooper/icemodel/blob/main/doc/IcemodelTechnicalDoc.pdf).

<!-- `icemodel` solves the unsteady one-dimensional heat equation:

$$
\dfrac{\partial H}{\partial t} = \dfrac{\partial}{\partial z}\biggr(F\biggr) + S
$$

where *H* [J m$^{-3}$] is enthalpy, *t* [s] is time, *F* [W m$^{-2}$] is net heat flux along vertical dimension *z* [m], and *S* [W m$^{-3}$] is a source/sink term. -->

## Getting Started

Thanks for your interest. To get started, here's what we recommend:

- Check the [system requirements](#system-requirements) and [installation guide](#installation-guide).
- If you do not have a MATLAB license, you can run this software using a free MATLAB Online account: [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=mgcooper/icemodel&file=demo/demo.m)
- The main program is `icemodel/icemodel.m`. Open the function to get a sense for the model structure.
- Open and run Example 1 in `demo/demo.m`. This will run an `IceModel` simulation for the KAN_M weather station, located on the Greenland ice sheet, for year 2016 on a 1-hr timestep.
- Inspect the demo plot created by the call to `icemodel.plot.enbal`. The simulated energy fluxes should closely track the weather station values.
- Set `saveflag=true` and re-run Example 1. Notice how the `demo/output` directory is created, and the model output is saved there.
- Set `backupflag=true` and re-run Example 1. Notice how the files are backed-up. By default, saveflag and backupflag are both false.

The examples in `demo.m` run IceModel in its "SkinModel" surface energy balance configuration. Run-time will depend on your computer, but should take less than one minute. An IceModel configuration (which includes a full subsurface energy balance) should take between one and a few minutes to run. Initial run times may be longer due to JIT compilation.

## Advanced Use

### Global configuration: Specify workspace paths

The model input and output directories default to the top-level folders `input/` and `output/` (note that the `.gitignore` in this repo ignores these folders).

To specify custom input and output directories, use the configuration function `icemodel/+icemodel/config.m`. In your matlab terminal:

- Type `edit icemodel.config` and press enter.
- Read the detailed documentation to understand the model input and output directory structure, and how to set them programmatically.

<!-- Note that in `demo.m` the `casename` argument is passed to the configuration function: `cfg = icemodel.config(casename="demo")`. This sets the input and output folders to `demo/input` and `demo/output`. The `casename` argument to `icemodel.config` currently does not serve any other purpose. -->

### Runtime configuration: Specify model options

To set run-specific model options and parameters, open and edit the function `icemodel/+icemodel/setopts.m`. In your matlab terminal:

- Type `edit icemodel.setopts` and press enter.
- Edit the options and resave the function.
- Run the model with the new options (see `demo.m` for an example of how to call `icemodel.run.point` to run the model).

## Input Data

Example input files are in `demo/inputs`. These include the meteorological forcing data in `inputs/met`, inputs to the two-stream spectral model in `inputs/spectral`, and optional "user data" in `inputs/userdata`.

The `spectral` directory contains values for the absorption coefficient of pure ice from Warren et al. 2008, in-situ absorption coefficients for glacier ice from Cooper et al. 2021, a downwelling solar spectrum for the Arctic atmosphere generated with `ATRAN` (Lord, 1991), and a library of mie-scattering coefficients as described in Cooper et al. 2021.

## User Data

The `inputs/userdata` directory contains alternative model forcings that can be "swapped out" with the standard model forcings to test hypotheses about processes and model sensitivity. For instance, users can prepare input forcing data generated from observations, climate model output, or satellite remote sensing, and place these files in `userdata`. Unlike the forcing files in `inputs/met`, these files do not need to contain the complete set of model forcings.

To swap out a variable in the input met file with a variable in a userdata file, set the `userdata` and `uservars` configuration parameters (see `demo/demo.m`). For example, setting `userdata="modis"` and `uservars="albedo"` would replace the albedo values in the input meteorological forcing file with modis albedo for the same time and location. This would require placing a file named `<sitename>_modis_<year>` (see [Naming Conventions](#naming-conventions)) in the `userdata` directory, containing a timetable named `Data` with a variable (column) named `albedo`.

## Summary

1. Install this repo and place it on your Matlab path
2. Open and run `demo/demo.m`
3. For advanced use, edit and run the following files:

| Function Name | Description | How to Run |
| --- | --- | --- |
| `icemodel.config` | Set global configuration (model input and output paths). | Type `edit icemodel.config` then press enter. Set the environment variables programmatically as needed. |
| `icemodel.setopts` | Set run-specific model configuration. | Type `edit icemodel.setopts` then press enter. Edit the model options and save the function. |
| `icemodel.run.point` | Run the model at a point. | See the example in `demo.m`, and the function arguments in `icemodel.run.point` for additional configuration. |
| `demo.m` | Script to run and evaluate the model output. | Place this repo on your matlab path, edit and run the script. |

## Naming Conventions

### 1. General conventions

- File names are in lowercase for operating system compatibility.
- Non-standard characters in filename parts (e.g., "-", "&") are replaced with blanks. For example, sitename "KAN-M" becomes "kanm".

<!-- - Version 7.3 MAT-files, which are based on the HDF5 format, are used for input and (by default) for output. For advanced usage -->

### 2. Met files

The met (forcing data) file naming convention is:

- `met_SITENAME_FORCINGS_YYYY_TIMESTEP`

Examples:

- `met_kanm_kanm_2016_1hr.mat` specifies a met (forcing) data file for site KAN-M with KAN-M forcings for year 2016 at a 1-hour timestep.
- `met_kanm_kanm_2016_15m.mat` specifies a met (forcing) data file for site KAN-M with KAN-M forcings for year 2016 at a 15-minute timestep.
- `met_kanm_merra_2016_15m.mat` specifies a met (forcing) data file for site KAN-M with MERRA-2 forcings for year 2016 at a 15-minute timestep.

Each met file must contain a timetable object named `met` with one column for each forcing variable. See the example met file.

### 3. User data files

The "userdata" file naming convention is:

- `SITENAME_FORCINGS_YYYY`

Note: at this time, hourly userdata files are supported, thus unlike the met file naming convention, there is no `TIMESTEP` file part.

Examples:

- `kanm_merra_2016.mat` specifies a user data file with MERRA-2 climate model forcings for the KAN-M weather station location for year 2016 on a 1-hr timestep.
- `kanm_modis_2016.mat` specifies a user data file with MODIS satellite albedo values for the KAN-M weather station location for year 2016 on a 1-hr timestep.

Each userdata file must contain a timetable named `Data` with column names matching the met file column-naming conventions. See the example met file in `demo/input/`.

### 4. Output files

IceModel produces two outputs: `ice1` and `ice2`, representing 1-dimensional and 2-dimensional data, respectively, with dimensions time and depth. The `ice1` data are saved in a `timetable` object and `ice2` are saved in a `struct` object. An IceModel simulation with `opts.saveflag` set true will save these two objects as well as the `opts` struct according to the following convention:

```sh
ICEMODEL_OUTPUT_PATH/SITENAME/SMBMODEL/YYYY/ice1_FORCINGS_forcings_USERDATA_USERVARS.mat
ICEMODEL_OUTPUT_PATH/SITENAME/SMBMODEL/YYYY/ice2_FORCINGS_forcings_USERDATA_USERVARS.mat
ICEMODEL_OUTPUT_PATH/SITENAME/SMBMODEL/opts/opts_FORCINGS_forcings_USERDATA_USERVARS.mat
```

Here, `ICEMODEL_OUTPUT_PATH` is an environment variable set by the `icemodel.config` function, the lowercase "forcings" is a string literal used to join the `FORCINGS` and `USERDATA` string variables, and `SITENAME`, `SMBMODEL`, `FORCINGS`, `USERDATA`, and `USERVARS` are parameters passed to the `icemodel.setopts` function (either directly or indirectly via the helper function `icemodel.run.point`). One `YYYY` folder is created for each year in the `SIMYEARS` parameter passed to `icemodel.setopts`.

Note that the `ICEMODEL_OUTPUT_PATH/SITENAME/SMBMODEL/YYYY` subfolders are generated automatically in the `icemodel.setopts` function if they do not exist.

Example: An IceModel simulation for the KAN-M weather station location for years 2015:2016 using MERRA forcings, with `userdata='modis'` and `uservars='albedo'` will produce the following output files:

```sh
ICEMODEL_OUTPUT_PATH/kanm/icemodel/2015/ice1_merra_forcings_modis_albedo.mat
ICEMODEL_OUTPUT_PATH/kanm/icemodel/2015/ice2_merra_forcings_modis_albedo.mat
ICEMODEL_OUTPUT_PATH/kanm/icemodel/2016/ice1_merra_forcings_modis_albedo.mat
ICEMODEL_OUTPUT_PATH/kanm/icemodel/2016/ice2_merra_forcings_modis_albedo.mat
ICEMODEL_OUTPUT_PATH/kanm/icemodel/opts/opts_merra_forcings_modis_albedo.mat
```

If `userdata` and `uservars` are not supplied to `icemodel.setopts` or are supplied as empty arrays `[]`, or empty chars `''` (they are optional arguments but must be supplied as empty if subsequent arguments `saveflag` and/or `backupflag` are supplied), the filenames would be:

```sh
ICEMODEL_OUTPUT_PATH/kanm/icemodel/2015/ice1_merra_forcings_merra_albedo.mat
ICEMODEL_OUTPUT_PATH/kanm/icemodel/2015/ice2_merra_forcings_merra_albedo.mat
```

Here, the `userdata` parameter takes on the value of the `forcings` parameter (`'merra'`), whereas the `uservars` parameter takes the value `'albedo'`, which is the default value set in `icemodel.setopts`. This reflects the focus of IceModel on studying the sensitivity of SMB models to ice albedo.

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

Cooper, M G, et al. *Greenland Ice Sheet runoff reduced by meltwater refreezing in bare ice* -->

## References

<!-- Anderson E A 1976 *A point energy and mass balance model of a snow cover* ed N W S United States [Online](https://repository.library.noaa.gov/view/noaa/6392) -->

Brandt R E and Warren S G 1993 *Solar-heating rates and temperature profiles in Antarctic snow and ice* Journal of Glaciology 39 99–110

Clark M P, Zolfaghari R, Green K R, Trim S, Knoben W J M, Bennett A, Nijssen B, Ireson A and Spiteri R J 2021 *The Numerical Implementation of Land Models: Problem Formulation and Laugh Tests* Journal of Hydrometeorology 22 1627–48 [Online](https://journals.ametsoc.org/view/journals/hydr/22/6/JHM-D-20-0175.1.xml)

Cooper M G, Smith L C, Rennermalm Å K, Tedesco M, Muthyala R, Leidman S Z, Moustafa S E and Fayne J V 2021 *Spectral attenuation coefficients from measurements of light transmission in bare ice on the Greenland Ice Sheet* The Cryosphere 15 1931–53, [Online](https://tc.copernicus.org/articles/15/1931/2021/)

Jordan R 1991 *A One-Dimensional Temperature Model For a Snowpack* (Hanover, NH: Cold Regions Research and Engineering Laboratory) [Online](http://hdl.handle.net/11681/11677)

Liston G E, Bruland O, Elvehøy H and Sand K 1999 *Below-surface ice melt on the coastal Antarctic ice sheet* Journal of Glaciology 45 273–85, [Online](https://doi.org/10.3189/002214399793377130)

Lord S D 1992 *A new software tool for computing Earth’s atmospheric transmission of near-and far-infrared radiation* vol 103957 (Ames Research Center)

Monin A S and Obukhov A M 1954 *Basic laws of turbulent mixing in the surface layer of the atmosphere* Contrib. Geophys. Inst. Acad. Sci. USSR 151 e187

Schlatter T W 1972 *The Local Surface Energy Balance and Subsurface Temperature Regime in Antarctica* J. Appl. Meteor. 11 1048–62 [Online](http://journals.ametsoc.org/doi/abs/10.1175/1520-0450%281972%29011%3C1048%3ATLSEBA%3E2.0.CO%3B2)

Swaminathan C R and Voller V R 1993 *ON THE ENTHALPY METHOD* International Journal of Numerical Methods for Heat & Fluid Flow 3 233–44

Warren S G and Brandt R E 2008 *Optical constants of ice from the ultraviolet to the microwave: A revised compilation* J. Geophys. Res. 113 D14220, [Online](http://onlinelibrary.wiley.com/doi/10.1029/2007JD009744/abstract)

## System Requirements

- Requires MATLAB&reg; version >=9.2* (R2017a).
  - R2016b is a hard limit on compatibility: IceModel forcing files are currently stored as `timetable` objects which were introduced (along with `string`) in R2016b.
- Developed and tested on MacOS Sonoma (Intel silicon), using MATLAB R2022b.
- Runs in MATLAB Online (a linux-based system), tested on R2022b.
- Runs on Windows 10, tested on R2017a.
<!-- - Runs in Octave on MacOS Sonoma (Intel silicon), using Octave version 9.2. -->

*Note that the main program `icemodel/icemodel.m` and core IceModel functions (UPPERCASE filenames in `icemodel/`), are written in a minimalist matlab style: all functions are compatible with code generation, all numerical methods employ custom hand-written solvers, and there are no toolbox dependencies or modern matlab conveniences such as `arguments` input parsers.

Exceptions to this style include namespace functions (e.g. `icemodel/+icemodel`), which by convention are helper functions not required by the numerical model. The `demo.m` script uses a modern matlab approach including name=value syntax (requires >=R2021a). The `arguments` parser used in `icemodel.run.point` requires >=R2019b. For users running pre-R2019b, see `demo/demo_pre_R2019.m`.

If you encounter incompatibilities, please [open an issue](https://github.com/mgcooper/icemodel/issues).

<!-- datetime, R2014b
string, R2016b
timetable, R2016b
isfolder / isfile, R2017b
isStringScalar, R2017b
convertStringsToChars, R2017b
arguments, R2019b
renamevars, R2020a
name=value, R2021a -->

## Installation Guide

Download this repo and place it on the matlab path.

In a terminal:

```sh
git clone https://github.com/mgcooper/icemodel.git
```

In your matlab terminal:

```matlab
cd('/path/to/this/repo')
setup()
```

Installation should only take a few seconds. If you encounter any issues, please [open an issue](https://github.com/mgcooper/icemodel/issues).

## Contribute

If you find a bug, have a question, or want to contribute, feel free to open an [issue](https://github.com/mgcooper/icemodel/issues) or start a [discussion](https://github.com/mgcooper/icemodel/discussions).

## How do I cite this?

If you find this software useful, please consider citing the software release in `CITATION.cff` (click the link in the "About" section along the right sidebar of this repo).
