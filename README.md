# icemodel

**A 1-d radiative-thermodynamic heat transfer model that simulates meltwater runoff from glacier ice.**

## Background
`icemodel` is loosely based on models described in Liston et al. (1999) and Jordan (1991) (*SNTHERM*). `icemodel` follows Liston et al. (1999) in adopting the two-stream radiative transfer model described in Schlatter (1972), the water vapor diffusion parameterization described in Anderson (1976), and a turbulent flux parameterization based on Monin-Obukhov similarity theory (Monin & Obukhov, 1954). `icemodel` loosely adopts *SNTHERM* thermodynamics, including the melt-zone preconditioner, but implements a numerical enthalpy conservation scheme described in Clark et al. (2021). State variables are ice temperature and liquid water content, from which various hydrologic and thermodynamic quantities are computed, namely meltwater runoff.

## Usage

The main program is `icemodel.m`, which calls functions saved in `functions/`. To get started, open the file `drive/icemodel_config.m` and set the global configuration options (environment variables). Then open `drive/icemodel_opts.m` which is a function that sets run-specific model options and parameters, such as which input forcing dataset you want to use, and whether you want to run `icemodel` or `skinmodel`. Then run the model by calling `drive/icemodel_run.m`. There are options in that script to save the model output and evaluate the results against observations.

The required inputs are in `inputs/`. These include the meteorological forcing data in `inputs/met/`, user data in `inputs/userdata/`, and the inputs to the two-stream model in `inputs/spectral/`. The spectral data inputs include values for the absorption coefficient of pure ice from Warren et al. 2008, in-situ absorption coefficients for glacier ice from Cooper et al. 2021, a proto-typical downwelling solar spectrum for the Arctic atmosphere generated with `ATRAN` (Lord, 1991), and a library of mie-scattering coefficients as described in Cooper et al. 2021. The "user data" directory contains alternative model forcing data that is "swapped out" with the standard model forcing to test hypotheses about processes and determine model sensitivity to forcings. For example, the user can select input forcing data generated from climate model output by setting the `opts.metfname` configuration parameter. This parameter points to the input meteorological forcing file in `inputs/met/`. The user can then override one or more of the forcing variables in this input file by setting `opts.userdata` and `opts.uservars` to point to an alternative forcing dataset. For example, `opts.userdata=modis` and `opts.uservars=albedo` would replace the albedo values in the input meteorological forcing file with modis albedo for the same time and location.

In summary:
1. Install this repo and place it on your Matlab path
2. Download the input data (model forcing) from [Input data](#input-data)
3. Run the following scripts in the `drive` directory to run the model:

| Script Name | Description | How to Run |
| --- | --- | --- |
| `icemodel_config.m` | global configuration settings | place this repo on your matlab path, edit and run the script |
| `icemodel_opts.m` | run-specific model configuration | place this repo on your matlab path, edit and run the script |
| `icemodel_run.m` | Script to run the model | place this repo on your matlab path, edit and run the script |
| `icemodel_eval.m` | Script to evaluate model output | place this repo on your matlab path, edit and run the script |

## Model input

### 1. met files
The met files are named using the following protocol:
`met_SITENAME_FORCINGDATA_YYYY_TIMESTEP`

For example:

`met_ak4_MAR_2009_1hr.mat` = met data for site ak4 from MAR forcings for year 2009 at a 1 hour timestep.

`met_KANM_KANM_2009_15m.mat` = met data for site KAN-M from KAN-M forcings for year 2009 at a 15 minute timestep.

`met_KANM_MAR_2009_15m.mat` = met data for site KAN-M from MAR forcings for year 2009 at a 15 minute timestep.

### 2. user data
The "userdata" files are named using the following protocol:
`FORCINGDATA_SITENAME_YYYY`

For example:

`KANM_behar_2015.mat` = forcing data from the KAN-M weather station for site 'behar' for year 2015.

`modis_behar_2015.mat` = modis albedo data for site 'behar' for year 2015.

Say you want to run a simulation at the KAN-M weather station for year 2018 at a 15-m timesep (15 m is required for icemodel, 1 hr is suitable for skinmodel), but you want to use MODIS albedo instead of KAN-M albedo. You would use the following settings:

    opts.metfname = 'met_KANM_KANM_2018_15m.mat';
    opts.userdata = 'modis';
    opts.uservars = 'albedo';

The function `METINIT.m` will then swap out the KAN-M albedo data in the met forcing data with the modis albedo.

In practice, however, these options are set programmatically by passing the `sitename`, `forcingdata`, `userdata`, `uservars`, `meltmodel`, `startyear`, and `endyear` variables to the `icemodel_opts.m` function, which builds the `opts.metfname` string. The function `icemodel_run.m` is used to set these variables, pass them to `icemodel_opts.m`, and then pass `opts` to the main program `icemodel.m`.

### Example model configuration:

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
    

## Code reference
forthcoming

## Data reference

### Input data

| Data | Source| Download website | Usage |
|-------|---------|-----------------|-----|
| Climate data | MAR3.11 |  | model forcing | 
| Climate reanalysis | MERRA2 |  | model forcing | 

### Output data

Cooper, M G (2022). icemodel output [Data set]. DataHub. https://doi.org/some-doi-number

| Data | Format| Content | Usage |
|-------|---------|-----------------|-----|
| Runoff | netcdf | runoff | model output | 


## Contributing modeling software
| Model | Version | Repository Link | DOI |
|-------|---------|-----------------|-----|
| ATRAN | version | link to code repository | link to DOI dataset |

## Journal reference
Cooper, M G, et al. *Greenland Ice Sheet runoff reduced by meltwater refreezing in bare ice*, in revision (preprint: https://assets.researchsquare.com/files/rs-842710/v1_covered.pdf?c=1631877020)

# References

Anderson E A 1976 A  point energy and mass balance model of a snow cover ed N W S United States Online: https://repository.library.noaa.gov/view/noaa/6392

Clark M P, Zolfaghari R, Green K R, Trim S, Knoben W J M, Bennett A, Nijssen B, Ireson A and Spiteri R J 2021 The Numerical Implementation of Land Models: Problem Formulation and Laugh Tests Journal of Hydrometeorology 22 1627–48 https://journals.ametsoc.org/view/journals/hydr/22/6/JHM-D-20-0175.1.xml

Cooper M G, Smith L C, Rennermalm Å K, Tedesco M, Muthyala R, Leidman S Z, Moustafa S E and Fayne J V 2021 *Spectral attenuation coefficients from measurements of light transmission in bare ice on the Greenland Ice Sheet* The Cryosphere 15 1931–53, https://tc.copernicus.org/articles/15/1931/2021/

Jordan R 1991 *A One-Dimensional Temperature Model For a Snowpack* (Hanover, NH: Cold Regions Research and Engineering Laboratory) http://hdl.handle.net/11681/11677

Liston G E, Bruland O, Elvehøy H and Sand K 1999 *Below-surface ice melt on the coastal Antarctic ice sheet* Journal of Glaciology 45 273–85, https://doi.org/10.3189/002214399793377130

Lord S D 1992 A new software tool for computing Earth’s atmospheric transmission of near-and far-infrared radiation vol 103957 (Ames Research Center)

Monin A S and Obukhov A M 1954 Basic laws of turbulent mixing in the surface layer of the atmosphere Contrib. Geophys. Inst. Acad. Sci. USSR 151 e187

Schlatter T W 1972 *The Local Surface Energy Balance and Subsurface Temperature Regime in Antarctica* J. Appl. Meteor. 11 1048–62 http://journals.ametsoc.org/doi/abs/10.1175/1520-0450%281972%29011%3C1048%3ATLSEBA%3E2.0.CO%3B2

Warren S G and Brandt R E 2008 *Optical constants of ice from the ultraviolet to the microwave: A revised compilation* J. Geophys. Res. 113 D14220, http://onlinelibrary.wiley.com/doi/10.1029/2007JD009744/abstract
