Python code for the evaluation of total emissions of CH4 and CO (as well as other chemical species after some changes) from a netcdf emissions database (e.g. CH4 EDGAR https://edgar.jrc.ec.europa.eu/gallery?release=v60ghg&substance=CH4&sector=TOTALS). The spatial selection is performed over provinces that are selected from the naturalhearth admin1.geojson file (https://www.naturalearthdata.com/downloads/10m-cultural-vectors/).

Only for Italy, also the evaluation of total emissions from the ISPRA database (https://www.isprambiente.gov.it/it) can be performed.

The code need to access to the admin1.geojson file from the naturalhearth database as well as to the netcdf emission database from EDGAR. Both the geojson file and the EDGAR databases must be stored locally.

The code is divided in a main program (select_emissions.py) and other modules that contain specific functions for the data elaboration.

The file **select_emissions_param.py** is used to specify the region over wich the evaluation of total emissions must be performed, as well as other parameters.
The provinces over which to perform the evaluation are defined using a list:
    prov_numname = [(3, 'Novara'),(12, 'Varese'),(13, 'Como')]
    For each province a touple with a number and a string is required. The number correspond to the province number in the ISPRA database (thus the province number is required only for italian provinces in the case of the evaluation of ISPRA emissions), while the string correspond to the province name on the admin1.geojson database. the case reported above will then select emissions from the provinces of Novara, Varese and Como. For non-italian provinces a random number can be used as province number. NB1 in case a province is divided in more parts (e.g. some islands are included in the province boundary) only the biggest part of the province is retained. NB2 provinces must be adjacent (e.g the evaluation of the emissions from the provinces of Paris and Beijing is not possible).
The evaluation of ISPRA emission is performed by using the IPR variable. If True, the evaluation of ISPRA emissions is performed, if False not.
Other parameters in the select_emissions_param.py file are used only to define filenames and plot titles.

The script can perform different tasks: 
1) evaluate yearly total emissions over some provinces with their relative error. An algorithm for the error evaluation was defined: the error is evaluated by selecting an inner and an outer boundaries and by evaluating the differences in the emission estimation for the "standard" case and the inner/outer ones. The error that is obtained with this routine is than compared with the EDGAR relative error and the larger error is retained. 
2) evaluate monthly emissions. This is obtained by multiplying the yearly EDGAR emissions by the seasonal profile that is provided by EDGAR. A profile is defined for each "region" (see EDGAR database for more details). The region profile must be specified with a number (profile_region) in the select_emissions_param file.
3) estimation of the emissions in the years after the last EDGAR netcdf database. Since EDGAR databases for CO and CH4 are up to 2015 and 2018 respectively, a routine for the estimation of emissions also in the following years was developed. Multiple polynomial fits at different orders are performed in order to evaluate future emissions and their error.

Results for total emissions are then stored in tabular files. Some plots that show the total emissions versus time as well as a 2D map reporting the bins that have been used for the emissions estimation are produced.
