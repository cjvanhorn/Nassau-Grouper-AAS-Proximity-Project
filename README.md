# Nassau Grouper AAS Proximity Project

# Purpose
The purpose of this repository is to provide the finalized data, Bayesian model framework, and code to produce figures presented in _Hydrophone placement yields high variability in detection of_ Epinephelus striatus _calls at a spawning site_, to be published in _Ecological Applications_. 

# Organization
The R script file 'model_and_figures.R' contains the Bayesian model and all code to construct each figure presented in the manuscript. The R data file 'data_for_figures_and_model.RData' contains all finalized data formatted for use in the script file and for figure generation. This includes:
  - **anthro_data** - a binary categorization of minutes in which significant interference is present ('INTERFERENCE_PRESENT') in the acoustic data recorded at a given hydrophone. The data is organized long with each hydrophone's dataset stacked on top of each other and labeled by a factored column termed 'HYDROPHONE'. These are time stamped in three ways:

    (1) DATE_TIME: a POSIXct class column recording the day and time (to the resolution of minute of hour);

    (2) TIME: a character class column recording only the time of day (to the resolution of minute of hour); and

    (3) HOUR: a factor class column of 24 levels recording only the hour of the day.
  - **array_data** - the primary dataset containing all acoustic and visual data at each hydrophone, with columns marked to indicate which hydrophone its data represents. For each hydrophone, 5 forms of data are included:

    (1) HX_AAS_COUNT: the aggregated count of AAS ('HX_AAS_COUNT', where X reflects the hydrophone station number) recorded for each hour (labeled by the 'HOUR_OF_DAY' column) originally extracted by the software FADAR for H1, H5, and H6 and by a human on a 1 minute for every 5 minute duty cycle for H2 and H4;

    (2) HX_FISH_PROX: the categorized proximity of fish to each hydrophone where 0 = not present, 1 = present within 100m of the hydrophone, and 2 = present within 20m of the hydrophone;

    (3) HX_MINUTES_OBSERVED: the number of observable minutes each hour, an observation necessary due to the inequal method of AAS count extraction across hydrophones and the scattered presence of interference in the acoustic data across hydrophones that masked potentially present AAS;

    (4) HX_AAS_PER_MIN: the caculated rate of AAS per minute within each hour, derived simply by dividing 'HX_AAS_COUNT' by 'HX_MINUTES_OBSERVED'; and

    (5) HX_NORMALIZED_AAS: the effort normalized rate of AAS detections per hour derived sipmly by taking 'HX_AAS_PER_MIN' and multiplying by 60 to reflect the estimated number of AAS per hour normalized by the effort possible in observable minutes within each hour.

    Included for reference in this data are temporal data divided in five ways:

    (1) DATE_TIME: a POSIXct class column recording the day and time (to the resolution of hour of day)

    (2) DAY_OF_MONTH: the day of the month the data was recorded in

    (3) HOUR_OF_DAY: the hour of the day the data was recorded in

    (4) DAFM: the number of days before or after the full moon

    (5) DAFS: the number of days before or after divers first observed Nassau Grouper spawning
  - **hydro_coordinates** - a simple data frame including the longitude and latitude coordinates for each hydrophone.
  - **hydro_distance_matrix** - a simple matrix of the distance between each FUNCTIONING HYDROPHONE (hence: H1, H2, H4, H5, H6). Each column from left to right reflects hydrophones ordered from lowest station number (or southernmost station) to highest station number (or northernmost station). Each row reflects this as well, organized from top to bottom (south to north). The distances in this matrix are in meters, and were calculated using coordinates presented in 'hydro_coordinates' and the 'geosphere' package in R.
  - **landmark_coordinates** - a similar data frame to hydro_coordinates only used in the creation of Figure S1A and B. Included are longitude and latitude coordinates for hydrophones and other landmarks used to estimate the position of the bulk of Nassau Grouper around Little Cayman's west end.
  - **lcbathy, worldLLfull, and carib_bathy** - shapefiles used to generate maps of the Caribbean and bathymetry around Little Cayman.
  - **model_formatted_data** - array_data pivoted longer with two new columns to differentiate hydrophone stations ('HYDRO_STATION') and hydrophone brand ('HYDRO_TYPE'). However, the dataset only contains observations up until February 15 due to a malfunction in recording at ST4.  
  - **visual_acoustic_data** - array_data pivoted longer, only including acoustic data (AAS_NORMALIZED_RATES) and bulk Nassau Grouper proximity to a hydrophone ('FISH_PROX').

# Important Notes
The third hydrophone station (LS3) malfunctioned early during recording and did not collect acoustic data. ST4 also experienced an abrupt malfunction near midnight of February 15, thus it did not collect data for the remainder of the recording period. In array_data, all acoustic data for ST4 is NA upon February 15 and after. However, model_formatted_data only includes data up until February 15 due to the necessity of the model to have consistent observations across all stations. 

Our model runs using the rethinking package developed by McElreath (2020) in R. This package, while remarkable, requires significant preparation in order to load properly. Details into this preparation can be found commented in the code provided. Any questions that arise from loading this package should not be directed here.
