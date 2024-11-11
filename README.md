# Bachelor_Thesis_Montserrat
Study range-shifts in Montserrat forests birds using dynamic occupancy models in 'unmarked'. 

Here, you can find all scripts used for this study, the data is available upon request: 

1. DataExploration.R to explore data and check for any mistakes.
2. SPECIES_Montserrat_bird_counts_analysis.R (where SPECIES equals the 4 letter code of one of the analysed forest bird species) to finally prepare the data for each species, fit the models, base inference on the best model chosen by AIC and afterwards export all data needed for visualisation script pretty_data.R
3. GOF_test_as_loop.R loads the best model objects from local files which have been exported after the modelling workflow for each species (see 2.) to run computation and time intensive MacKenzie-Bailey Goodness-of-Fit test during night in a loop and save the output.
4. pretty_data.R loads all output and prepares or arranges figures as well as produces tables in gt as html and picture for for the manuscript.
