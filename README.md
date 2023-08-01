# Code
There are three folders containing all the codes for my thesis, "Estimation of case fatality rate of the 2014-2016 Ebola Virus Disease outbreak in the presence of missing data: a comparison between complete case analysis, multiple imputation and inverse probability weighting".

Among these, two folders named "Simulation for confirmed cases" and "Simulation for probable cases" are used for simulation model. Within these two folders, there are 2 data files (.csv) and 10 code files (.R) designed for simulations. In details, 

    2 data files (.csv): 
               full_data.csv & original_data.csv;
    1 code file (.R) for data cleaning and preparation: 
               full_data.R;
    9 code files (.R) for simulations:
               10_MAR (representing 10% of outcomes are MAR)
               10_MCAR (representing 10% of outcomes are MCAR)
               10_MNAR (representing 10% of outcomes are MNAR)
               30_MAR (representing 30% of outcomes are MAR)
               30_MCAR (representing 30% of outcomes are MCAR)
               30_MNAR (representing 30% of outcomes are MNAR)
               50_MAR (representing 50% of outcomes are MAR)
               50_MCAR (representing 50% of outcomes are MCAR)
               50_MNAR (representing 50% of outcomes are MNAR).

The last folder named "Application to EDP" is used to apply the simulation results to the real Ebola dataset. There are one data file (.csv) and one code file (.R) in it.
