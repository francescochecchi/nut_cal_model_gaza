#...............................................................................
### + GENERIC MODEL OF ACUTE MALNUTRITION AS A FUNCTION OF CALORIC INTAKE ++ ###
#...............................................................................

#...............................................................................
## ------ R SCRIPT TO LOAD PACKAGES AND SOURCE OTHER ANALYSIS SCRIPTS  ------ ##
#...............................................................................



#...............................................................................
### Preparatory steps
#...............................................................................

  #...................................      
  ## Install or load required R packages
  if (!"pacman" %in% rownames(installed.packages())){install.packages("pacman")}
  
  pacman::p_load(
    anthro,        # Compute anthropometric scores
    doParallel,    # Support parallelised implementation
    foreach,       # Parallelise loops
    gamlss,        # Fit growth curves
    ggplot2,       # Data visualisation
    ggpubr,        # Arrange multiple plots into a single plot
    ggrepel,       # Improve labelling of plots
    gtools,        # Assist various programming tasks
    lubridate,     # Work with dates and times
    MASS,          # Various statistical methods
    mgcv,          # Fit GAM models
    pammtools,     # Produce uncertainty bands for step plots    
    parallel,      # Parallelise functions
    readxl,        # Read Excel files
    scales,        # Scale and format data
    tidyverse,     # Tidyverse suite of packages
    viridis,       # Colour palettes
    zoo)           # Compute running means

  #...................................      
  ## Starting setup

    # Clean up from previous code / runs
    rm(list=ls(all=TRUE) )
  
    # Set font for Windows or Mac
    suppressWarnings(windowsFonts(Arial = windowsFont("Arial")))
    suppressWarnings(par(family = "Arial"))
    
    # Set working directory to where this file is stored
    dir_path <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/")
    setwd(dir_path)
    print( getwd() )
    dir_path <- gsub("/code", "", dir_path)
    suppressWarnings(dir.create(paste0(dir_path, "out")))
    
    # Initialise random numbers
    set.seed(123)
    
    # Colour-blind palette for graphing
    palette_gen <- viridis(16)
    show_col(palette_gen)
          

#...............................................................................
### Sourcing dependent scripts
#...............................................................................
    
  #...................................      
  ## Run the following code always:
    
    # Specify functions (running time on a standard laptop: <1 min)
    source(paste0(dir_path, "code/01_specify_functions.r") )

    # Read and set parameters and data (<1 min)
    source(paste0(dir_path, "code/02_read_prepare_inputs.r") )

  #...................................      
  ## The following code needs to be run only once (CAREFUL: takes a long time!):
    
    # Estimate growth curves (up to 10-20 hours depending on 
        # size of growth monitoring dataset being modelled)
    source(paste0(dir_path, "code/11_estimate_growth_curves.r") )

    # Calibrate weight loss/gain model (20-40 hours)
    source(paste0(dir_path, "code/12_calibrate_wt_model.r") )
    
    # Validate weight loss/gain model (<5 min; OPTIONAL CODE)
    source(paste0(dir_path, "code/13_validate_wt_model.r") )
    
    # Model relationship of infectious disease prevalence and incidence (20 min)
    source(paste0(dir_path, "code/14_model_infectious_disease.r") )
    
    # Model caloric sacrifice by adults in favour of children (20 min)
    source(paste0(dir_path, "code/15_model_caloric_sacrifice.r") )
          
  #...................................      
  ## Run the following code always:
    
    # Prepare timeline of protective/risk factors for simulation (<1 min) 
    source(paste0(dir_path, "code/21_prepare_timeline.r") )
   
    # Run retrospective and, if desired, scenario-based simulation (*) 
    source(paste0(dir_path, "code/22_run_scenarios.r") )
        # (*) computation time depends on number of children and runs,
        # but is expected to be high and may require a high-performance cluster
    
    
#...............................................................................
### Code structure
#...............................................................................

      # code sequence:
      # 00 ----> 01 -> 02 ----> 11 -> 12 (-> 13) -> 14 -> 15 ----> 21 -> 22
      # or, if codes 11-15 have been run once already:
      # 00 ----> 01 -> 02 ---------------------------------------> 21 -> 22
    
      # code dependence (which rows depend on which columns):
      #             | 00 | 01 | 02 | 11 | 12 | 13 | 14 | 15 | 21 | 22 |
      # 01          | X  |    |    |    |    |    |    |    |    |    |
      # 02          | X  |    |    |    |    |    |    |    |    |    |
      # 11          | X  | X  |    |    |    |    |    |    |    |    |
      # 12          | X  | X  | X  | X  |    |    |    |    |    |    |
      # 13          | X  | X  | X  | X  | X  |    |    |    |    |    |
      # 14          | X  | X  | X  |    |    |    |    |    |    |    |
      # 15          | X  | X  | X  |    |    |    |    |    |    |    |
      # 21          | X  |    | X  |    |    |    |    |    |    |    |
      # 22          | X  | X  | X  | X  | X  |    | X  | X  | X  |    |
    
      # function dependence (which functions are needed in which code):
      #             | 00 | 01 | 02 | 11 | 12 | 13 | 14 | 15 | 21 | 22 |
      # f_bw_fm     |    |    |    |    | X  | X  |    |    |    | X  |
      # f_co        |    |    |    |    |    |    |    |    |    | X  |
      # f_hall_kids |    |    |    |    | X  | X  |    |    |    | X  |
      # f_hall_pred |    |    |    |    | X  | X  |    |    |    | X  |
      # f_prev      |    |    |    |    |    |    | X  |    |    |    |
      # f_sacr      |    |    |    |    |    |    |    | X  |    |    |
      # f_sim       |    |    |    |    |    |    |    |    |    | X  |
      
      # function cascade (which functions call other functions):
      #
      # f_hall_pred    calls  f_hall_kids, f_bw_fm
      # f_co           calls  f_hall_pred           calls   f_hall_kids, f_bw_fm
      # f_sim          calls  f_hall_kids, f_bw_fm

#...............................................................................  
### ENDS
#...............................................................................
     