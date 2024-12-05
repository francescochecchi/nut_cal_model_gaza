#...............................................................................
### + GENERIC MODEL OF ACUTE MALNUTRITION AS A FUNCTION OF CALORIC INTAKE ++ ###
#...............................................................................

#...............................................................................
## ------ R SCRIPT TO READ DATASETS AND PARAMETERS AND MANAGE DATASETS  ----- ##
#...............................................................................


#...............................................................................  
### Reading in and/or setting general and crisis-specific parameters/data
#...............................................................................

  #...................................      
  ## Read generic parameters

    # Identify file name
    filename <- paste0(dir_path, "in/generic_pars.xlsx")
    
    # Read each worksheet
    x <- excel_sheets(filename)
    for (i in x[! x == "contents"]) {
      assign(i, data.frame(read_excel(filename, sheet = i)))}

  #...................................      
  ## Read crisis-specific data and parameters

    # Read each worksheet of the main file
    filename <- paste0(dir_path, "in/crisis_specs.xlsx")
    x <- excel_sheets(filename)
    for (i in x[! x %in% c("lists", "readme")]) {
      
      # read
      df <- data.frame(read_excel(filename, sheet = i))
    
      # fix dates
      x <- grep("date", colnames(df), value = T)
      if (length(x) != 0) { 
        for (j in x) {df[, j] <- suppressWarnings(as.Date(df[, j], "%d-%m-%Y"))}
      }
      
      # assign object
      assign(i, df)
    }
    
    # Individual parameters: read each as its own object, with the right format
    for (i in 1:nrow(pars)) {
      if (grepl("date", pars[i, "parameter"])) {
        assign(pars[i, "parameter"], as.Date(pars[i, "value"], "%d%b%Y") )
      } 
      else {
        if (! is.na(suppressWarnings(as.numeric(pars[i, "value"])))) {
          assign(pars[i, "parameter"], as.numeric(pars[i, "value"]))
        } else {
          assign(pars[i, "parameter"], pars[i, "value"])
        }
      }
    }

    # Also read Kcal intake output from separate analysis, if file is provided
    if (! is.na(kcal_capita_file)) {
      intake <- readRDS(paste0(dir_path, "in/", kcal_capita_file))
      intake$intake <- intake[, kcal_capita_var]
    }
    
    
  #...................................      
  ## Set or calculate some derived/other parameters

    # # Number of months to date
    # months_todate <- as.integer((date_start - date_crisis) / 30.44)

    # Geographical areas and their colours
    areas <- trimws(unlist(strsplit(areas, ",")))
    x <- palette_gen[c(4,12,2,10,14,6)]
    colours_areas <- x[1:length(areas)]

    # Risk/protective factor names
    factors <- data.frame(factor = c("intake", "tx_mam", "tx_sam", 
      "ari","dia","prop_formula", "req_ch_protected", "req_ad_protected"),
      factor_name=c("Kcal intake", "MAM treatment", "SAM treatment",
        "excess ARI period prevalence", "excess diarrhoea period prevalence", 
        "excess proportion of infants formula-fed",
        "proportion of child's intake safeguarded", 
        "proportion of adults' intake safeguarded"))

    # Proportion of pregnant and lactating women
    prop_preg <- prop_age_sex[which(prop_age_sex$age == "pregnant"), "female"]
    prop_lact <- prop_age_sex[which(prop_age_sex$age == "lactating"), "female"]
    
    # Extra intake requirement for pregnant and lactating women
    req_preg <- intake_req[which(intake_req$age == "pregnant"), "female"]
    req_lact <- intake_req[which(intake_req$age == "lactating"), "female"]
    
    # Convert age-sex proportion to long
##### SHOULD GENERALISE TO ALLOW FOR VARYING INPUT AGE CATEGORIES
    prop_age_sex <- reshape(
      subset(prop_age_sex, ! age %in% c("pregnant", "lactating") ), 
      direction = "long", varying = c("male", "female"), 
      times = c("male", "female"), v.names = "prop", idvar = "age", 
      timevar = "sex")
    prop_age_sex$id <- 1:nrow(prop_age_sex)
    
    # Maximum household size
    max_size <- max(hh_size$n_members)

    # Model outcomes to be tracked
    outcomes <- c("sam", "mam", "sam_inc", "mam_inc",
      "n_0659", "sam_0659", "mam_0659", 
      "n_0623", "sam_0623", "mam_0623", "whz_mean", "sam_tx", "mam_tx",
      "tx_sam_coverage", "tx_mam_coverage",
      "n_dia", "n_ari", "formula_fed", "n_u6mo", "intake")    
    
  #...................................      
  ## Set infectious disease parameters
    
    # Incidence and duration of disease gamma distribution parameters
    inc_scale_dia <- (2.97+6.75+6.88+5.46+2.34)/5
    dur_shape_dia <- (0.79+0.62+0.98)/3
    dur_scale_dia <- (2.69+3.07+5.98)/3

    inc_scale_ari <- (2.05+3.61+3.29+2.89)/4
    dur_shape_ari <- (1.26+1.04+0.79+0.86)/4
    dur_scale_ari <- (4.57+6.51+2.45+6.75)/4

    # Linear coefficient of correlation between episodes and duration
    b_corr <- 2/35
    # c_corr <- 2 # intercept (not needed)

    # Mean ratios of Kcal intake with/without diarrhoea / ARI
    kcal_ratio_dia <- mean(effect_illness[
      which(effect_illness$disease == "dia"), "kcal_ratio"])
    kcal_ratio_ari <- mean(effect_illness[
      which(effect_illness$disease == "ari"), "kcal_ratio"])

  #...................................      
  ## Generate dataset of population 6-59mo by area over time (if provided)
  if (exists("pop")) {  
    # Interpolate overall population to whole timeline
    df_po <- expand.grid(area = areas, date = as.Date(date_crisis : date_end))
    df_po$pop <- NA
    df_po <- df_po[order(df_po$area, df_po$date), ]
    for (i in areas) {
      x <- which(df_po$area == i)
      df_po[x, "pop"] <- approx(x = pop[which(pop$area == i), "date"], 
        y = pop[which(pop$area == i), "pop"], xout = df_po[x, "date"], rule=2)$y
    }

    # Calculate population 6-59mo
    x <- sum(prop_age_sex[which(prop_age_sex$age %in% c("0mo", "1 to 11mo",
      "12 to 59mo")), "prop"])
    df_po$pop_0059 <- df_po$pop * x
  }
        

#...............................................................................  
### Preparing metabolic parameters and data for children's weight loss model
    # (see Hall et al. 2013, http://dx.doi.org/10.1016/)
#...............................................................................

  #...................................      
  ## Prepare data from Fomon, Haschke et al. (1982): 10.1093/ajcn/35.5.1169
        # becomes base dataframe to which other variables are added below
    
    # Specify data
    df_ha <- fomon_weight_composition

    # Rescale some variables
    df_ha$sex <- factor(df_ha$sex, levels = c("female", "male"))
    df_ha$fm <- df_ha$fm / 1000
    df_ha$ffm <- df_ha$ffm / 1000

  #...................................      
  ## Compute delta parameter (relationship between age and 
        # energy expenditure from activity)
  
    # Necessary parameters
    delta_min <- 10 # Kcal/Kg/day - minimum activity level
    delta_max <- c(19, 17) # Kcal/Kg/day - max activity level for boys and girls

    # Compute activity level by age and sex; assume min during first 6mths, 
      # progressing to max from age 24mths; interpolate linearly otherwise
    df_ha[which(df_ha$age < 6), "delta"] <- delta_min
    df_ha[which(df_ha$sex == "male" & df_ha$age >= 24), "delta"] <- delta_max[1]
    df_ha[which(df_ha$sex =="female" & df_ha$age >= 24), "delta"]<- delta_max[2]
    x <- by(df_ha, df_ha$sex, function(x) {approx(x = x$age, y = x$delta, 
      xout <- sort(unique(x$age)), method = "linear")$y})
    df_ha$delta <- as.vector(unlist(x))

  #...................................      
  ## Compute g parameter (age-dependent adjustment to fat/lean mass 
        # partitioning ratio)
  
    # Compute age-specific changes in fat and fat-free mass, per day
    x <- by(df_ha, df_ha$sex, function(x) {c(NA, diff(x$fm))})
    df_ha$delta_fm <- as.vector(unlist(x))
    x <- by(df_ha, df_ha$sex, function(x) {c(NA, diff(x$ffm))})
    df_ha$delta_ffm <- as.vector(unlist(x))
    df_ha$age <- df_ha$age * 30.44
    x <- by(df_ha, df_ha$sex, function(x) {c(NA, diff(x$age))})
    df_ha$delta_age <- as.vector(unlist(x))
    df_ha$delta_fm <- df_ha$delta_fm / df_ha$delta_age
    df_ha$delta_ffm <- df_ha$delta_ffm / df_ha$delta_age
  
    # Other needed parameters
    df_ha$rho_fm <- 9400 # energy density of fat mass change
    df_ha$rhohat_ffm <- 4.3 * df_ha$ffm + 837 
      # effective energy density of ffm change
    df_ha$C <- 10.4 * df_ha$rhohat_ffm / df_ha$rho_fm 
      # Forbes body composition parameter
    df_ha$p <- df_ha$C / (df_ha$C + df_ha$fm) # partitioning ratio
    df_ha$K <- ifelse(df_ha$sex == "male", 800, 700) # expenditure constant
    gamma_fm <- 4.5 # metabolic rate of fat mass
    gamma_ffm <- 22.4 # metabolic rate of fat-free mass
    nu_fm <- 180 # cost of fat synthesis
    nu_ffm <- 230 # cost of fat-free tissue synthesis
    
    # Compute age-specific reference energy imbalance gap
    df_ha$eb_ref <- df_ha$rhohat_ffm * df_ha$delta_ffm + 
      df_ha$rho_fm * df_ha$delta_fm
    
    # Slope of the ratio in the changes of fm and ffm
    df_ha$alpha <- df_ha$delta_ffm / df_ha$delta_fm
    
    # g parameter
    df_ha$g <- (((df_ha$rhohat_ffm * df_ha$alpha * (1 - df_ha$p)) - 
      (df_ha$rho_fm * df_ha$p)) * df_ha$eb_ref) /
      (df_ha$rho_fm + df_ha$alpha * df_ha$rhohat_ffm)

  #...................................      
  ## Compute age-specific recommended kcal intake associated with normal growth

    # Recommended kcal intake associated with normal growth
    df_ha$intake_ref <- df_ha$eb_ref + df_ha$K + (gamma_ffm + df_ha$delta) * 
      df_ha$ffm + (gamma_fm + df_ha$delta) * df_ha$fm + 
      (nu_ffm / df_ha$rhohat_ffm) * (df_ha$p * df_ha$eb_ref + df_ha$g) +
      (nu_fm / df_ha$rho_fm) * ((1 - df_ha$p) * df_ha$eb_ref - df_ha$g)
    
    # Add extra energy requirement if formula-fed
    adj_formula$age <- adj_formula$age * 30.44
    df_ha <- merge(df_ha, adj_formula[, c("age", "sex", "adj_formula_fed")],
      by = c("age", "sex"), all.x = T)
    df_ha[which(df_ha$age > (18 * 30.44)), "adj_formula_fed"] <- 1
        
  #...................................      
  ## Interpolate key output parameters over all ages
  
    # Initialise output dataset
    age <- 0:(as.integer(10 * 365.25))
    out <- expand.grid(sex = c("female", "male"), age = age)
    out[, c("g", "intake_ref", "fm", "ffm", "delta", "adj_formula_fed")] <- NA
    out <- out[order(out$sex, out$age), ]
    
    # Interpolate, by sex
    df_ha_old <- df_ha
    for (i in c("g", "intake_ref", "fm", "ffm", "delta", "adj_formula_fed")) {
      df <- na.omit(df_ha_old)
      x <- by(df, df$sex, function(x) {predict(smooth.spline(x[, c("age", i)]),
        x = age)$y})
      out[, i] <- as.vector(unlist(x))
    }
  
    # Save new data
    df_ha <- out[, c("age", "sex", "delta", "g","adj_formula_fed","intake_ref")]
         

  #...................................      
  ## Visualise some metabolic quantities
    
    # Prepare data for plotting
    x <- df_ha_old[, c("age", "sex", "g", "intake_ref", "adj_formula_fed")]
    colnames(x) <- c("age","sex","g_obs","intake_ref_obs","adj_formula_fed_obs")
    x$age <- as.integer(x$age)
    df <- merge(out, x, by = c("age", "sex"), all.x = TRUE)
    
    # G parameter
    ggplot(df, aes(x = age / 30.44, y = g, colour = sex)) +
      geom_line(alpha = 0.5, linewidth = 1) +
      geom_point(aes(y = g_obs), alpha = 0.75, size = 3) +
      theme_bw() +
      scale_colour_manual(values = palette_gen[c(9,3)]) +
      scale_x_continuous("age (months)", breaks = seq(0, 60, 5), 
        limits = c(0, 60)) +
      scale_y_continuous("fat-free mass accretion preference (g)") +
      theme(legend.position = "top")
    ggsave(paste0(dir_path, "out/02_g_smooth.png"),
      dpi = "print", units = "cm", height = 15, width = 22)       
  
    # Recommended intake
    ggplot(df, aes(x = age / 30.44, y = intake_ref, colour = sex)) +
      geom_line(alpha = 0.5, linewidth = 1) +
      geom_point(aes(y = intake_ref_obs), alpha = 0.75, size = 3) +
      theme_bw() +
      scale_colour_manual(values = palette_gen[c(9,3)]) +
      scale_x_continuous("age (months)", breaks = seq(0, 60, 5),
        limits = c(0, 62), expand = c(0,0)) +
      scale_y_continuous("recommended intake (Kcal/day)", limits = c(NA, 1600)) +
      theme(legend.position = "top")
    ggsave(paste0(dir_path, "out/02_intake_ref_smooth.png"),
      dpi = "print", units = "cm", height = 15, width = 22)       
      
    # Extra energy requirement if formula-fed, by age
    ggplot(df, aes(x = age / 30.44, y = adj_formula_fed - 1,
      colour = sex, fill = sex)) +
      geom_point(aes(y = adj_formula_fed_obs - 1), alpha = 0.5, size = 3) +
      geom_line(linetype = "22", linewidth = 1) +
      theme_bw() +
      theme(legend.position = "top") +
      scale_x_continuous("age (months)", breaks = seq(0, 60, 6),
        limits = c(0, 18)) +
      scale_y_continuous(
        "added caloric requirement, relative to breastfed child", 
        labels = scales::percent_format()) +
      scale_colour_manual(values = palette_gen[c(9,3)]) +
      scale_fill_manual(values = palette_gen[c(9,3)])
    ggsave(paste0(dir_path, "out/02_formula_feeding_adj.png"),
      dpi = "print", units = "cm", height = 12, width = 20)


#...............................................................................  
### Computing recommended caloric intake by age-sex
#...............................................................................
    
  #...................................      
  ## Compute recommended Kcal intake for the population, by age-sex 

    # Convert age-sex recommended caloric intake to long
    intake_req <- reshape(
      subset(intake_req, ! age %in% c("pregnant", "lactating") ), 
      direction = "long", varying = c("male", "female"), 
      times = c("male", "female"), v.names = "req", idvar = "age", 
      timevar = "sex")

    # Merge with age-sex proportion
#### GENERALISE FOR VARYING AGE-SEX PROPORTION AGE CATEGORIES
    x <- prop_age_sex
    x$age_cat <- x$age
    x[which(x$id %in% c(1:3, 20:22)), "age_cat"] <- "0 to 4yo"
    x[which(x$id %in% c(7:14, 26:33)), "age_cat"] <- "20 to 59yo"
    x[which(x$id %in% c(15:19, 34:38)), "age_cat"] <- "60 to 100yo"
    x <- aggregate(list(prop = x$prop), by = x[, c("age_cat", "sex")], sum)
    colnames(x)[1] <- "age"
    hh_pars <- merge(intake_req[, c("age", "sex", "req")], x, 
      by = c("age", "sex"), all.x = TRUE)
    hh_pars$age <- factor(hh_pars$age, levels = c("0 to 4yo", "5 to 9yo", 
      "10 to 14yo", "15 to 19yo", "20 to 59yo", "60 to 100yo"))
    hh_pars <- hh_pars[order(hh_pars$age, hh_pars$sex), ]
    hh_pars$id <- 1:nrow(hh_pars)
    
    # Store parameter values among adults only (needed later)
    hh_pars_ad <- hh_pars
    hh_pars_ad$prop <- ifelse(!hh_pars_ad$age 
      %in% c("20 to 59yo", "60 to 100yo"), 0, hh_pars_ad$prop)
    hh_pars_ad$prop <- hh_pars_ad$prop / sum(hh_pars_ad$prop)

    # Amend recommended intake for children <10y to reflect Hall model estimates
    x <- df_ha
    x$age_cat <- cut(x$age, breaks = c(0, 30.44, 365.25, 365.25*5, 365.25*10),
      include.lowest = TRUE, right = FALSE, labels = c("0mo", "1 to 11mo",
        "12 to 59mo", "5 to 9yo"))
    x <- aggregate(list(intake_ref = x$intake_ref), 
      by = x[, c("sex", "age_cat")], FUN = mean)
    colnames(x) <- c("sex", "age", "intake_ref")
    x <- merge(x, prop_age_sex, by = c("sex", "age"), all.x = TRUE)
    x1 <- aggregate(list(x$prop), by = list(x$sex), FUN = sum)
    colnames(x1) <- c("sex", "prop_sum")
    x <- merge(x, x1, by = "sex", all.x = TRUE)
    x$wt <- x$prop / x$prop_sum
    x$age_cat <- ifelse(x$age %in% c("0mo", "1 to 11mo", "12 to 59mo"), 
      "0 to 4yo", "5 to 9yo")
    x1 <- aggregate(list(wt_sum = x$wt), by = x[, c("sex", "age_cat")], 
      FUN = sum)
    x <- merge(x, x1, by = c("sex", "age_cat"), all.x = TRUE)
    x$wt <- x$wt / x$wt_sum
    x$wt_prod <- x$intake_ref * x$wt    
    x <- aggregate(list(intake_ref = x$wt_prod), by = x[, c("sex", "age_cat")],
      FUN = sum)
    colnames(x) <- c("sex", "age", "intake_ref")
    hh_pars <- merge(hh_pars, x, by = c("sex", "age"), all.x = TRUE)
    hh_pars$req <- ifelse(is.na(hh_pars$intake_ref), hh_pars$req, 
      hh_pars$intake_ref)
    hh_pars_ad <- merge(hh_pars_ad, x, by = c("sex", "age"), all.x = TRUE)
    hh_pars_ad$req <- ifelse(is.na(hh_pars_ad$intake_ref), hh_pars_ad$req, 
      hh_pars_ad$intake_ref)
    
    # Compute mean recommended intake per capita
    req_intake <- weighted.mean(hh_pars$req, w = hh_pars$prop) + 
      prop_preg * req_preg + prop_lact * req_lact
    
    # Store recommended intake for infants < 6mths (needed later)
    x <- subset(df_ha, age < (30.44 * 6))
    req_u6mo <- aggregate(list(kcal = x$intake_ref), by = list(sex = x$sex), 
      FUN = mean)
    

#...............................................................................  
### ENDS
#...............................................................................
     