#...............................................................................
### + GENERIC MODEL OF ACUTE MALNUTRITION AS A FUNCTION OF CALORIC INTAKE ++ ###
#...............................................................................

#...............................................................................
## ----- R SCRIPT TO SPECIFY BESPOKE FUNCTIONS NEEDED FOR THE ANALYSIS  ----- ##
#...............................................................................


#...............................................................................
### Function to estimate fat mass at birth as a function of sex and birth weight
# based on regression models in Eriksson et al. (2010) - first visit data
# https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1651-2227.2009.01665.x
#...............................................................................

f_bw_fm <- function(sex_f, birth_weight_f) {
  
  # Estimate percent fat mass
    # girls
    if (sex_f %in% c("female", "f", "F", 2, "2")) {
      prop_fm <- (-4.60 + 5.12 * birth_weight_f) / 100
      fm <- birth_weight_f * prop_fm
    }
    
    # boys
    if (sex_f %in% c("male", "m", "M", 1, "1")) {
      prop_fm <- (-5.04 + 4.64 * birth_weight_f) / 100
      fm <- birth_weight_f * prop_fm
    }
  
  # Return fat mass in Kg
  return(fm)
}


    
#...............................................................................  
### Function to initialise a new cohort of children up to the crisis' start date
#...............................................................................

f_co <- function(obj_co_f = obj_co, parallel = F) {

  #...................................      
  ## Unpack objects needed within the function
  for (i in names(obj_co_f)) {assign(i, obj_co_f[[i]])}
  
  #...................................      
  ## Initialise cohort with random ages, sex and growth centiles, by area
  
    # Number of children (n_kids per area)
    x <- n_kids * length(areas)
  
    # Initialise cohort
    co <- data.frame(
      area = sort(rep(areas, n_kids)),
      child = 1:x,  
      age = sample(0:as.integer(5*365.25 - 1), x, replace = T),
      sex = sample(c("female", "male"), x, replace = T),
      centile_weight = rep(paste0("c", 0:99), times = x/100),
      weight = NA, 
      ffm = NA, 
      fm = NA, 
      height = NA,
      formula_fed = F,
      dia = F,
      ari = F
    )
    co$centile_weight <- as.character(co$centile_weight)
    co$age_m <- floor(co$age / 30.44)
    
####IMPROVE AGE DISTRIBUTION? (MINOR ISSUE)

    # Add weight growth centile calibration parameters
    co <- merge(co, cal_ml[, c("sex", "centile_weight", 
      grep("par", colnames(cal_ml), value = T))], 
      by = c("sex", "centile_weight"), all.x = T)
    co <- co[order(co$child), ]

    # Select random height trajectories from fitted correlation model
    for (i in 1:nrow(co)) {
      x <- corr_mat[which(corr_mat$sex == co[i, "sex"] & 
        corr_mat$centile_weight == co[i, "centile_weight"]), ]
      co[i, "centile_height"]  <- sample(x$centile_height, size = 1, 
        prob = x$prob_fitted)
    }
        
  #...................................      
  ## Compute weight and height for each child up to current time (single loop)
  if (! parallel) {
    
    # Progress bar
    pb <- txtProgressBar(min = 1, max = nrow(co), style = 3)

    for (i in 1:nrow(co)) {

    # Update progress
    setTxtProgressBar(pb, i)
    
    # Predict weight at child's age
    x <- f_hall_pred(data_f = co[i, ], par_f = co[i, grep("par", colnames(co))], 
      intake_f = "actual")
    co[i, c("ffm", "fm", "weight")] <- x[which(x$age == co[i, "age"]), 
      c("ffm", "fm", "weight")]
    
    # Report height at child's age
    co[i, "height"] <- centiles[which(centiles$sex == co[i, "sex"] &
      centiles$measure== "height" & centiles$age == co[i, "age"]), 
      co[i, "centile_height"]]
    }
    close(pb)
  } 
    
  #...................................      
  ## Compute weight and height for each child up to current time (parallelised)
  if (parallel) {
    # Prepare for parallel processing
    sets <- split(co, seq(nrow(co))) 
    
    # Parallel function
    f_parallel <- function(set_i, centiles_f = centiles,
      f_hall_pred_f = f_hall_pred){
      
      # identify calibration parameters
      par_i <- set_i[grep("par", names(set_i))]
      
      # prepare output dataset
      df_i <- centiles_f[which(centiles_f$sex == set_i[["sex"]] &
        centiles_f$measure== "weight"), c("age", set_i[["centile_weight"]])]
      colnames(df_i) <- c("age", "weight")
      df_i[which(df_i$age > 0), "weight"] <- NA
      df_i$sex <- set_i[["sex"]]
      df_i$formula_fed <- F
      df_i$dia <- F
      df_i$ari <- F
      
      # predict fat-free mass, fat mass and weight for all ages
      x <- f_hall_pred_f(data_f = df_i, f_hall_kids_f = f_hall_kids, 
        df_ha_f = df_ha, f_bw_fm_f = f_bw_fm, par_f = par_i, 
        intake_f = "actual")
      df_i[, c("ffm", "fm", "weight")] <- x[, c("ffm", "fm", "weight")]

      # report FFM, FM and weight at child's age
      set_i[c("ffm", "fm", "weight")] <- 
        df_i[which(df_i$age == set_i[["age"]]), c("ffm", "fm", "weight")]

      # report height at child's age
      set_i[["height"]] <- centiles_f[which(centiles_f$sex == set_i[["sex"]] &
        centiles_f$measure== "height" & centiles_f$age == set_i[["age"]]), 
        set_i[["centile_height"]]]
      
      # return output
      return(set_i)
    }
    
    # Prepare parallel cores
    cl <- parallel::makeCluster(detectCores())
    x <- c("df_ha", "centiles", "f_hall_kids", "f_bw_fm","f_hall_pred","cal_ml")
    clusterExport(cl, x, envir=environment())  
    clusterEvalQ(cl, library("gamlss"))    
    
    # Run and collect output
    co <- parallel::parLapply(cl, sets, f_parallel) # 100 children = ~1 min 
    stopCluster(cl)
    co <- do.call(rbind, co)
  }

  #...................................      
  ## Add crisis-relevant variables to cohort and return
    
    # Daily caloric intake
    co$intake <- NA
    
    # Starting anthropometric status
      # anthropometric indices
      co$sex_anthro <- ifelse(co$sex == "female", "f", "m")
      x <- with(co, anthro_zscores(sex = sex_anthro, age = age, 
        lenhei= height, weight = weight))
      co <- cbind(co, x)

      # SAM/MAM status now
      co[, c("mam", "sam")] <- F
      co[which(co$zwfl < (-2) & co$zwfl >= (-3)), "mam"] <- T
      co[which(co$zwfl < (-3)), "sam"] <- T
    
      # SAM/MAM time since becoming a case
      co[, c("mam_time", "sam_time")] <- NA

      # SAM/MAM random treatment delays
      co$mam_tx_delay <- rlnorm(nrow(co), delay_tx_mam, delay_tx_mam/2)
      co$sam_tx_delay <- rlnorm(nrow(co), delay_tx_sam, delay_tx_sam/2)
      
    # SAM/MAM treatment 
      # treatment status now
      co[, c("sam_tx", "mam_tx")] <- F
      
      # time since reverting from MAM to OK, or from SAM to MAM
      co[, c("recovery_time", "sam_to_mam_time")] <- NA
      
      # therapeutic food dosage
      co$therapeutic_food <- 0

    # Diarrhoea/ARI occurrence
      # illness status now
      co[, c("dia", "ari")] <- F

      # expected episodes per year, relative probability weight of falling ill,
        # time since symptom onset and episode duration
      co[, c("dia_time", "ari_time")] <-NA
      co[, c("dia_wt", "ari_wt", "dia_y", "ari_y", "dia_dur", "ari_dur")] <- 0
      
  # Formula feeding - determined by uniform random number x: any child with x
        # below the excess population prevalence is formula-fed (area-specific)
    co$formula_fed <- F
    co$formula_fed_p <- runif(nrow(co))
      
    # Return starting cohort
    return(co)         
}
    

#...............................................................................  
### Function encoding child weight loss model: adapted from
  # Hall KD, Butte N, Swinburn BA, et al.
  # Dynamics of childhood growth and obesity: development and validation of 
  # a quantitative mathematical model.
  # Lancet 2013 https://doi.org/10.1016/S2213-8587(13)70051-2 (see webappendix)
#...............................................................................

f_hall_kids <- function(df_f, par_f, kcal_ratio_dia_f = kcal_ratio_dia, 
  kcal_ratio_ari_f = kcal_ratio_ari,
    # other fixed parameters:
    forbes = 10.4, # Forbes' constant
    rho_fm = 9400, # energy density of fat mass change (Kcal/Kg)
    gamma_fm = 4.5, # metabolic rate of fat mass (Kcal/kg/day)
    gamma_ffm = 22.4, # metabolic rate of fat-free mass (Kcal/kg/day)
    nu_fm = 180, # cost of fat synthesis (Kcal/kg)
    nu_ffm = 230, # cost of fat-free tissue synthesis  (Kcal/kg)
    beta = 0.24 # adaptive thermogenesis
  ) {  

  #...................................      
  ## Read or compute parameters from data
    
    # Read parameters from data
    weight <- df_f[["weight"]] # weight now
    fm <- df_f[["fm"]] # fat mass now
    ffm <- df_f[["ffm"]] # fat-free mass now
    intake <- df_f[["intake"]] # intake (Kcal/day)
    sex <- df_f[["sex"]] # sex
    age <- df_f[["age"]] # age in days
    formula_fed <- df_f[["formula_fed"]]
    delta <- df_f[["delta"]] # mobility
    g <- df_f[["g"]] # growth parameter
    intake_ref <- df_f[["intake_ref"]] # recommended intake
    adj_formula_fed <- df_f[["adj_formula_fed"]] # adjustment for formula fed
    dia <- df_f[["dia"]]
    ari <- df_f[["ari"]]

    # Read individual variability parameter values
    par0 <- as.numeric(par_f[1])
    par1 <- as.numeric(par_f[2])

    # Adjust reference intake for formula feeding
    if (age < as.integer(30.44 * 6) & formula_fed) {
      intake_ref <- intake_ref * adj_formula_fed
    }
    
    # Set intake to be equal to the reference if intake is missing 
      # (for testing model)
    if (is.na(intake)) {intake <- intake_ref}

    # Reduce intake if the child has diarrhoea and/or ARI
    if (ari) {intake <- intake * kcal_ratio_ari_f}
    if (dia) {intake <- intake * kcal_ratio_dia_f}

    # Derived metabolic parameters
    K <- ifelse(sex == "female", 700, 800) # expenditure constant (Kcal/day)
    rhohat_ffm <- 4.3 * ffm + 837 # effective energy density of ffm change
    C <- forbes * rhohat_ffm / rho_fm # Forbes body composition parameter
    p <- C / (C + fm) # partitioning ratio
    intake_diff <- intake - intake_ref # difference between actual and reference
      # intake
    
  #...................................      
  ## Compute energy expenditure
    
    # Uncalibrated energy expenditure
    energy <- (K + (gamma_ffm + delta) * ffm + (gamma_fm + delta) * fm + 
      beta * intake_diff + 
      ((nu_ffm / rhohat_ffm) * p + (nu_fm / rho_fm) * (1 - p)) * intake +
      g * ((nu_ffm / rhohat_ffm) - (nu_fm / rho_fm)) ) /
      (1 + (nu_ffm / rhohat_ffm) * p + (nu_fm / rho_fm) * (1 - p))

    # Calibrate energy expenditure for individual variability
    mm <- par0 + par1
    energy <- energy + energy * as.numeric(mm)
    
  #...................................      
  ## Compute and output change in fat mass, fat-free mass and thus weight
    
    # Change
    delta_ffm <- (p * (intake - energy) + g) / rhohat_ffm
    delta_fm <- ((1 - p) * (intake - energy) - g) / rho_fm
    delta_weight <- delta_ffm + delta_fm

    # Output new masses / weight
    return(c(ffm + delta_ffm, fm + delta_fm, weight + delta_weight))
}  


#...............................................................................  
### Function to predict weight trajectory for a single child across age / time
#...............................................................................

f_hall_pred <- function(data_f, f_hall_kids_f = f_hall_kids, df_ha_f = df_ha,
  f_bw_fm_f = f_bw_fm, centiles_f = centiles,
  par_f, # individual calibration parameters for boys and girls 
  intake_f = "reference" # options: "reference", "actual"
  ) {  
 
  #...................................      
  ## Prepare the timeline
  
    # Set up timeline for child
    timeline <- data.frame(
      age = 0:max(data_f$age),
      sex = unique(data_f$sex), 
      weight = NA,
      fm = NA,
      ffm = NA,
      formula_fed = unique(data_f$formula_fed),
      dia = unique(data_f$dia),
      ari = unique(data_f$ari)
    )
    
    # Merge in reference and actual intake parameters
    timeline <- merge(timeline,df_ha_f[which(df_ha_f$sex==unique(data_f$sex)),], 
      by = c("age", "sex"), all.x = T)
    if (intake_f == "reference") {timeline$intake <- timeline$intake_ref}
    if (intake_f == "actual") {timeline$intake <- timeline$intake_act}
    
    # Starting weight, FM and FFM at birth
    x <- which(timeline$age == 0)
    timeline[x, "weight"] <- centiles_f[
      which(centiles_f$sex == unique(data_f$sex) & 
        centiles_f$measure == "weight" & 
        centiles_f$age == 0), 
      unique(data_f$centile_weight)]
    timeline[x, "fm"] <- f_bw_fm_f(unique(timeline$sex), timeline[x, "weight"])
    timeline[x, "ffm"] <- timeline[x, "weight"] - timeline[x, "fm"]  
    
    # Sort
    timeline <- timeline[order(timeline$age), ]
 
    # Prepare parameter values
    x <- dLOGNO(timeline$age / 100, as.numeric(par_f["par1"]))
    par_timeline <- suppressWarnings(data.frame(par_f["par0"], 
      as.numeric(par_f["par0"]) * x/max(x) ))
      # change sign for centiles >= 50% (= negative energy penalty)
      if (as.integer(gsub("c", "", unique(data_f$centile_weight))) >= 50)
        {par_timeline <- -par_timeline}
      
  #...................................      
  ## Predict weight for each day in the timeline
  for (k in 1:max(nrow(timeline) -1, 1)) {
    
    # Identify data now
    x <- timeline[k, ]
    
    # Compute new weight, FM and FFM
    timeline[k + 1, c("ffm", "fm", "weight")] <- f_hall_kids_f(df_f = x,
      par_f = par_timeline[k, ])
  }
  
    # Output
    return(timeline)
}  



#...............................................................................  
### Function to predict longitudinal and period prevalence of disease 
    # as a function of incidence and disease duration; disease parameters from 
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2829935/#sec006
#...............................................................................

f_prev <- function(inc_shape_f = log(10), n_kids_f = 1000, period_f = 14,
  disease_f = "dia") {

  #...................................
  ## Set values of disease parameters
  
  if (disease_f == "dia") {
    # Diarrhoea: mean of Table 1, excluding all ages and dysentery (Thailand)
    inc_shape <- inc_shape_f
    inc_scale <- (2.97+6.75+6.88+5.46+2.34)/5
    dur_shape <- (0.79+0.62+0.98)/3
    dur_scale <- (2.69+3.07+5.98)/3
  }
  
  if (disease_f == "ari") {
    # ARI: mean of Table 1, other conditions
    inc_shape <- inc_shape_f    
    inc_scale <- (2.05+3.61+3.29+2.89)/4
    dur_shape <- (1.26+1.04+0.79+0.86)/4
    dur_scale <- (4.57+6.51+2.45+6.75)/4
  }
  
    # Intercept and linear coefficient of correlation 
      # between episodes and duration (eyeballed from Figure 3)
    b_corr <- 2/35
    c_corr <- 2
  
  
  #...................................      
  ## Set up a dataset of children's disease status over time
    
    # Matrix of children's status over time
    kids <- matrix(0, nrow = n_kids_f, ncol = 365 + 100)
  
    # Children's individual parameter values
    ind <- data.frame(id = 1:nrow(kids),
      n_ill = NA, dur_ill = I(vector("list", nrow(kids))) )
    
  #...................................      
  ## Attribute disease episodes and their duration to each child

    # Number of episodes per year
    ind$n_ill <- as.integer(
      rgamma(nrow(kids), shape = inc_shape, scale = inc_scale))
    
    # Episode duration, adjusted for correlation between episodes and duration
      # and intra-individual variability (see paper)
    for (i in 1:nrow(ind)) {
      dur_ill <- rgamma(ind[i, "n_ill"], shape = dur_shape, scale = dur_scale)
      dur_ill <- dur_ill + rnorm(1, mean = 1, sd = 0.5)
      dur_ill <- dur_ill + (ind[i, "n_ill"] - mean(ind$n_ill)) * b_corr
####NEED A SECOND LOOK (normalise? but likely a negligible source of error)    
      dur_ill[which(dur_ill <= 0)] <- 1
      ind[i, "dur_ill"][[1]] <- list(as.integer(dur_ill))
    }    

    # Position the episodes across the year
    for (i in 1:nrow(kids)) {
      if (ind[i, "n_ill"] > 0) {
        out <- c()
        x1 <- sample(1:365, ind[i, "n_ill"], replace = F)
        x2 <- unlist(ind[i, "dur_ill"])
        for (j in 1:length(x1)) {
          out <- c(out, x1[j]:(x1[j] + x2[j] - 1))
        }
        kids[i, unique(out)] <- 1
      }
    }
      
  #...................................      
  ## Compute and output mean period and longitudinal prevalences

    # Period prevalence (only compute in 2nd half of year to ensure equilibrium)
    x <- t(apply(kids[, 183:365], 1, rollmax, period_f, align = "right"))
    period_prev <- mean(colSums(x)/nrow(kids))
    
    # Longitudinal/mean point prevalence
    point_prev <- mean(kids[, 1:365])

    # Output
    out <- c(period_prev, point_prev)
    names(out) <- c("period_prev", "point_prev")
    return(out)
}


#...............................................................................  
### Function to simulate households and adults' Kcal sacrifice to feed children
#...............................................................................

f_sacr <- function(hh_f = hh, max_size_f = max_size, hh_size_f = hh_size,
  req_ch_protected_f = req_ch_protected, req_ad_protected_f = req_ad_protected, 
  kcal_capita_f = kcal_capita_sim, hh_pars_f = hh_pars, 
  hh_pars_ad_f = hh_pars_ad, prop_lact_f = prop_lact, prop_preg_f = prop_preg,
  req_lact_f = req_lact, req_preg_f = req_preg, req_u6mo_f = req_u6mo) {
  
  #...................................      
  ## Attribute random household (HH) member characteristics
  
    # Size of each HH
    hh_f$size <- sample(x = max_size_f, size = nrow(hh_f), replace = T,
      prob = hh_size_f$prop)
    
    # Age-sex of HH "head" (can only be aged 20+)
    hh_f$m1 <- sample(hh_pars_ad_f$id, size = nrow(hh_f), replace = T,
      prob = hh_pars_ad_f$prop)
    
    # Now adjust age-sex distribution by taking away HH heads
    x <- data.frame(id = as.integer(names(table(hh_f$m1))), 
      prop_head = as.vector(table(hh_f$m1) / sum(hh_f$size)) )
    hh_pars_f <- merge(hh_pars_f, x, by = "id", all.x = T)
    hh_pars_f[which(is.na(hh_pars_f$prop_head)), "prop_head"] <- 0
    hh_pars_f$prop <- hh_pars_f$prop - hh_pars_f$prop_head
      # sometimes adults are oversampled, leading to a (small) negative value
      hh_pars_f[which(hh_pars_f$prop < 0), "prop"] <- 0
    hh_pars_f$prop <- hh_pars_f$prop / sum(hh_pars_f$prop)
    
    # Then attribute age-sex of further HH members
    x <- sapply(hh_f$size, function(xx) {sample(hh_pars_f$id, size = xx - 1, 
      replace = T, prob = hh_pars_f$prop)})
    hh_f[, paste0("m", 2:max_size_f)] <- 
      t(list2DF(lapply(x, `length<-`, max_size_f - 1)))
    hh_pars_f <- hh_pars[order(hh_pars$id), ]
    
  #...................................      
  ## Work out age, sex, whether infant <6mo and Kcal requirement of HH members
  for (i in 1:max_size_f) {
    # Age
    hh_f[, paste0("age", i)] <- sapply(hh_f[, paste0("m", i)],
      function(x) {hh_pars_f[x, "age"]})
    
    # Infant aged <6mo?
    hh_f[, paste0("u6mo", i)] <- 0
    x <- which(hh_f[, paste0("age", i)] == "0 to 4yo")
    if (length(x) > 0) {
      hh_f[x, paste0("u6mo", i)] <- rbinom(length(x), 1, 6/59) }

    # Sex
    hh_f[, paste0("sex", i)] <- sapply(hh_f[, paste0("m", i)],
      function(x) {hh_pars_f[x, "sex"]})
    
    # Kcal requirement (zero if infant <6mo, meaning no external food needed)
    hh_f[, paste0("req", i)] <- sapply(hh_f[, paste0("m", i)],
      function(x) {hh_pars_f[x, "req"]})
    x <- which(hh_f[, paste0("u6mo", i)] == 1)
    hh_f[x, paste0("req", i)] <- 0
  }
  
    # Identify adults
    ad <- apply(hh_f, 1, function(x) {
      which(x[paste0("age", 1:max_size_f)] %in% c("20 to 59yo","60 to 100yo"))})
    
    # Identify small children
    ch <- apply(hh_f, 1, function(x) {
      which(x[paste0("age", 1:max_size_f)] %in% c("0 to 4yo", "5 to 9yo"))})
    
    # How many pregnant women there are in the HH
    hh_f$n_preg <- rpois(hh_f$size, prop_preg_f)
    
    # How many lactating women there are in the HH
    hh_f$n_lact <- rpois(hh_f$size, prop_lact_f)
  
  #...................................      
  ## Work out required and actual Kcal available
  
    # Kcal requirement of HH
    hh_f$kcal_needed <- rowSums(hh_f[, paste0("req", 1:max_size_f)], 
      na.rm = T) + hh_f$n_preg * req_preg_f + hh_f$n_lact * req_lact_f
    hh_f$kcal_share_needed <- hh_f$kcal_needed / sum(hh_f$kcal_needed)
    
    # Available Kcal
    hh_f$kcal_avail <- kcal_capita_f * sum(hh_f$size) * hh_f$kcal_share_needed
  
  #...................................      
  ## Attribute Kcal among HH members
  
    # First distribute Kcal among all HH members based on their need
    x <- which(hh_f$kcal_needed != 0)
    for (i in 1:max_size_f) {
      hh_f[x, paste0("act", i)] <- hh_f[x, "kcal_avail"] * 
        hh_f[x, paste0("req", i)] / hh_f[x, "kcal_needed"]
      hh_f[-x, paste0("act", i)] <- hh_f[-x, paste0("req", i)]      
    }
    
    # Next, work out spare adult Kcal that can be sacrificed to feed children
    for (i in 1:nrow(hh_f)) {
      hh_f[i, "act_ad"] <- sum(hh_f[i, paste0("act", ad[[i]])])
      hh_f[i, "act_spare"] <- pmax(0, sum(hh_f[i, paste0("act", ad[[i]])]) - 
          sum(hh_f[i, paste0("req", ad[[i]])]) * req_ad_protected_f)
    }
    
    # Then work out gap, if any, to reach 'acceptable' intake in children
      # total HH deficit among children 6mo-9yo
      for (i in 1:nrow(hh_f)) {
        hh_f[i, "req_ch"] <- sum(hh_f[i, paste0("req", ch[[i]])])
        hh_f[i, "act_ch"] <- sum(hh_f[i, paste0("act", ch[[i]])])
      }
      hh_f$deficit <- hh_f$req_ch * req_ch_protected_f - hh_f$act_ch
      hh_f$deficit <- ifelse(hh_f$deficit < 0, 0, hh_f$deficit)
      
      # accordingly, how much of the spare must be used
      hh_f$spare_used <- pmin(hh_f$act_spare, hh_f$deficit)
    
    # Now use the spare
    for (i in 1:nrow(hh_f)) {
      
      #...use the spare to fill each child's deficit...
      if (length(ch[[i]]) > 0 & sum(hh_f[i, paste0("req", ch[[i]])]) != 0)  {
        x <- hh_f[i, paste0("req", ch[[i]])] / 
          sum(hh_f[i, paste0("req", ch[[i]])])
        hh_f[i, paste0("act", ch[[i]])] <- 
          hh_f[i, paste0("act", ch[[i]])] + hh_f[i, "spare_used"] * x
        
        #...and subtract it from the adults
        x <- hh_f[i, paste0("req", ad[[i]])] / 
          sum(hh_f[i, paste0("req", ad[[i]])])
        hh_f[i, paste0("act", ad[[i]])] <- 
          hh_f[i, paste0("act", ad[[i]])] - hh_f[i, "spare_used"] * x
      }
    }
  
  #...................................      
  ## Final summing up and outputs
  
    # Compute adjusted Kcal intake for kids and adults
    for (i in 1:nrow(hh_f)) {
      hh_f[i, "act_ch_adj"] <- sum(hh_f[i, paste0("act", ch[[i]])])
      hh_f[i, "act_ad_adj"] <- sum(hh_f[i, paste0("act", ad[[i]])])
    }
    
    # Return needed output  
    out <- c(        
      # ratio of Kcal intake among kids with and without adult Kcal sacrifice
      sum(hh_f$act_ch_adj) / sum(hh_f$act_ch),
      # ratio of Kcal intake among adults with and without adult Kcal sacrifice
      sum(hh_f$act_ad_adj) / sum(hh_f$act_ad)
    )
    return(out)
}



#...............................................................................  
### Function to take the simulated child cohort through the crisis period
      # one run only; does not include scenario-based projection
#...............................................................................

f_sim <- function(co_f = co, tl_f = rtl, obj_sim_f = obj_sim) {

  #...................................      
  ## Prepare inputs
  
    # Unpack objects needed within the function
    for (ii in names(obj_sim_f)) {assign(ii, obj_sim_f[[ii]])}
  
    # Select random per capita caloric intake, 
        # if not already specified
    if (all(is.na(tl_f$intake))) {
      tl_f <- subset(tl_f, select = -intake)
      intake <- data.frame(intake)
      x <- sample(1:max(intake$run), 1, replace = T)
      x <- intake[which(intake$run == x), ]    
      tl_f <- merge(tl_f, x[, c("date", "area", "intake")], 
        by = c("date", "area"), all.x = T)
      tl_f <- tl_f[order(tl_f$area, tl_f$date), ]
    }
  
  # >> ...................................      
  ## >> For each area (start of long loop ii)...
  for (ii in areas) {
      
      # Progress
      print(paste0("    area: ", ii))    
      
      # Subset timeline for this day and area and cohort for this area
      tl_ii <- tl_f[which(tl_f$area == ii), ]
      co_ii <- co_f[which(co_f$area == ii), ]
      
      # Progress bar
      pb <- txtProgressBar(min = 1, max = nrow(tl_ii), style = 3)

  # >>>>...................................      
  ## >>>> For each time unit (start of long sub-loop jj)...
    for (jj in 1:nrow(tl_ii)) {
        
        # Update progress
        setTxtProgressBar(pb, jj)
        
      #...................................      
      ## Work out caloric intake for each child
        
        # Update metabolic parameters to reflect child's new age
        co_ii <- co_ii[, ! colnames(co_ii) %in% 
          colnames(df_ha)[! colnames(df_ha) %in% c("age", "sex")]]
        co_ii <- merge(co_ii, df_ha, by = c("age", "sex"), all.x = T)
        co_ii <- co_ii[order(co_ii$child), ]
    
        # Implement caloric adjustment due to adult sacrifice
          # find right row in caloric sacrifice database
          x <- unlist(tl_ii[jj, c("intake", "req_ch_protected",
            "req_ad_protected")])
          x <- sacrifice[, c("kcal_capita", "req_ch_protected",
            "req_ad_protected")] - t(replicate(nrow(sacrifice), x))
          x <- abs(x)
          x <- as.data.frame(apply(x, 2, function(x) {x / sum(x)}))
          x$error <- as.vector(rowSums(x))
          
          # work out caloric adjustment multiplier
          if (tl_ii[jj, "intake"] < req_intake) {
            kcal_adj <- sacrifice[which.min(x$error), "kcal_ch_adj"]
          } else 
          {kcal_adj <- 1}
    
        # Compute caloric intake for each child
        co_ii$intake <- pmin(co_ii$intake_act, co_ii$intake_ref * kcal_adj * 
          tl_ii[jj, "intake"] / req_intake) + co_ii$therapeutic_food
      
        # For infants <6mo (breastfed or formula-fed), maintain baseline intake
        x <- which(co_ii$age < (30.44 * 6))
        co_ii[x,"intake"] <- co_ii[x,"intake_act"] + co_ii[x,"therapeutic_food"]
    
        
      #...................................      
      ## Update non-exclusive breastfeeding status
        
        # Reset and update
        co_ii$formula_fed <- F
        x <- which(co_ii$age < (30.44 * 6) & co_ii$formula_fed_p < 
          tl_ii[jj, "prop_formula"])
        co_ii[x, "formula_fed"] <- T

      #...................................      
      ## Model diarrhoea and ARI
      for (kk in c("ari", "dia")) {
        
        # Which children newly fall ill? For how long?
        if (tl_ii[jj, paste0(kk, "_point_prev")] > 0) {   
          
          # relative probability of disease per child, under current incidence
            # number of expected annual episodes
            co_ii[, paste0(kk, "_y")] <- rgamma(nrow(co_ii), 
              shape = tl_ii[jj, paste0(kk, "_incidence")], 
              scale = get(paste0("inc_scale_", kk)))
         
            # relative probability of falling ill among healthy children now
            x <- which(! co_ii[, kk])
            co_ii[x, paste0(kk, "_wt")] <- co_ii[x, paste0(kk, "_y")] / 
              mean(co_ii[x, paste0(kk, "_y")])
            co_ii[which(co_ii[, kk]), paste0(kk, "_wt")] <- 0  
              # if ill now, weight = 0
            
          # work out how many free excess disease 'slots' there are
          x <- 0
          if (sum(! co_ii[, kk]) > 0) {
            x <- max(round(nrow(co_ii) * tl_ii[jj,paste0(kk,"_point_prev")],0) -
              sum(co_ii[, kk]), 0)
          }
          if (x > 0) {
            # attribute the free slots and start symptom clock
            x <- sample(1:nrow(co_ii), x, prob = co_ii[, paste0(kk, "_wt")])
            co_ii[x, kk] <- T
            co_ii[x, paste0(kk, "_time")] <- 0
     
            # set disease episode durations:
              # set a random disease episode duration...
              co_ii[x, paste(kk, "_dur")] <- 5
               
              # as.numeric(rgamma(length(x), 
              #   shape = get(paste0("dur_shape_", kk)), 
              #   scale = get(paste0("dur_scale_", kk))))
   
              # ...adjusted for intra-individual variability...
              co_ii[x, paste0(kk, "_dur")] <- co_ii[x, paste0(kk, "_dur")] + 
                rnorm(length(x), mean = 1, sd = 0.5)
      
              #...and correlation between n of (expected) episodes and duration
              co_ii[x, paste0(kk, "_dur")] <- co_ii[x, paste0(kk, "_dur")] + 
                (co_ii[x, paste0(kk, "_y")] - 
                mean(co_ii[, paste0(kk, "_y")], na.rm = T)) * b_corr   
               
          }
        }

        # If symptom time reaches duration of episode, episode ends
        x <- which(co_ii[, paste0(kk, "_time")] >= co_ii[, paste0(kk, "_dur")])
        co_ii[x, kk] <- F
        co_ii[x, paste0(kk, "_time")] <- NA
        co_ii[x, paste0(kk, "_dur")] <- 0
        co_ii[x, paste0(kk, "_wt")] <- 0
        
        # Lastly, advance symptom time
        co_ii[, paste0(kk, "_time")] <- co_ii[, paste0(kk, "_time")] + 1
      }
        
      #...................................      
      ## Model weight and height change
      for (kk in 1:nrow(co_ii)) {
        
        # Identify calibration parameter values
        x <- cal_par[which(cal_par$sex == co_ii[kk, "sex"] &
          cal_par$centile_weight == co_ii[kk, "centile_weight"] &
          cal_par$age == co_ii[kk, "age"]), c("par0", "par1")]
    
        # Update weight, FFM and FM for each child    
        co_ii[kk, c("ffm","fm","weight")] <- 
          f_hall_kids(df_f = co_ii[kk, ], par_f = x)
        
        # Update height
        co_ii[kk, "height"] <- centiles[which(centiles$sex == co_ii[kk, "sex"] &
        centiles$measure == "height" & centiles$age == co_ii[kk, "age"]), 
          co_ii[kk, "centile_height"]]      
      }
        
        # Update anthropometry
        x <- with(co_ii, anthro_zscores(sex = sex_anthro, age = age, 
          lenhei = height, weight = weight))
        co_ii[, colnames(x)] <- x
        
        # Identify incident MAM/SAM cases and start MAM/SAM time
        inc_mam_now <- which(!co_ii$mam & co_ii$zwfl < (-2) & co_ii$zwfl >=(-3))
        inc_sam_now <- which(!co_ii$sam & co_ii$zwfl < (-3))
        co_ii[inc_mam_now, "mam_time"] <- 0
        co_ii[inc_sam_now, "sam_time"] <- 0
        
        # Identify MAM/SAM cases who are no longer cases and end MAM/SAM time
        end_mam_now <- which(co_ii$mam & co_ii$zwfl >= (-2))
        end_sam_now <- which(co_ii$sam & co_ii$zwfl >= (-3))
        co_ii[end_mam_now, "mam_time"] <- NA
        co_ii[end_sam_now, "sam_time"] <- NA
        
        # Update MAM/SAM status
        co_ii[, c("mam", "sam")] <- F
        co_ii[which(co_ii$zwfl < (-2) & co_ii$zwfl >= (-3)), "mam"] <- T
        co_ii[which(co_ii$zwfl < (-3)), "sam"] <- T

      #...................................      
      ## Model MAM/SAM treatment
    
        # If a child treated for MAM deteriorates to SAM, end MAM treatment
            # (following code will give child a chance to start SAM treatment)
        x <- which(co_ii[inc_sam_now, "mam_tx"])
        co_ii[x, "mam_tx"] <- F
        co_ii[x, "therapeutic_food"] <- 0
        co_ii[x, c("mam_time", "sam_to_mam_time", "recovery_time")] <- NA
          
        # How many free MAM/SAM treatment slots are there? Attribute these
        for (kk in c("mam", "sam")) {
          # if treatment coverage is non-missing, sample from binomial distr.
          x <- 0
          if (! is.na(tl_ii[jj, paste0("tx_", kk, "_coverage")]) ){
            x <- max(rbinom(1, sum(co_ii[, kk]),
              tl_ii[jj, paste0("tx_", kk, "_coverage")]) - 
              sum(co_ii[, paste0(kk, "_tx")]), 0)
          }
          # alternatively, use daily n of admissions
          if (! is.na(tl_ii[jj, paste0("tx_", kk, "_admissions")]) ) {
            x <- tl_ii[jj, paste0("tx_", kk, "_admissions")]
          }
        
          # attribute slots to untreated cases that have made it past
              # treatment delay, and set dosage
          kk_other <- ifelse(kk == "mam", "sam", "mam")
          dosage <- ifelse(kk == "mam", mam_tx_daily, sam_tx_daily)
          if (x > 0) {
            
            # eligible for treatment
            xx <- which(co_ii[, kk] & !co_ii[, paste0(kk, "_tx")] & 
              !co_ii[, paste0(kk_other, "_tx")] & 
              co_ii[, paste0(kk, "_time")] >= co_ii[, paste0(kk, "_tx_delay")])
            
            # attribute slots
            if (length(xx) > 0) {
              if (length(xx) > x) {xx <- sample(xx, x)} # if need > tx available
              co_ii[xx, paste0(kk, "_tx")] <- T
              co_ii[xx, "therapeutic_food"] <- dosage
            }
          }
        }
          
        # If children have been on SAM treatment until now, but are now not SAM,
            # start SAM -> MAM recovery clock:
        x <- which(co_ii$sam_tx & ! co_ii$sam & is.na(co_ii$sam_to_mam_time))
        co_ii[x, "sam_to_mam_time"] <- 0
        
        # If children have been on M/SAM treatment until now, but are now OK,
            # start MAM -> OK recovery clock:
        x <- which((co_ii$sam_tx | co_ii$mam_tx) & ! co_ii$sam & ! co_ii$mam &
          is.na(co_ii$recovery_time))
        co_ii[x, "recovery_time"] <- 0
    
        # First-week recovery check:
            # has the child regressed to MAM/SAM? If so, reset recovery clock
        x <- which(co_ii$recovery_time == 7 & (co_ii$sam | co_ii$mam))
        co_ii[x, "recovery_time"] <- NA
        x <- which(co_ii$sam_to_mam_time == 7 & co_ii$sam)
        co_ii[x, "sam_to_mam_time"] <- NA 
        
        # Second-week recovery check:
          # has the child regressed to MAM/SAM? If so, reset recovery clock
          x <- which(co_ii$recovery_time == 14 & (co_ii$sam | co_ii$mam))
          co_ii[x, "recovery_time"] <- NA 
          x <- which(co_ii$sam_to_mam_time == 14 & co_ii$sam)
          co_ii[x, "sam_to_mam_time"] <- NA 
        
          # if improvement is sustained at second week check, 
            # discharge cured cases... 
            x <- which(co_ii$recovery_time == 14 & ! co_ii$sam & ! co_ii$mam)
            co_ii[x, c("recovery_time", "sam_to_mam_time")] <- NA
            co_ii[x, "therapeutic_food"] <- 0
            co_ii[x, c("mam_tx", "sam_tx")] <- F
            
            #...or reduce food dosage for SAM on treatment
            x <- which(co_ii$sam_to_mam_time == 14 & ! co_ii$sam)
            co_ii[x, "therapeutic_food"] <- mam_tx_daily 
            co_ii[x, "sam_to_mam_time"] <- NA
    
        # Advance time
        co_ii$mam_time <- co_ii$mam_time + 1
        co_ii$sam_time <- co_ii$sam_time + 1
        co_ii$recovery_time <- co_ii$recovery_time + 1
        co_ii$sam_to_mam_time <- co_ii$sam_to_mam_time + 1
        
    #### SMALL PROBLEMS: 
    # a) CHILDREN WHO DETERIORATE FROM MAM TO SAM ARE NOT GUARANTEED TREATMENT
    # b) KIDS KEEP EATING NORMAL INTAKE WHILE UNDER TREATMENT: SHOULD THEY?
    # c) SHOULD WE RESTRICT TREATMENT TO 6mo AND ABOVE?
            
      #...................................      
      ## Advance age; replace children who age out with newborns with same sex,
        # weight centile and formula-feeding characteristics
        
        # Advance age
        co_ii$age <- co_ii$age + 1
        co_ii$age_m <- floor(co_ii$age / 30.44)
        
        # Replace children who age out with newborns and reset some variables
        x <- which(co_ii$age >= 1826)
        if (length(x) > 0) {
          # reset time-varying variables
          co_ii[x, c("age", "age_m", "therapeutic_food")] <- 0
          co_ii[x, c("mam", "sam", "sam_tx", "mam_tx", "dia", "ari")] <- F 
          co_ii[x, c("mam_time", "sam_time", "recovery_time", "sam_to_mam_time",
            "ari_time", "dia_time")] <- NA
          co_ii[x, c("ari_dur", "ari_wt", "ari_y", 
            "dia_dur", "dia_wt", "dia_y")] <- 0
          
          # attribute birth weight, FM and FFM
          for (j in x) {
            co_ii[j,"weight"] <- centiles[which(centiles$sex == co_ii[j,"sex"] &
              centiles$measure== "weight" & centiles$age == 0), 
              co_ii[j, "centile_weight"]]
            co_ii[j, "fm"] <- f_bw_fm(co_ii[j, "sex"], co_ii[j, "weight"]) 
          }
          co_ii[x, "ffm"] <- co_ii[x, "weight"] - co_ii[x, "fm"]
          
          # select random height centiles for newborns from correlation model
          for (j in x) {
            xx <- corr_mat[which(corr_mat$sex == co_ii[j, "sex"] & 
              corr_mat$centile_weight == co_ii[j, "centile_weight"]), ]
            co_ii[j, "centile_height"]  <- sample(xx$centile_height, size = 1, 
              prob = xx$prob_fitted)
          }
        }
    
      #...................................      
      ## Tally outcomes for the time unit and return output
        
        # Outcomes for the area
        tl_ii[jj, c("n_dia", "n_ari")] <- colSums(co_ii[, c("dia", "ari")])
        tl_ii[jj, c("mam", "sam")] <- colSums(co_ii[, c("mam", "sam")])
        tl_ii[jj, "mam_inc"] <- length(inc_mam_now)
        tl_ii[jj, "sam_inc"] <- length(inc_sam_now)
        tl_ii[jj, "tx_mam_coverage"] <- sum(co_ii$mam & co_ii$mam_tx) /
          sum(co_ii$mam)
        tl_ii[jj, "tx_sam_coverage"] <- sum(co_ii$sam & co_ii$sam_tx) /
          sum(co_ii$sam)
        tl_ii[jj, "intake"] <- mean(co_ii$intake)
        
        x <- co_ii[which(co_ii$age_m >= 6), ]
        tl_ii[jj, "n_0659"] <- nrow(x)
        tl_ii[jj, "whz_mean"] <- mean(x$zwfl)
        tl_ii[jj, c("mam_0659", "sam_0659")] <- colSums(x[, c("mam", "sam")])
        tl_ii[jj, "mam_tx"] <- sum(x$mam_tx)
        tl_ii[jj, "sam_tx"] <- sum(x$sam_tx)
        
        x <- x[which(x$age_m < 24), ]
        tl_ii[jj, "n_0623"] <- nrow(x)
        tl_ii[jj, c("mam_0623", "sam_0623")] <- colSums(x[, c("mam", "sam")])
        
        x <- co_ii[which(co_ii$age_m < 6), ]
        tl_ii[jj, "formula_fed"] <- mean(x$formula_fed)
        tl_ii[jj, "n_u6mo"] <- nrow(x)
        
    } # <<<< (close time unit sub-loop jj)
    close(pb)
    
    #...................................      
    ## Add to overall output for the area
    tl_f[which(tl_f$area == ii), ] <- tl_ii
    co_f[which(co_f$area == ii), ] <- co_ii[, colnames(co_f)]
    
  } # << (close area loop ii)

  #...................................      
  ## Return output as list: timeline and cohort as of the last day of timeline
  out <- list(tl = tl_f, co = co_f)  
  return(out)
}  

 

#...............................................................................  
### ENDS
#...............................................................................
     