#...............................................................................
### + GENERIC MODEL OF ACUTE MALNUTRITION AS A FUNCTION OF CALORIC INTAKE ++ ###
#...............................................................................

#...............................................................................
## - R SCRIPT TO GENERATE A TIMELINE AND SCENARIOS OF CRISIS FACTOR VALUES -- ##
#...............................................................................


#...............................................................................  
### Generating an interpolated timeline of values over the retrospective period
#...............................................................................

if (retro == "yes") {
  
  #...................................      
  ## Preparatory steps

    # Establish timespan
      # if projection is desired but no date_start is specified, 
          # set it to 60 days before date_end
      if (project == "yes" & is.na(date_start)) {
        date_start <- as.Date(date_end - 60)
      }
      
      # work out maximum date for retrospective timeline
      x <- ifelse(is.na(date_start), date_end, date_start - 1)

    # Initialise retrospective timeline: one day for each area
    rtl <- expand.grid(date = as.Date(date_crisis : x), area = areas)
    
    # Add total population of under 5yos, if provided
    if (exists("pop")) {
      rtl <- merge(rtl, df_po, by = c("area", "date"), all.x = T)
      rtl <- rtl[order(rtl$area, rtl$date), ]
    }
    
    # Convert number of SAM/MAM admissions to daily mean during the date ranges
    x <- which(retro_data$indicator == "admissions")
    if (length(x) > 0) {
      retro_data[x, "value"] <- retro_data[x, "value"] / 
        as.numeric(retro_data[x, "date_2"] - retro_data[x, "date_1"] + 1)
    }
      
    # Identify exact variables (combination of factors and indicators)
    retro_data$variable <- paste(retro_data$factor,retro_data$indicator,sep="_")
    retro_data$variable <- gsub("_NA", "", retro_data$variable)
        
for (i in factors[["factor"]]) {
  #...................................      
  ## More preparatory steps
  
    # Select specific variables needed
    x <- i
    if(i %in% c("tx_sam", "tx_mam")) {
      x <- paste(i, c("coverage", "admissions"), sep = "_")
    }
    if (i %in% c("ari", "dia")) {
      x <- paste(i, c("period_prev", "point_prev", "incidence"), sep = "_")
    }
  
    # Initalise variables in timeline, including baseline values
    rtl[, x] <- NA
    rtl[, paste0(x, "_base")] <- 0
    if (any(x %in% c("req_ch_protected", "req_ad_protected"))) {
      rtl[, x] <- 1
      rtl[, paste0(x, "_base")] <- 1
    }

    # Select data needed for this factor
    df <- retro_data[which(retro_data$variable %in% x), 
      c("period", "area", "date_1", "date_2", "variable", "value")]
          
  #...................................      
  ## If any data have been entered...
  if (nrow(df) != 0) {

    # Identify baseline values for each protective/risk factor, if specified
      # any baseline values?
      base <- which(df$period == "baseline")
      
      # if yes, add them to the database
      if (length(base) != 0) {
        # set baseline values aside and compute their period-weighted mean
          # just in case multiple baseline sub-periods are specified
        df_base <- df[base, ]
        df_base$n_days <- as.numeric(df_base$date_2 - df_base$date_1 + 1)
        df_base[which(is.na(df_base$n_days)), "n_days"] <- 1
        
        # for each area...
        for (j in unique(df_base$area)) {
          for (k in x) {
            # take weighted mean baseline value, in case multiple baselines are 
                # entered for different sub-periods; weight = duration of period
            xx <- df_base[which(df_base$area == j & df_base$variable == k), ]
            xx <- weighted.mean(xx$value, xx$n_days/sum(xx$n_days))
            if (!is.na(xx)) {rtl[which(rtl$area == j), paste0(k,"_base")] <- xx}
          }
        }
      }
      
    # Add crisis values to timeline
    if (length(base) != 0) {df <- df[-base, ]}
    for (j in 1:nrow(df)) {
        var <- df[j, "variable"]
        rtl[which(rtl$area == df[j, "area"] & rtl$date %in% 
          as.Date(df[j, "date_1"]:df[j, "date_2"])), var] <- df[j, "value"]
      }
    }
    
  #...................................      
  ## Compute excess, and linearly interpolate all factor values
  
    # Excess (making sure the value is never negative)
    if (all(! x %in% c("req_ch_protected", "req_ad_protected"))) {
      for (j in x) {rtl[, j] <- pmax(rtl[, j] - rtl[, paste0(j, "_base")], 0)}
    }
    
    # Linear interpolation, by area, if there are at least two unique values
    rtl <- rtl[order(rtl$area, rtl$date), ]
    for (j in areas) {
      for (k in x) {
        xx <- na.omit(rtl[which(rtl$area == j), c("date", k)])
        if (length(unique(xx[, k])) >= 2) {
          rtl[which(rtl$area == j), k] <- suppressWarnings(approx(xx$date, xx[, k], 
            xout = as.Date(rtl[which(rtl$area == j), "date"]), rule = 2)$y)
        }
      }      
    }
}


  #...................................      
  ## Further prepare the dataset

    # Specify scenario as retrospective
    rtl$scenario <- "retrospective"
    
    # Scale SAM/MAM admissions to size of population being modelled
    rtl$tx_sam_admissions <- rtl$tx_sam_admissions * n_kids / rtl$pop_0059
    rtl$tx_mam_admissions <- rtl$tx_mam_admissions * n_kids / rtl$pop_0059

    # Establish equivalence of prevalence and incidence of ARI and diarrhoea
      # read equivalence table
      equi <- readRDS(paste0(dir_path, "out/14_disease_equi.rds"))
      equi$incidence <- equi$inc_shape

    indicators <- c("period_prev", "point_prev", "incidence")
    for (i in c("ari", "dia")) {
      
      # select equivalence table for this disease
      equi_i <- equi[which(equi$disease == i), ]
      
      # establish which indicator is non-missing during the timeline
      which_ind <- NA
      for (j in indicators) {
        if (! all(is.na(rtl[, paste(i, j, sep = "_")]))) {which_ind <- j} }
      
      # fill in the others based on equivalence  
      for (j in indicators[indicators != which_ind]) {
        var <- paste(i, j, sep = "_")
        rtl[, var] <- na.replace(approx(equi_i[, which_ind], equi_i[, j], 
          xout = rtl[, paste(i, which_ind, sep = "_")])$y, 0)
        rtl[, paste0(var, "_base")] <- na.replace(
          approx(equi_i[, which_ind], equi_i[, j], 
          xout = rtl[, paste(i, which_ind, "base", sep = "_")])$y, 0)
      }
    }

    
  #...................................      
  ## Visualise trends in factors
        
    # Prepare data for graphing
    df <- rtl
    x <- sapply(colnames(df), function(xx) {length(which(is.na(df[, xx])))==0})
    x <- names(x[x])
    x <- x[! x %in% c("date", "area", "req_ch_protected", "req_ad_protected",
      "scenario", grep("_base|incidence|point_prev|pop", x, value = T))]
    for (i in x) {df[, i] <- df[, i] + df[, paste0(i, "_base")]}
    df <- df[, c("area", "date", x)]
    df <- reshape(df, direction = "long", varying = x, 
      idvar = c("area", "date"), timevar = "variable", times=x, v.names="value")
    df$factor <- regmatches(df$variable, regexpr(
      paste(factors[["factor"]], collapse = "|"), df$variable))
    r <- regexpr(paste(c(indicators, "coverage", "admissions"), collapse = "|"), 
      df$variable)
    df$indicator <- sapply(regmatches(df$variable, r, invert=NA), `[`, 2)
    df <- merge(df, factors, by = "factor", all.x = T)
    x <- which(df$indicator %in% c("admissions", "coverage"))
    df[x, "factor_name"] <- paste(df[x, "factor_name"], df[x, "indicator"], 
      sep = " ")
    
    # Plot
    pl <- ggplot(df, aes(x = date, y = value, colour = area, 
      linetype = factor_name)) +
      geom_step(linewidth = 1) +
      scale_x_date("date", date_labels = "%b-%Y", breaks = "1 month", 
        expand = c(0,0)) +
      scale_y_continuous("interpolated value", 
        expand = expansion(mult = c(0, 0.1), add = 0), 
        limits = c(0, NA)) +
      scale_colour_manual(values = colours_areas) +
      scale_linetype_discrete() +
      theme_bw() +
      theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      facet_grid(factor_name ~ area, scales = "free_y")
    ggsave(paste0(dir_path, "out/21_timeline_factors_retro.png"), dpi = "print",
      units = "cm", height = 28, width = (length(areas)) * 10)

    # Add/remove additional columns, save
      # remove baseline values
      rtl <- rtl[, ! colnames(rtl) %in% grep("_base", colnames(rtl), value = T)]
      rtl[, outcomes] <- NA
    
      # save
      saveRDS(rtl, paste0(dir_path, "out/21_rtl.rds"))
}
      
      
#...............................................................................
### Generating a list of projection scenarios
#...............................................................................

if (project == "yes") {
  
  #...................................
  ## Harmonise scenarios dataset

    # Convert disease indicators to point prevalence, if needed
    equi <- readRDS(paste0(dir_path, "out/14_disease_equi.rds"))
    for (i in c("dia", "ari")) {
      
      # if period prevalence...
      x <- which(scenarios$factor == i & scenarios$indicator == "period_prev" & 
        ! is.na(scenarios$abs_value))
      if (length(x) > 0) {
        scenarios[x, "abs_value"] <- 
          na.replace(approx(equi[which(equi$disease == i), "period_prev"], 
            equi[which(equi$disease == i), "point_prev"], 
            xout = scenarios[x, "abs_value"])$y, 0)
      }
      
      # if incidence...
      x <- which(scenarios$factor == i & scenarios$indicator == "incidence" & 
        ! is.na(scenarios$abs_value))
      if (length(x) > 0) {
        scenarios[x, "abs_value"] <- 
          na.replace(approx(equi[which(equi$disease == i), "inc_shape"], 
            equi[which(equi$disease == i), "point_prev"], 
            xout = scenarios[x, "abs_value"])$y, 0)
      }
    }
    scenarios$indicator <- gsub("period_prev","point_prev", scenarios$indicator)
    scenarios$indicator <- gsub("incidence", "point_prev", scenarios$indicator)

    # If relative risk is specified but reference period is missing, assume 30d
    x <- which(! is.na(scenarios$rel_risk) & is.na(scenarios$days_ref))
    scenarios[x, "days_ref"] <- 30

    # Add any factors that do not feature in one or more scenarios/areas
        # assume mean of the previous 30 days (i.e. rel_risk = 1, days_ref = 30)
        # only add treatment coverage if neither coverage nor 
        # treatment admissions are already included; never add admissions
    for (i in unique(scenarios$scenario)) {
      for (j in areas) {
        
        # select observations
        df <- scenarios[which(scenarios$scenario == i & scenarios$area == j), ]
        
        # which factors are missing
        x <- base::setdiff(factors[["factor"]], unique(df$factor))
        
        # add
        for (k in x) {
          scenarios[nrow(scenarios) + 1, ] <- 
            data.frame(i, k, NA, j, date_start, date_end, NA, 1, 30)
        }
      }
    }
    x <- which(scenarios$factor %in% c("tx_sam", "tx_mam") 
      & is.na(scenarios$indicator))
    scenarios[x, "indicator"] <- "coverage"
    x <-which(scenarios$factor %in% c("ari","dia") & is.na(scenarios$indicator))
    scenarios[x, "indicator"] <- "point_prev"
        
    # Identify exact variables (combination of factors and indicators)
    scenarios$variable <- paste(scenarios$factor, scenarios$indicator,sep="_")
    scenarios$variable <- gsub("_NA", "", scenarios$variable)
  
    
  #...................................
  ## Initialise and save list of scenarios

    # Create list, each component of which contains a scenario
    x <- unique(scenarios$scenario)
    scens <- vector("list", length = length(x))
    names(scens) <- unique(scenarios$scenario)
    for (i in names(scens)) {
      scens[[i]] <- scenarios[which(scenarios$scenario == i), ]
    }

    # Save
    saveRDS(scens, paste0(dir_path, "out/21_scens.rds"))
}
      

#...............................................................................  
### ENDS
#...............................................................................