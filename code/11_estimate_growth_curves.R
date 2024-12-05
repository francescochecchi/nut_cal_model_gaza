#...............................................................................
### + GENERIC MODEL OF ACUTE MALNUTRITION AS A FUNCTION OF CALORIC INTAKE ++ ###
#...............................................................................

#...............................................................................
## ------ R SCRIPT TO ESTIMATE GROWTH CURVES FROM ANTHROPOMETRIC DATA  ------ ##
   # will need to be adapted for each crisis unless we can generalise the script
#...............................................................................

if (! file.exists(paste0(dir_path, "out/11_centiles.rds")))
{
  
#...............................................................................  
### Preparing growth monitoring dataset to estimate growth curves
#...............................................................................

  #...................................      
  ## Read individual data from 2014 [g]rowth mo[n]itoring of children in Gaza
    
    # Read data
    df_gn <- data.frame(read_rds(paste0(dir_path, "in/gm_anthro_ind.rds")))
    
    # Track attrition
    attrition <- data.frame(group = c("all observations", "non-missing",
      ">= 3 observations", "no flags"), number = NA, percentage = NA)
    attrition[which(attrition$group=="all observations"),"number"] <-nrow(df_gn)

    # Explore missingness
    for (i in colnames(df_gn)) {
      print(i)
      print(prop.table(table(is.na(df_gn[, i]))))
    }
      # all complete
    attrition[which(attrition$group == "non-missing"), "number"] <-
      nrow(na.omit(df_gn))
    
    # Fix some variables
    df_gn$sex <- ifelse(df_gn$sex == "F", "female", "male")
    df_gn$age <- as.integer(df_gn$date - df_gn$dob)

  #...................................      
  ## Describe and select individual growth monitoring dataset
  
    # How many observations, by sex
    nrow(df_gn)
    table(df_gn$sex)
    
    # How many observations per child
      
      # tabulate
      df_gn$obs <- 1
      x <- aggregate(list(n_obs = df_gn$obs),by = list(id = df_gn$id), FUN= sum)
      x <- table(x$n_obs)
      x <- data.frame(n_obs = as.integer(names(x)), n_children = as.vector(x) )
      
      # plot
      plot1 <- ggplot(x, aes(x = n_obs, y = n_children)) +
        geom_bar(fill = palette_gen[10], alpha = 0.75, stat = "identity",
          colour = palette_gen[10]) +
        theme_bw() +
        scale_y_continuous("number of children", labels = comma) +
        scale_x_continuous("number of measurements", limits = c(0, 30),
          breaks = seq(0, 100, 5))
    
      # how many < 3 and >= 30
      x$prop <- x$n_children / sum(x$n_children)
      sum(x[which(x$n_obs < 3), "n_children"])
      sum(x[which(x$n_obs < 3), "n_children"]) / sum(x$n_children)
      sum(x[which(x$n_obs >= 30), "n_children"])
      sum(x[which(x$n_obs >= 30), "n_children"]) / sum(x$n_children)
      
    # Plot timing of monitoring visits
    plot2 <- ggplot(df_gn, aes(x = age / 30.44)) +
      geom_histogram(binwidth = 1, alpha = 0.5, fill = palette_gen[7],
        colour = palette_gen[7]) +
      theme_bw() +
      scale_x_continuous("age (months)", breaks = seq(0, 60, 3), 
        expand = c(0,0)) +
      scale_y_continuous("number of measurements", 
        breaks = seq(0, 1000000, 20000), expand = c(0,5000), labels = comma)
    
    # Combo plot
    ggarrange(plot1, plot2, ncol = 1, nrow = 2, labels = c("A", "B"),
      font.label = list(size = 10))
    ggsave(paste0(dir_path, "out/11_n_obs_combi.png"),
      dpi = "print", units = "cm", height = 20, width = 25)    
    
    # Exclude children with < 3 observations
    x <- aggregate(list(n_obs = df_gn$obs),by = list(id = df_gn$id), FUN= sum)  
    df_gn <- merge(df_gn, x, by = "id", all.x = T)
    attrition[which(attrition$group == ">= 3 observations"), "number"] <-
      length(which(df_gn$n_obs >= 3))
    df_gn <- subset(df_gn, n_obs >= 3)
    
    
  #...................................      
  ## Anonymise data and eliminate outliers (WAZ, WFH flags)

    # Anonymise IDs
    x <- data.frame(id = sort(unique(df_gn$id)))
    x$id_new <- paste0("c", 1:length(x$id))
    df <- merge(df_gn, x, by = "id", all.x = T)
    df <- subset(df, select = -id)
    colnames(df)[colnames(df) == "id_new"] <- "id"
    df$id <- factor(df$id)
        
    # Eliminate flags
    df$sex_anthro <- ifelse(df$sex == "female", "f", "m")
    x <- with(df, anthro::anthro_zscores(sex = sex_anthro, weight = weight,
      lenhei = height, age = age, is_age_in_month = F))
    table(x$fwei, x$flen)
    df <- cbind(df, x)
    attrition[which(attrition$group == "no flags"), "number"] <-
      length(which(df$fwei == 0 & df$flen == 0))
    df <- subset(df, fwei == 0 & flen == 0)
    
    # Check that all children still have >=3 observations
    table(df$n_obs, df$sex)

    # Figure out ranges of height and weight (for plots)
    height_lims <- quantile(df$height, c(0.01, 0.99))
    weight_lims <- quantile(df$weight, c(0.01, 0.99))
    
    # Remove uncleaned individual database
    rm(df_gn)
    
    # Save attrition
    attrition$percentage <- scales::percent(
      attrition$number / attrition[1, "number"], 0.1)
    write.csv(attrition, paste0(dir_path, "out/11_attrition.csv"), row.names=F)

            
#...............................................................................  
### Estimating height for age curves
#...............................................................................
          
  #...................................      
  ## Estimate height by age for girls
    
    # Fit growth model (height as outcome)
    x <- na.omit(df[which(df$sex == "female"), c("age", "height")])
    fit_hf <- lms(y = height, x = age, data = x, families = c("BCPEo"), k = 4)

    # Compute median and 95% percentiles at each age
    out_hf <- centiles.pred(fit_hf, cent = c(2.5, 25, 50, 75, 97.5), 
      xvalues = 0:as.integer(5*365.25), xname = "age")
    colnames(out_hf) <- c("age", "lci025", "lci250", 
      "median", "uci750", "uci975")
    
    # Visualise growth curve
    plot_hf_curve <- ggplot() +
      geom_line(data = subset(df, sex == "female"), 
        aes(x = age / 30.44, y = height, group = id), alpha = 0.1) +
      geom_line(data = out_hf, aes(x = age / 30.44, y = median), 
        colour = palette_gen[9], linewidth = 1.5, linetype = "21") +
      geom_ribbon(data = out_hf, 
        aes(x = age / 30.44, ymin = lci025, ymax = uci975),
        fill = palette_gen[9], colour = NA, alpha = 0.15) +
      geom_ribbon(data = out_hf, 
        aes(x = age / 30.44, ymin = lci250, ymax = uci750),
        fill = palette_gen[9], colour = NA, alpha = 0.3) +
      geom_line(data = subset(who_haz_standards, sex == "female"), 
        aes(y = SD0, x = age / 30.44), colour = palette_gen[15], linewidth = 1,
        alpha = 0.5) +
      theme_bw() +
      scale_y_continuous("height (cm)", breaks = seq(0, 200, 10),
        limits = c(40, 120)) +
      scale_x_continuous("age (months)", breaks = seq(0, 60, 5), 
        expand = c(0, 0))

  #...................................      
  ## Estimate height by age for boys
    
    # Fit growth model (height as outcome)
    x <- na.omit(df[which(df$sex == "male"), c("age", "height")])
    fit_hm <- lms(y = height, x = age, data = x, families = c("BCPEo"), k = 4)
    
    # Compute median and 95% percentiles at each age
    out_hm <- centiles.pred(fit_hm, cent = c(2.5, 25, 50, 75, 97.5), 
      xvalues = 0:as.integer(5*365.25), xname = "age")
    colnames(out_hm) <- c("age", "lci025", "lci250", 
      "median", "uci750", "uci975")
    
    # Visualise growth curve
    plot_hm_curve <- ggplot() +
      geom_line(data = subset(df, sex == "male"), 
        aes(x = age / 30.44, y = height, group = id), alpha = 0.1) +
      geom_line(data = out_hm, aes(x = age / 30.44, y = median), 
        colour = palette_gen[3], linewidth = 1.5, linetype = "21") +
      geom_ribbon(data = out_hm, 
        aes(x = age / 30.44, ymin = lci025, ymax = uci975),
        fill = palette_gen[3], colour = NA, alpha = 0.15) +
      geom_ribbon(data = out_hm, 
        aes(x = age / 30.44, ymin = lci250, ymax = uci750),
        fill = palette_gen[3], colour = NA, alpha = 0.3) +
      geom_line(data = subset(who_haz_standards, sex == "male"), 
        aes(y = SD0, x = age / 30.44), colour = palette_gen[15], linewidth = 1,
        alpha = 0.5) +
      theme_bw() +
      scale_y_continuous("height (cm)", breaks = seq(0, 200, 10),
        limits = c(40, 120)) +
      scale_x_continuous("age (months)", breaks = seq(0, 60, 5), 
        expand = c(0, 0))
    
  #...................................      
  ## Combined graph
    
    # Combined graph
    ggarrange(
      plot_hf_curve + theme(plot.margin = margin(1,0.1,0.1,0.1, "cm")), 
      plot_hm_curve + theme(plot.margin = margin(1,0.1,0.1,0.1, "cm")), 
      labels = c("girls", "boys"), ncol = 2, nrow = 1, 
      align = "hv", font.label = list(size = 11))
    ggsave(paste0(dir_path, "out/11_height_curves_combi.png"),
      dpi = "print", units = "cm", height = 15, width = 30)    
    
    

#...............................................................................  
### Estimating weight for age curves
#...............................................................................
          
  #...................................      
  ## Estimate weight by age for girls
    
    # Fit growth model (weight as outcome)
    x <- na.omit(df[which(df$sex == "female"), c("age", "weight")])
    fit_wf <- lms(y = weight, x = age, data = x, families = c("BCPEo"), k = 4)
    
    # Compute median and 95% percentiles at each age
    out_wf <- centiles.pred(fit_wf, cent = c(2.5, 25, 50, 75, 97.5), 
      xvalues = 0:as.integer(5*365.25), xname = "age")
    colnames(out_wf) <- c("age", "lci025", "lci250", 
      "median", "uci750", "uci975")
    
    # Visualise growth curve
    plot_wf_curve <- ggplot() +
      geom_line(data = subset(df, sex == "female"), 
        aes(x = age / 30.44, y = weight, group = id), alpha = 0.1) +
      geom_line(data = out_wf, aes(x = age / 30.44, y = median), 
        colour = palette_gen[9], linewidth = 1.5, linetype = "21") +
      geom_ribbon(data = out_wf, 
        aes(x = age / 30.44, ymin = lci025, ymax = uci975),
        fill = palette_gen[9], colour = NA, alpha = 0.15) +
      geom_ribbon(data = out_wf, 
        aes(x = age / 30.44, ymin = lci250, ymax = uci750),
        fill = palette_gen[9], colour = NA, alpha = 0.3) +
      geom_line(data = subset(who_waz_standards, sex == "female"), 
        aes(y = SD0, x = age / 30.44), colour = palette_gen[15], linewidth = 1,
        alpha = 0.5) +
      theme_bw() +
      scale_y_continuous("weight (cm)", breaks = seq(0, 24, 2),
        limits = c(0, 24)) +
      scale_x_continuous("age (months)", breaks = seq(0, 60, 5), 
        expand = c(0, 0))
    
  #...................................      
  ## Estimate weight by age for boys
  
    # Fit growth model (weight as outcome)
    x <- na.omit(df[which(df$sex == "male"), c("age", "weight")])
    fit_wm <- lms(y = weight, x = age, data = x, families = c("BCPEo"), k = 4)
    
    # Compute median and 95% percentiles at each age
    out_wm <- centiles.pred(fit_wm, cent = c(2.5, 25, 50, 75, 97.5), 
      xvalues = 0:as.integer(5*365.25), xname = "age")
    colnames(out_wm) <- c("age", "lci025", "lci250", 
      "median", "uci750", "uci975")
    
    # Visualise growth curve
    plot_wm_curve <- ggplot() +
      geom_line(data = subset(df, sex == "male"), 
        aes(x = age / 30.44, y = weight, group = id), alpha = 0.1) +
      geom_line(data = out_wm, aes(x = age / 30.44, y = median), 
        colour = palette_gen[3], linewidth = 1.5, linetype = "21") +
      geom_ribbon(data = out_wm, 
        aes(x = age / 30.44, ymin = lci025, ymax = uci975),
        fill = palette_gen[3], colour = NA, alpha = 0.15) +
      geom_ribbon(data = out_wm, 
        aes(x = age / 30.44, ymin = lci250, ymax = uci750),
        fill = palette_gen[3], colour = NA, alpha = 0.3) +
      geom_line(data = subset(who_waz_standards, sex == "male"), 
        aes(y = SD0, x = age / 30.44), colour = palette_gen[15], linewidth = 1,
        alpha = 0.5) +
      theme_bw() +
      scale_y_continuous("weight (cm)", breaks = seq(0, 24, 2),
        limits = c(0, 24)) +
      scale_x_continuous("age (months)", breaks = seq(0, 60, 5), 
        expand = c(0, 0))
    
  #...................................      
  ## Combined graph
    
    # Combined graph
    ggarrange(
      plot_wf_curve + theme(plot.margin = margin(1,0.1,0.1,0.1, "cm")), 
      plot_wm_curve + theme(plot.margin = margin(1,0.1,0.1,0.1, "cm")), 
      labels = c("girls", "boys"), ncol = 2, nrow = 1, 
      align = "hv", font.label = list(size = 11))
    ggsave(paste0(dir_path, "out/11_weight_curves_combi.png"),
      dpi = "print", units = "cm", height = 12, width = 25)    
    
    
#...............................................................................  
### Computing growth percentiles and height distributions by age, sex and weight
#...............................................................................

  #...................................      
  ## Compute growth centiles for each sex and weight / height
    
    # Compute centiles
    out_centiles <- data.frame()
    for (i in c("h", "w")) {
      for (j in c("f", "m")) {
        
        # get fit
        fit <- get(paste0("fit_", i, j))
        
        # predict centiles
        x <- centiles.pred(fit, cent = c(0.1, 1:99, 99.9),
          xvalues =  0:as.integer(5*365.25), xname = "age")
        
        # add to output
        x$sex <- ifelse(j == "f", "female", "male")
        x$measure <- ifelse(i == "h", "height", "weight")
        out_centiles <- rbind(out_centiles, x)
      }
    }
    
    # Rename columns
    colnames(out_centiles) <- c("age", paste0("c", 0:100), "sex", "measure")
    
    
  # #...................................      
  # ## Collect empirical distribution of height by sex, age and weight centile
  #     # needed for model simulation; age in months, weight deciles
  # 
  #   # Generate needed variables
  #   df$age_m <- floor(df$age / 30.44)
  #   df$centile_weight_cat <- cut(df$centile_weight, breaks = seq(0, 100, 10),
  #     include.lowest = T, right = F)
  #   
  #   # Collect height vectors for each stratum
  #   height_vecs <- aggregate(list(height = df$height), 
  #     by = df[, c("sex", "age_m", "centile_weight_cat")], FUN = as.vector)
    
 
#...............................................................................
### Investigating height-weight correlation
#...............................................................................

  #...................................
  ## Identify nearest weight and height centiles over each child's life course

    # Output: height and weight centile positions along age
    df$centile_height <- NA
    df$centile_weight <- NA

    # For each child, figure out nearest centile at each measurement age
    for (i in 1:nrow(df)) {

      # find nearest height centile
      df[i, "centile_height"] <- findInterval(df[i, "height"],
        out_centiles[which(out_centiles$measure == "height" &
          out_centiles$sex == df[i, "sex"] & out_centiles$age == df[i, "age"]),
          paste0("c", 0:100)], all.inside = T)

      # find nearest weight centile
      df[i, "centile_weight"] <- findInterval(df[i, "weight"],
        out_centiles[which(out_centiles$measure == "weight" &
            out_centiles$sex == df[i, "sex"] & out_centiles$age ==df[i, "age"]),
          paste0("c", 0:100)], all.inside = T)
    }
    df$centile_height <- df$centile_height - 1
    df$centile_weight <- df$centile_weight - 1


  #...................................
  ## Visualise correlation between height and weight centiles

    # Compute correlation coefficient (Spearman's rho), by 3-month period in age
    df$age_cat <- cut(df$age, breaks = seq(0, 2000, by = 90))
    table(df$age_cat)
    x <- by(df, df[, c("age_cat", "sex")],
      function(x) {data.frame(
        age_cat = unique(x$age_cat),
        sex = unique(x$sex),
        srho = cor(x$centile_height, x$centile_weight, method = "spearman")
      )})
    x <- do.call(rbind, x)

    # Visualise correlation by age
    ggplot(x, aes(x = age_cat, y = srho, colour = sex, group = sex)) +
      geom_point() +
      geom_line(alpha = 0.5) +
      scale_colour_manual(values = palette_gen[c(9,3)]) +
      theme_bw() +
      theme(legend.position = "top") +
      scale_y_continuous("Spearman's correlation coefficient") +
      scale_x_discrete("age (days)", labels = seq(90, 2000, by = 90))
    ggsave(paste0(dir_path, "out/11_corr_by_age.png"),
      dpi = "print", units = "cm", height = 12, width = 25)

    # Compute correlation matrix, with each row = one combination of centiles
    x <- df
    x <- aggregate(x[, c("centile_height", "centile_weight")],
      by = x[, c("id", "sex")], mean)
    x[, c("centile_height", "centile_weight")] <-
      round(x[, c("centile_height", "centile_weight")], digits = 0)
    x$n <- 1
    x <- aggregate(list(n = x$n),
      by = x[, c("centile_height", "centile_weight", "sex")], sum)
    x1 <- expand.grid(centile_height = 0:99, centile_weight = 0:99,
      sex = c("female", "male"))
    x <- merge(x1, x, by = c("centile_height", "centile_weight", "sex"),
      all.x = T)
    x$n <- na.replace(x$n, 0)
    x1 <- aggregate(list(n_bin = x$n), by = x[, c("centile_weight", "sex")], 
      FUN = sum)
    x <- merge(x, x1, by = c("centile_weight", "sex"), all.x = T)
    x$prob <- x$n / x$n_bin

    # Plot correlation matrix
    ggplot(x, aes(x = centile_height, y = centile_weight, fill = n)) +
      geom_tile(colour = "black") +
      scale_fill_viridis_c() +
      scale_x_continuous("mean height percentile during growth (%)",
        breaks = seq(0, 100, 10), expand = c(0, 0)) +
      scale_y_continuous("mean weight percentile during growth (%)",
        breaks = seq(0, 100, 10), expand = c(0, 0)) +
      theme_bw() +
      labs(fill = "frequency")
    ggsave(paste0(dir_path, "out/11_weight_height_corr.png"),
      dpi = "print", units = "cm", height = 20, width = 23)

    
  #...................................
  ## Model correlation matrix between height and weight centiles

    # Fit zero-inflated neg-bin model with tensor smoothing product
    fit_corr <- mgcv::bam(n ~ sex + te(centile_height, centile_weight) + 
        offset(log(n_bin)), data = x, family = nb())
    summary(fit_corr)
    x$fitted <- predict(fit_corr, newdata = x, type = "response")
    x1 <- aggregate(list(n_bin_fitted = x$fitted), 
      by = x[, c("centile_weight", "sex")], FUN = sum)
    x <- merge(x, x1, by = c("centile_weight", "sex"), all.x = T)
    x$prob_fitted <- x$fitted / x$n_bin_fitted
    x <- x[order(x$sex, x$centile_weight, x$centile_height), ]

    # Plot fitted values
    ggplot(x, aes(x = centile_height, y = centile_weight, fill = fitted)) +
      geom_tile(colour = "black") +
      scale_fill_viridis_c() +
      scale_x_continuous("mean height percentile during growth (%)",
        breaks = seq(0, 100, 10), expand = c(0, 0)) +
      scale_y_continuous("mean weight percentile during growth (%)",
        breaks = seq(0, 100, 10), expand = c(0, 0)) +
      theme_bw() +
      labs(fill = "frequency")
    ggsave(paste0(dir_path, "out/11_weight_height_corr_fitted.png"),
      dpi = "print", units = "cm", height = 20, width = 23)



#...............................................................................  
### Saving data and model fits needed for other scripts
#...............................................................................
    
  #...................................      
  ## Save model fits
    
    # # Save height model fits
    # saveRDS(fit_hf, paste0(dir_path, "out/11_fit_hf.rds"))
    # saveRDS(fit_hm, paste0(dir_path, "out/11_fit_hm.rds"))
    # 
    # # Save weight model fits
    # saveRDS(fit_wf, paste0(dir_path, "out/11_fit_wf.rds"))
    # saveRDS(fit_wm, paste0(dir_path, "out/11_fit_wm.rds"))
    
    # Save weight-height correlation fit
    saveRDS(fit_corr, paste0(dir_path, "out/11_fit_corr.rds"))
    
  #...................................      
  ## Save other output

    # Growth percentiles
    saveRDS(out_centiles, paste0(dir_path, "out/11_centiles.rds"))
    
    # # Height vectors, by sex, age and weight
    # saveRDS(height_vecs, paste0(dir_path, "out/11_height_vectors.rds"))
    
    # Correlation matrix
    x$centile_height <- paste0("c", x$centile_height)
    x$centile_weight <- paste0("c", x$centile_weight)
    x$sex <- as.character(x$sex)
    saveRDS(x, paste0(dir_path, "out/11_ht_wt_correlation.rds"))
      
    # Individual growth curves
    x <- c("id", "sex", "n_obs", "age", "weight", "height",
      "centile_height", "centile_weight")
    saveRDS(df[, x],paste0(dir_path, "out/11_ind_growth_curves.rds"))
    
}


 
#...............................................................................  
### Comparing crisis-specific growth curves to WHO standards (weight-for-age)
#...............................................................................

  #...................................      
  ## Compare the standard deviation of the GAMLSS centiles with WHO 2006
  
    # Read growth centiles from GAMLSS model
    out_centiles <- readRDS(paste0(dir_path,"out/11_centiles.rds"))

    # Approximate SD of GAMLSS centiles
    out_centiles$sd <- (abs(out_centiles[, "c50"] - out_centiles[, "c16"]) + 
        abs(out_centiles[, "c50"] - out_centiles[, "c84"])) / 2     

    # Compare HAZ
    who_haz_standards$sd_who <- 
      (abs(who_haz_standards$SD0-who_haz_standards$SD1neg) 
      + abs(who_haz_standards$SD0-who_haz_standards$SD1)) / 2 
    df <- merge(out_centiles[which(out_centiles$measure == "height"), 
      c("age", "sex", "sd")], who_haz_standards[, c("age", "sex", "sd_who")],
      by = c("age", "sex"), all.x = T)
    df <- reshape(df, direction = "long", varying = c("sd", "sd_who"),
      timevar = "source", times = c("GAMLSS fit", "WHO 2006"), 
      idvar = c("age", "sex"), v.names = "sd")
    ggplot(df, aes(x = age/30.44, y = sd, linetype = source, colour = sex)) +
      geom_line() +
      facet_wrap(sex ~., nrow = 2) +
      theme_bw() +
      scale_linetype_manual(values = c("22", "solid")) +
      scale_colour_manual(values = palette_gen[c(9, 3)]) +
      theme(legend.position = "top") +
      scale_y_continuous("standard deviation of height-for-age") +
      scale_x_continuous("age (months)", breaks = seq(0, 60, 6)) +
      guides(colour = "none")
    ggsave(paste0(dir_path, "out/13_haz_sd_gamlss_vs_who.png"),
      dpi = "print", units = "cm", height = 10, width = 22)  
    
    # Compare WAZ
    who_waz_standards$sd_who <- 
      (abs(who_waz_standards$SD0-who_waz_standards$SD1neg) 
        + abs(who_waz_standards$SD0-who_waz_standards$SD1)) / 2 
    df <- merge(out_centiles[which(out_centiles$measure == "weight"), 
      c("age", "sex", "sd")], who_waz_standards[, c("age", "sex", "sd_who")],
      by = c("age", "sex"), all.x = T)
    df <- reshape(df, direction = "long", varying = c("sd", "sd_who"),
      timevar = "source", times = c("GAMLSS fit", "WHO 2006"), 
      idvar = c("age", "sex"), v.names = "sd")
    ggplot(df, aes(x = age/30.44, y = sd, linetype = source, colour = sex)) +
      geom_line() +
      facet_wrap(sex ~., nrow = 2) +
      theme_bw() +
      scale_linetype_manual(values = c("22", "solid")) +
      scale_colour_manual(values = palette_gen[c(9, 3)]) +
      theme(legend.position = "top") +
      scale_y_continuous("standard deviation of weight-for-age") +
      scale_x_continuous("age (months)", breaks = seq(0, 60, 6)) +
      guides(colour = "none")
    ggsave(paste0(dir_path, "out/13_waz_sd_gamlss_vs_who.png"),
      dpi = "print", units = "cm", height = 10, width = 22)  
    


#...............................................................................  
### ENDS
#...............................................................................
  