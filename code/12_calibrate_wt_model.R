#...............................................................................
### + GENERIC MODEL OF ACUTE MALNUTRITION AS A FUNCTION OF CALORIC INTAKE ++ ###
#...............................................................................

#...............................................................................
## ----- R SCRIPT TO CALIBRATE WEIGHT MODEL TO REPRODUCE GROWTH CURVES  ----- ##
#...............................................................................


if (! file.exists(paste0(dir_path, "out/12_df_ha.csv"))) {
  
#...............................................................................  
### Calibrating Kcal intake to reproduce median weight growth
#...............................................................................

  #...................................      
  ## Read or generate necessary inputs
  
    # Read growth centiles from GAMLSS model
    out_centiles <- readRDS(paste0(dir_path, "out/11_centiles.rds"))    
    
    # Select only median weight curves
    medians <- out_centiles[which(out_centiles$measure == "weight"), 
      c("age", "sex", "c50")]
    medians <- medians[order(medians$sex, medians$age), ]
  
    # Blank set of individual calibration parameters (see below)
    par_blank <- c(
      "par0" = 0, # minimum adjustment over lifespan
      "par1" = 0 # sigma (standard deviation of log-normal function)
    )
    
    
  #...................................      
  ## Estimate actual Kcal intake by brute force
    
    # Candidate values of ratio of actual to reference Kcal to explore over
    cand <- seq(0.80, 1.20, 0.01)
    
    # Initialise output
    out_act <- expand.grid(cand = cand, sex = c("female", "male"))
    out_act$rmse <- NA
    
    # For each candidate value, and for each sex, compute RMSE of predictions
        # versus medians
    for (i in cand) {
      # set candidate actual Kcal intake
      df_ha_i <- df_ha
      df_ha_i$intake_act <- df_ha_i$intake_ref * i
        
      for (j in c("female", "male")) {
        
        # select dataset
        df_j <- subset(medians, sex == j)
        df_j$weight <- df_j$c50
        df_j$formula_fed <- F
        df_j$centile_weight <- "c50"
        df_j$dia <- F
        df_j$ari <- F
        
        # generate weight predictions
        x <- f_hall_pred(data_f = df_j, df_ha_f = df_ha_i, intake_f = "actual",
          par_f = par_blank)
        
        # compute RMSE of predictions vs actual medians
        out_act[which(out_act$cand == i & out_act$sex == j), "rmse"] <-
          sqrt(mean((df_j$weight - x$weight)^2))
      }  
    }
    
    # Estimate actual intake based on lowest RMSE
    cand <- by(out_act, out_act$sex, function(x) {x$cand[which.min(x$rmse)]})
    cand <- data.frame(sex = c("female", "male"), act = as.vector(unlist(cand)))
    df_ha <- merge(df_ha, cand, by = "sex", all.x = T)
    df_ha$intake_act <- df_ha$intake_ref * df_ha$act
    
    # Plot RMSE  
    ggplot(out_act, aes(x = cand, y = rmse, colour = sex, shape = sex)) +
      geom_point(alpha = 0.5) +
      geom_line() +
      scale_colour_manual(values = palette_gen[c(9,3)]) +
      scale_shape_manual(values = c(1,2)) +
      theme_bw() +
      theme(legend.position = "inside", legend.position.inside = c(0.5, 0.8)) +
      scale_y_continuous("root mean square error") +
      scale_x_continuous("ratio of actual to recommended intake", 
        breaks = seq(0.8, 1.2, 0.05))
    
    ggsave(paste0(dir_path, "out/12_intake_act_rmse.png"),
      dpi = "print", units = "cm", height = 8, width = 20)               
    
        
  #...................................      
  ## Visualise model fit vs. median curves with/without Kcal intake adjustment
  
    # Generate unadjusted and adjusted plots
    for (i in c("reference", "actual") ) {
      
      # for girls and boys...
      for (j in c("female", "male")) {
        
        # select dataset
        df_j <- subset(medians, sex == j)
        df_j$weight <- df_j$c50
        df_j$centile_weight <- "c50"
        df_j$formula_fed <- F
        df_j$dia <- F
        df_j$ari <- F
        
        # generate weight predictions
        x <- f_hall_pred(data_f = df_j, intake_f = i, par_f = par_blank)
        
        # collect predictions
        medians[which(medians$sex == j), "hall"] <- x$weight
      }  
      
      # plot model fit
      x <- reshape(medians, direction = "long", varying = c("c50", "hall"),
        idvar = c("age", "sex"), timevar = "source", 
        times = c("median", "weight model"), v.names = "weight")
    
      plot <- ggplot(x, aes(x = age / 30.44, y = weight, colour = sex, 
        linetype = source)) +
        geom_line(alpha = 0.5, linewidth = 1.5) +
        theme_bw() +
        theme(legend.position = "none") +
        facet_wrap(. ~ sex, ncol = 2) +
        scale_colour_manual(values = palette_gen[c(9,3)]) +
        scale_linetype_manual(values = c("solid", "11")) +
        guides(colour = "none") +
        scale_y_continuous("weight (Kg)", breaks = seq(0, 24, 2),
          limits = c(0, 20), expand = c(0, 0)) +
        scale_x_continuous("age (months)", breaks = seq(0, 60, 5), 
          expand = c(0, 0))
      assign(paste("plot", i, sep = "_"), plot)
    }

    # Combined plot
    ggarrange(
      plot_reference + theme(plot.margin = margin(1,0.1,0.1,0.1, "cm")), 
      plot_actual + theme(plot.margin = margin(1,0.1,0.1,0.1, "cm")), 
      labels = c("recommended intake", "actual intake"), ncol = 1, nrow = 2, 
      align = "hv", hjust = 0, font.label = list(size = 11)
    )
    ggsave(paste0(dir_path, "out/12_w_wo_act_intake_combi.png"),
      dpi = "print", units = "cm", height = 15, width = 20)    
    
    # Save actual and reference intakes
    write.csv(df_ha, paste0(dir_path, "out/12_df_ha.csv"), row.names = F)
}
 

if (! file.exists(paste0(dir_path, "out/12_cal_ml.csv"))) {
#...............................................................................
### Calibrating Hall et al. model to reproduce individual growth variability
      # by grid search combined with least-squares method
#...............................................................................
    
  #...................................
  ## Prepare inputs
    
    # Read reference and actual intake
    df_ha <- read.csv(paste0(dir_path,"out/12_df_ha.csv"))

    # Read growth centiles from GAMLSS model
    out_centiles <- readRDS(paste0(dir_path, "out/11_centiles.rds"))
    out_centiles <- out_centiles[order(out_centiles$measure, out_centiles$sex,
      out_centiles$age), ]
    
    # Prepare data
    df_ls <- reshape(subset(out_centiles, measure == "weight"),
      direction = "long", idvar = c("sex", "age"), timevar = "centile_weight",
      varying = grep("^q[[:digit:]]", colnames(out_centiles), value = T),
      times = grep("^q[[:digit:]]", colnames(out_centiles), value = T),
      v.names = "weight", drop = "measure")
    df_ls$centile_weight <- factor(df_ls$centile, levels = paste0("c", 0:100))
    df_ls <- df_ls[order(df_ls$centile_weight, df_ls$sex, df_ls$age),
      c("centile_weight", "age", "sex", "weight")]
    df_ls$formula_fed <- F
    df_ls$pred <- NA
    df_ls$dia <- F
    df_ls$ari <- F
    
    # Parameter names
    par_names <- c("par0", "par1")
    
    # Blank set of individual calibration parameters
    par_blank <- c(
      "par0" = 0, # minimum adjustment over lifespan
      "par1" = 0 # sigma (standard deviation of log-normal function)
    )
    
    # Prior distribution ranges
    priors_min <- c(0, 0)
    priors_max <- c(0.15, 15)
    
    # Create parameter grid
    pars <- expand.grid(
      seq(priors_min[1], priors_max[1], 
        abs(priors_max[1] - priors_min[1]) / 50),
      seq(priors_min[2], priors_max[2], 
        abs(priors_max[2] - priors_min[2]) / 50)
    )
    colnames(pars) <- par_names
    
    # Initialise output
    out_ls <- expand.grid(sex=c("female","male"),
      centile_weight = unique(df_ls$centile_weight))
    out_ls <- merge(out_ls, pars)
    out_ls$rss <- NA
    out_ls <- out_ls[order(out_ls$sex, out_ls$centile_weight), ]
  
  #...................................
  ## Parallel implementation of loop through all sex-centile-grid combinations
  if (! file.exists(paste0(dir_path, "out/12_ls.rds"))) {
      
    # Prepare inputs
    out_ls$sex <- as.character(out_ls$sex)
    out_ls$centile_weight <- as.character(out_ls$centile_weight)
    df_ls$centile_weight <- as.character(df_ls$centile_weight)
    sets <- split(out_ls, seq(nrow(out_ls))) 

    # Parallel function
    f_parallel <- function(set_i, df_ls_f = df_ls, f_hall_pred_f = f_hall_pred, 
      par_names_f = par_names){
      
      # unlist set values
      set_i <- unlist(set_i)
      
      # subset data
      df_i <- subset(df_ls_f, sex == set_i["sex"] & 
        centile_weight == set_i["centile_weight"])
      
      # identify parameter values
      par_i <- as.numeric(set_i[par_names_f])
      names(par_i) <- par_names_f

      # predict weight for each age
      pred <- f_hall_pred_f(data_f = df_i, par_f = par_i, intake_f = "actual")
      df_i$pred <- pred$weight
      
      # compute and store residual sum of squares
      set_i["rss"] <- sum((df_i$weight - df_i$pred)^2)
      return(set_i)
    }
  
    # Prepare parallel cores
    cl <- parallel::makeCluster(detectCores())
    x <- c("df_ha", "df_ls", "f_hall_kids", "f_bw_fm", "f_hall_pred", 
      "par_names", "kcal_ratio_dia", "kcal_ratio_ari")
    clusterExport(cl, x, envir = environment())
    clusterEvalQ(cl, library("gamlss"))    

    # Run and collect output
    system.time(out <- parallel::parLapply(cl, sets, f_parallel))
      # 100 sets ~ 30 sec
    stopCluster(cl)
    x <- do.call(rbind, out)
    saveRDS(x, paste0(dir_path, "out/12_ls.rds"))
  }
  
    
  if (file.exists(paste0(dir_path, "out/12_ls.rds"))) {
    x <- readRDS(paste0(dir_path, "out/12_ls.rds"))
  }  
  
    
  #...................................
  ## Compute maximum likelihood parameter estimates
    
    # Find maximum likelihood sets of parameter values    
    x$centile_n <- as.numeric(gsub("q", "", x$centile_weight))
    par_names <- c("par0", "par1")
    out_ml <- by(x, x[, c("sex", "centile_n")], FUN = function(x) {
      x[which.min(x$rss), c("sex", "centile_n", "centile_weight", par_names)]})
    out_ml <- do.call(rbind, out_ml)
    
    # Save
    write.csv(out_ml, paste0(dir_path, "out/12_cal_ml.csv"), row.names = F)
      

  #...................................      
  ## Prepare calibration parameter values for each age, sex and weight centile
    
    # Initialise output
    cal_par <- expand.grid(sex = c("female", "male"), 
      age = sort(unique(centiles$age)), 
      centile_weight = unique(out_ml$centile_weight))
    cal_par <- cal_par[order(cal_par$sex, cal_par$centile_weight, cal_par$age),]
    cal_par$sex <- as.character(cal_par$sex)
    cal_par$centile_weight <- as.character(cal_par$centile_weight)
    cal_par[, c("par0", "par1")] <- NA
    cal_par <- cal_par[order(cal_par$centile_weight, cal_par$sex, cal_par$age),]
        
    # Prepare parameter values for each age-sex-centile combination
    for (i in 1:nrow(out_ml)) {
      
      # identify combination
      x <- cal_par[which(cal_par$sex == out_ml[i, "sex"] &
        cal_par$centile_weight == out_ml[i, "centile_weight"]), ]
      
      # identify parameter values
      par0 <- as.numeric(out_ml[i, "par0"])
      par1 <- as.numeric(out_ml[i, "par1"])
      
      # compute age-specific parameter values
      x <- dLOGNO(x$age / 100, par1)
      par_i <- suppressWarnings(data.frame(par0 = par0, 
        par1 = as.numeric(par0 * x/max(x))))
      
      # change sign for centiles >= 50% (= negative energy penalty)
      if (out_ml[i, "centile_n"] >= 50) {par_i <- -par_i}
      
      # add to output
      cal_par[which(cal_par$sex == out_ml[i, "sex"] &
        cal_par$centile_weight == out_ml[i, "centile_weight"]), 
        c("par0", "par1")] <- par_i
    }
    
    # Save
    saveRDS(cal_par, paste0(dir_path, "out/12_cal_par.rds"))
            
  #...................................
  ## Plot relationship between ML values and centile
  
    # Prepare for plotting
    df <- reshape(out_ml, direction = "long", varying = par_names,
      idvar = c("sex", "centile_weight", "centile_n"), timevar = "par",
      times = par_names, v.names = "value")
    df[which(df$centile_n >= 50 & df$par == "par0"), "value"] <- 
      - df[which(df$centile_n >= 50 & df$par == "par0"), "value"]
    panel_labels <- c("par0" = paste0("\u03B8", "_0"), 
      "par1" = paste0("\u03B8", "_1"))
      
    # Plot
    ggplot(df, aes(x = centile_n, y = value, colour = sex)) +
      geom_point() +
      geom_smooth() +
      facet_wrap(par ~., ncol = 3, scales = "free_y", 
        labeller =  labeller(par = panel_labels)) +
      theme_bw() +
      scale_colour_manual(values = palette_gen[c(9,3)]) +
      theme(legend.position = "top") +
      scale_x_continuous("growth percentile (weight)") +
      scale_y_continuous("best-fitting value")
    ggsave(paste0(dir_path, "out/12_cal_ml_values.png"),
      dpi = "print", units = "cm", height = 15, width = 20)    
    

  #...................................
  ## Predict and plot uncalibrated model, for certain centiles
  
    # Prepare output
    df_ls <- df_ls[order(df_ls$sex, df_ls$centile_weight, df_ls$age), ]
    df <- subset(df_ls, centile_weight %in% paste0("c", c(5, 50, 95)))
    df$unc <- NA
  
    # For each sex...
    for (i in c("female", "male")) {
      
      # for each centile...
      for (j in unique(df$centile_weight)) {
        
        # progress
        print(paste0("sex: ", i, ", percentile: ", j))
        
        # subset data
        df_i <- subset(df, sex == i & centile_weight == j)
        
        # predict weight for each age (uncalibrated)
        x <- f_hall_pred(data_f = df_i,par_f = par_blank,intake_f = "reference")
        df[which(df$sex == i & df$centile_weight == j), "unc"] <- x$weight
      }
    }
    
    # Prepare data for plotting
    df <- reshape(df, direction = "long", varying = c("weight", "unc"),
      idvar = c("sex", "centile_weight", "age"), timevar = "estimate", 
      times = c("fitted growth curve", "model (uncalibrated)"), 
      v.names = "value", drop = c("formula_fed", "pred"))
    df <- subset(df, sex == "female")
    df$centile_weight <- factor(df$centile_weight, 
      levels = paste0("c", c(5, 50, 95)), labels = paste0(c(5, 50, 95), "%"))

    # Plot
    ggplot(df, aes(x = age/30.44, y = value, colour = centile_weight, 
      linetype = estimate, linewidth = estimate)) +
      geom_line(alpha = 0.5) +
      scale_linetype_manual(values = c("solid", "31")) +
      scale_linewidth_manual(values = c(1.15, 0.75)) +
      theme_bw() +
      scale_colour_manual(values = palette_gen[c(2,8,14)]) +
      scale_y_continuous("weight (Kg)", breaks = seq(0, 24, 2),
        limits = c(0, 24), expand = c(0, 0)) +
      scale_x_continuous("age (months)", breaks = seq(0, 60, 12),
        limits = c(0, 60), expand = c(0, 0)) +
      theme(legend.position = "top")
    ggsave(paste0(dir_path, "out/12_curves_uncal_examples.png"),
      dpi = "print", units = "cm", height = 15, width = 20)    
    
  #...................................
  ## Predict and plot calibrated and uncalibrated model, by centile and sex
  
    # Prepare output
    df_ls <- df_ls[order(df_ls$sex, df_ls$centile_weight, df_ls$age), ]
    df_ls[, c("cal", "unc")] <- NA
  
    # For each sex...
    for (i in c("female", "male")) {
      
      # for each centile...
      for (j in unique(df_ls$centile_weight)) {
        
        # progress
        print(paste0("sex: ", i, ", percentile: ", j))
        
        # subset data
        df_i <- subset(df_ls, sex == i & centile_weight == j)
        
        # identify best-fit parameter values
        par_i <- out_ml[which(out_ml$sex == i & 
          out_ml$centile_weight == j), par_names]
        names(par_i) <- par_names
        
        # predict weight for each age (calibrated)
        x <- f_hall_pred(data_f = df_i, par_f = par_i, intake_f = "actual")
        df_ls[which(df_ls$sex==i & df_ls$centile_weight == j),"cal"] <- x$weight
        
        # predict weight for each age (uncalibrated)
        x <- f_hall_pred(data_f = df_i, par_f = par_blank, intake_f = "actual")
        df_ls[which(df_ls$sex==i & df_ls$centile_weight == j),"unc"] <- x$weight
      }
    }
  
    # Prepare data for plotting
    df <- reshape(df_ls, direction = "long", varying = c("weight", "cal","unc"),
      idvar = c("sex", "centile_weight", "age"), timevar = "estimate", 
      times = c("fitted growth curve", "model (calibrated)", 
        "model (uncalibrated)"), v.names = "value", 
      drop = c("formula_fed", "pred"))
    df$centile_weight <- factor(df$centile_weight, levels = paste0("c", 0:100), 
      labels = paste0(0:100, "%"))
    df <- subset(df, centile_weight != "100%") # to avoid extra row in plot
    
    # Plot
    ggplot(df, aes(x = age/30.44, y = value, colour = sex, 
      linetype = estimate)) +
      geom_line(linewidth = 0.5, alpha = 0.5) +
      scale_linetype_manual(values = c("solid", "11", "31")) +
      facet_wrap(centile_weight ~., ncol = 10, scales = "free_y") +
      theme_bw() +
      scale_colour_manual(values = palette_gen[c(9,3)]) +
      scale_y_continuous("weight (Kg)", breaks = seq(0, 24, 2),
        limits = c(0, 24), expand = c(0, 0)) +
      scale_x_continuous("age (months)", breaks = seq(0, 60, 12),
        limits = c(0, 60), expand = c(0, 0)) +
      theme(legend.position = "top")
    ggsave(paste0(dir_path, "out/12_curves_cal_vs_uncal.png"),
      dpi = "print", units = "cm", height = 50, width = 30)    
  
    # Plot with finer detail in early months of life
    ggplot(df, aes(x = age/30.44, y = value, colour = sex, 
      linetype = estimate)) +
      geom_line(linewidth = 0.5, alpha = 0.5) +
      scale_linetype_manual(values = c("solid", "11", "31")) +
      facet_wrap(centile_weight ~., ncol = 10, scales = "free_y") +
      theme_bw() +
      scale_colour_manual(values = palette_gen[c(9,3)]) +
      scale_y_continuous("weight (Kg)", breaks = seq(0, 24, 2),
        limits = c(0, 14), expand = c(0, 0)) +
      scale_x_continuous("age (months)", breaks = seq(0, 12, 2),
        limits = c(0, 12), expand = c(0, 0)) +
      theme(legend.position = "top")
    ggsave(paste0(dir_path, "out/12_curves_cal_vs_uncal_12mths.png"), 
      dpi = "print", units = "cm", height = 50, width = 30)    

    
    # Plot only a selection of centiles
    df <- subset(df, centile_weight %in% paste0(c(0,5,20,50,80,95,99), "%"))
    ggplot(df, aes(x = age/30.44, y = value, colour = sex, 
      linetype = estimate)) +
      geom_line(linewidth = 0.5, alpha = 0.5) +
      scale_linetype_manual(values = c("solid", "11", "31")) +
      facet_wrap(centile_weight ~., ncol = 4, scales = "free_y") +
      theme_bw() +
      scale_colour_manual(values = palette_gen[c(9,3)]) +
      scale_y_continuous("weight (Kg)", breaks = seq(0, 24, 2),
        limits = c(0, 24), expand = c(0, 0)) +
      scale_x_continuous("age (months)", breaks = seq(0, 60, 12),
        limits = c(0, 60), expand = c(0, 0)) +
      theme(legend.position = "top")
    ggsave(paste0(dir_path, "out/12_curves_cal_vs_uncal_subset.png"),
      dpi = "print", units = "cm", height = 15, width = 25)
}
        
#...............................................................................  
### ENDS
#...............................................................................

