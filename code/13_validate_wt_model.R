#...............................................................................
### + GENERIC MODEL OF ACUTE MALNUTRITION AS A FUNCTION OF CALORIC INTAKE ++ ###
#...............................................................................

#...............................................................................
## ------------ R SCRIPT TO VALIDATE WEIGHT MODEL AMONG CHILDREN  ----------- ##
#...............................................................................



#...............................................................................  
### Validating model predictions against WHO growth standards (weight-for-age)
#...............................................................................

  #...................................      
  ## Initialise test dataset: female and male, each starting out at 0SD weight

    # Initialise dataset
    test <- expand.grid(age = 0:as.integer(365.35*5), sex = c("female", "male"))
    test$formula_fed <- F
    test$weight_0sd <- NA
    
    # Set starting weights (birth) based on WHO standards
    for (i in c("female", "male")) {
      test[which(test$age == 0 & test$sex == i),"weight_0sd"] <- 
        who_waz_standards[which(who_waz_standards$age == 0 & 
          who_waz_standards$sex == i), "weight_median"]
    }
    
    # Sort
    test <- test[order(test$sex, test$age), ]
    
    # Assumed proportion of fat mass at birth
    x <- fomon_weight_composition[which(fomon_weight_composition$sex == "female" 
      & fomon_weight_composition$age == 0), ]
    prop_fm <- x$fm / (x$fm + x$ffm)
    x <- fomon_weight_composition[which(fomon_weight_composition$sex == "male" 
      & fomon_weight_composition$age == 0), ]
    prop_fm <- c(prop_fm, x$fm / (x$fm + x$ffm))
    names(prop_fm) <- c("female", "male")


  #...................................      
  ## Predict weight based on Hall model, by day
    
    # Read various outputs from previous scripts
    centiles <- readRDS(paste0(dir_path, "out/11_centiles.rds"))
    corr_mat <- readRDS(paste0(dir_path, "out/11_ht_wt_correlation.rds"))
    df_ha <- read.csv(paste0(dir_path, "out/12_df_ha.csv"))
    cal_ml <- read.csv(paste0(dir_path, "out/12_cal_ml.csv"))
    cal_par <- readRDS(paste0(dir_path, "out/12_cal_par.rds"))
    
    # For each sex...
    for (i in c("female", "male")) {
      
      # identify data      
      df_i <- test[which(test$sex == i), ]
      df_i[, c("weight", "fm", "ffm")] <- NA
      x <- which(df_i$age == 0)
      df_i[x, "weight"] <- df_i[x, "weight_0sd"] 
      df_i[x, "fm"] <- df_i[x, "weight"] * prop_fm[[i]]
      df_i[x, "ffm"] <- df_i[x, "weight"] - df_i[x, "fm"]
      df_i$centile_weight <- "c50"
      df_i$dia <- F
      df_i$ari <- F
      
      # compute weight by age
      x <- f_hall_pred(data_f = df_i, par_f = cal_ml[which(cal_ml$sex == i & 
        cal_ml$centile_n == 50), grep("par", colnames(cal_ml))], 
        intake_f = "reference")
      test[which(test$sex == i), "pred"] <- x$weight      
    }
  
  #...................................      
  ## Plot predictions as WAZ and as deviation from WHO standard WAZ
    
    # Compute WAZ scores
      # dataset to compute WAZ
      df <- test[, c("age", "sex", "pred")]
      df$weight <- df$pred
      df$sex <- ifelse(df$sex == "female", 2, 1)

      # compute WAZ
      x <- with(df, anthro_zscores(sex = sex, age = age, weight = weight))
      test$waz_0sd <- x$zwei
    
    # Add WHO reference weights
    x <- who_waz_standards[, c("age", "sex", "SD0", "SD2neg", "SD2")]
    test <- merge(test, x, by = c("age", "sex"), all.x = T)

    # Plot predictions as weight, compared to WHO reference weight
    ggplot(test, aes(x = age / 30.44, y = pred, colour = sex)) +
      geom_line(linewidth = 1, alpha = 0.5, linetype = "22") +
      geom_line(aes(y = SD0), col = palette_gen[15], linewidth = 1, 
        alpha = 0.5) +
      geom_ribbon(aes(ymin = SD2neg, ymax = SD2), alpha = 0.2,
        fill = palette_gen[15], colour = NA) +
      scale_x_continuous("age (months)", breaks = seq(0, 60, 5)) +
      scale_y_continuous("weight (Kg)") +
      theme_bw() +
      facet_wrap(sex ~., ncol = 2) +
      scale_colour_manual(values = palette_gen[c(7,2)]) +
      theme(legend.position = "none")
    ggsave(paste0(dir_path, "out/13_wt_model_vs_wt_who.png", sep=""),
      dpi = "print", units = "cm", height = 10, width = 22)  
    
    # Plot predictions as deviation from WHO WAZ
    ggplot(test, aes(x = age / 30.44, y = waz_0sd, colour = sex)) +
      geom_line(linewidth = 1, alpha = 0.5, linetype = "22") +
      scale_x_continuous("age (months)", breaks = seq(0, 60, 5)) +
      scale_y_continuous("difference in weight-for-age Z-scores") +
      geom_hline(aes(yintercept = 0), colour = palette_gen[15], alpha = 0.5) +
      theme_bw() +
      facet_wrap(sex ~., ncol = 2) +
      scale_colour_manual(values = palette_gen[c(7,2)]) +
      theme(legend.position = "none")
    ggsave(paste0(dir_path, "out/13_waz_model_vs_waz_who.png", sep=""),
      dpi = "print", units = "cm", height = 10, width = 22)  
 
    
#...............................................................................  
### Comparing calibrated model predictions to observations used to construct
        # pre-crisis growth curves
#...............................................................................

  #...................................      
  ## Read required objects
    
    # Various outputs from previous scripts
    centiles <- readRDS(paste0(dir_path, "out/11_centiles.rds"))
    corr_mat <- readRDS(paste0(dir_path, "out/11_ht_wt_correlation.rds"))
    df_ha <- read.csv(paste0(dir_path, "out/12_df_ha.csv"))
    cal_ml <- read.csv(paste0(dir_path, "out/12_cal_ml.csv"))
    cal_par <- readRDS(paste0(dir_path, "out/12_cal_par.rds"))
    indices <- c("zwei", "zlen", "zwfl")
    sources <- c("co", "obs")

    # Individual growth monitoring observations (one observation per child)
    if (! file.exists(paste0(dir_path, "out/13_ind_growth_curves_unique.rds"))){
      
      # read all observations    
      ind <- readRDS(paste0(dir_path, "out/11_ind_growth_curves.rds"))
      
      # choose a single observation per child (to avoid clustering)
      x <- unique(ind[, c("id", "n_obs")])
      x$which_n <- apply(x, 1, function(x) {sample.int(x["n_obs"], 1)})
      x$cum_row <- cumsum(x$n_obs)
      x$which_row <- x$cum_row + x$n_obs
      ind_unique <- ind[x$which_row, ]
      
      # add anthropometric indices
      ind_unique$sex_anthro <- ifelse(ind_unique$sex == "female", "f", "m")
      x <- with(ind_unique, anthro_zscores(sex=sex_anthro, age = age, 
        lenhei= height, weight = weight))
      ind_unique <- cbind(ind_unique, x)  

      # save
      saveRDS(ind_unique,paste0(dir_path,"out/13_ind_growth_curves_unique.rds"))
      rm(ind)
    }
    if (file.exists(paste0(dir_path, "out/13_ind_growth_curves_unique.rds"))) {
      ind_unique <- readRDS(
        paste0(dir_path,"out/13_ind_growth_curves_unique.rds"))
    }
    obs <- ind_unique
    
  #...................................      
  ## Initialising a cohort of children and running it the crisis' start
  if (! file.exists(paste0(dir_path, "out/13_test_cohort.rds"))) {
    
    # Objects for cohort initiation function (f_co)
    x <- c("areas", "cal_ml", "centiles", "corr_mat", "df_ha", "f_bw_fm",
      "f_hall_kids", "f_hall_pred", "n_kids")
    obj_co <- vector("list", length = length(x))
    names(obj_co) <- x
    for (i in x) {obj_co[[i]] <- get(i)}
    
    # Initiate cohort of 10000 children and save it
    obj_co$n_kids <- 10000
    co <- f_co()
    saveRDS(co, paste0(dir_path, "out/13_test_cohort.rds"))
  }
  if (file.exists(paste0(dir_path, "out/13_test_cohort.rds"))) 
    {co <- readRDS(paste0(dir_path, "out/13_test_cohort.rds"))}
      
      
  #...................................      
  ## Compare modelled and observed anthropometry statistics

    # Initiate table
    tab <- expand.grid(source = sources, index = indices)
    tab[, c("mean", "sd", "below2sd", "below3sd", "missing")] <- NA
    tab <- tab[order(tab$source, tab$index), ]
 
    # Compute values
    for (i in sources) {
      df <- get(i)
      for (j in indices) {
        x <- which(tab$source == i & tab$index == j)
        tab[x, "missing"] <- length(which(is.na(df[, j]))) / length(df[, j])
        df_c <- df[, j]
        df_c <- df_c[! is.na(df_c)]
        tab[x, "mean"] <- mean(df_c)
        tab[x, "sd"] <- sd(df_c)
        tab[x, "below2sd"] <- length(which(df_c < -2)) / length(df_c)
        tab[x, "below3sd"] <- length(which(df_c < -3)) / length(df_c)
      }
    }
    
    # Format and save
    for (i in c("mean", "sd")) {tab[, i] <- format(round(tab[, i], 2),nsmall=2)}
    for (i in c("below2sd", "below3sd", "missing")) {
      tab[, i] <- percent(tab[, i ], 0.1)}
    tab$source <- rep(c("model", "observed"), each = length(unique(tab$index)))   
    tab$index <- rep(c("WAZ", "HAZ", "WHZ"), length(unique(tab$source)))
    write.csv(tab, paste0(dir_path, "out/13_model_vs_data_pre_crisis.csv"),
      row.names = F)

  #...................................      
  ## Compare modelled and observed anthropometry distributions

    # Plot the distributions of each index
    for (i in sources) {
      df <- get(i)
      for (j in indices) {
        # x axis title
        if (j == "zwei") {xlab <- "weight for age Z-score"}
        if (j == "zlen") {xlab <- "height/length for age Z-score"}
        if (j == "zwfl") {xlab <- "weight for height/length Z-score"}
        
        # plot
        df$variable <- df[, j]
        pl <- ggplot(df, aes(x = variable, colour = sex)) +
          geom_density(linewidth = 1) +
          scale_colour_manual(values = palette_gen[c(9,3)]) +
          theme_bw() +
          theme(legend.position = "top") +
          scale_x_continuous(xlab, limits = c(-5, 5), expand = c(0,0)) +
          scale_y_continuous(limits = c(0, 0.55), 
            expand = c(0,0))

        # assign plot name
        assign(paste("pl", i, j, sep = "_"), pl)
      }
    }

    # Combined plots
    for (i in indices) {
      x <- mget(paste("pl", sources, i, sep = "_"))
      pl <- ggarrange(plotlist = x, nrow=2, labels = c("modelled", "observed"),
        hjust = 0, font.label = list(size = 11), common.legend = T)
      ggsave(paste0(dir_path, "out/13_model_vs_data_", i, "_combi.png"),
        units = "cm", dpi = "print", height = 25, width = 15)
    }
    x <- mget(paste("pl", paste(rep(sources, each = length(indices)), indices,
      sep = "_"), sep="_"))    
    pl <- ggarrange(plotlist = x, nrow = length(sources), ncol = length(indices),
      labels = c("modelled", rep("", length(indices) - 1), 
        "observed", rep("", length(indices) - 1)), hjust = 0, 
      font.label = list(size = 11))
    ggsave(paste0(dir_path, "out/13_model_vs_data_all_combi.png"),
      units = "cm", dpi = "print", height = 25, width = 45)


  #...................................      
  ## Compare mean modelled indices by age with Albelbeisi et al. 
      # (https://doi.org/10.26719/2018.24.3.302)
    
    # Paper's estimates
    paper <- data.frame(
      age_m = c(0,6,9,12,15,18,24),
      waz = c(-0.07,-0.35,-0.23,0.05,-0.09,-0.03,-0.11),
      haz = c(0.06,0.17,-0.06,-0.18,-0.06,-0.25,-0.85),
      whz = c(-0.17,-0.40,-0.09,0.18,0.17,0.18,0.34)
    )

    # Model's estimates
    model <- paper
    for (i in model$age_m) {
      model[which(model$age_m == i), "waz"] <- 
        mean(co[which(floor(co$age / 30.44) == i), "zwei"], na.rm = T)
      model[which(model$age_m == i), "haz"] <- 
        mean(co[which(floor(co$age / 30.44) == i), "zlen"], na.rm = T)
      model[which(model$age_m == i), "whz"] <- 
        mean(co[which(floor(co$age / 30.44) == i), "zwfl"], na.rm = T)
    }  
    
    # Plot both estimates
    df <- rbind(model, paper)
    df$source <- rep(c("model", "Albelbeisi et al. (2018)"), each = nrow(df)/2)
    df <- reshape(df, direction = "long", varying = c("waz", "haz", "whz"),
      idvar = c("age_m", "source"), timevar = "index", 
      times = c("weight for age Z-score", "height/length for age Z-score", 
        "weight for height/length Z-score"), v.names = "mean")
    ggplot(df, aes(x = age_m, y = mean, group = source, colour = source,
      fill = source, linetype = source)) +
      geom_point() +
      geom_line() +
      theme_bw() +
      theme(legend.position = "top") +
      facet_wrap(index ~., ncol = 3) +
      scale_x_continuous("age (months", expand = c(0,0)) +
      scale_y_continuous("mean Z-score") +
      scale_colour_manual("source", values = palette_gen[c(4,12)]) +
      scale_fill_manual("source", values = palette_gen[c(4,12)]) +
      scale_linetype_manual("source", values = c("solid", "11"))      
    ggsave(paste0(dir_path, "out/13_model_vs_albelbeisi.png"),
      units = "cm", dpi = "print", height = 15, width = 30)
  
            
#...............................................................................  
### ENDS
#...............................................................................
     