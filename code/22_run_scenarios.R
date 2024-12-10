#...............................................................................
### + GENERIC MODEL OF ACUTE MALNUTRITION AS A FUNCTION OF CALORIC INTAKE ++ ###
#...............................................................................

#...............................................................................
## --------- R SCRIPT TO RUN RETROSPECTIVE AND SCENARIO SIMULATION  --------- ##
#...............................................................................

  
#...............................................................................  
### Preparing for the simulation
#...............................................................................

  #...................................      
  ## Read necessary inputs and prepare list of objects that go into functions
    
    # Various outputs from previous scripts
    centiles <- readRDS(paste0(dir_path, "out/11_centiles.rds"))
    corr_mat <- readRDS(paste0(dir_path, "out/11_ht_wt_correlation.rds"))
    df_ha <- read.csv(paste0(dir_path, "out/12_df_ha.csv"))
    cal_ml <- read.csv(paste0(dir_path, "out/12_cal_ml.csv"))
    cal_par <- readRDS(paste0(dir_path, "out/12_cal_par.rds"))
    equi <- readRDS(paste0(dir_path, "out/14_disease_equi.rds"))
    sacrifice <- readRDS(paste0(dir_path, "out/15_sacrifice.rds"))
    if (retro == "yes") {rtl <- readRDS(paste0(dir_path, "out/21_rtl.rds"))}
    if (project == "yes") {scens <-readRDS(paste0(dir_path,"out/21_scens.rds"))}

    # Objects for cohort initiation function (f_co)
    x <- c("areas", "cal_ml", "centiles", "corr_mat", "df_ha", "f_bw_fm",
      "f_hall_kids", "f_hall_pred", "n_kids")
    obj_co <- vector("list", length = length(x))
    names(obj_co) <- x
    for (i in x) {obj_co[[i]] <- get(i)}
    
    # Objects for crisis period simulation function (f_sim)
    x <- c("areas", "b_corr", "cal_par", "centiles", "corr_mat", "df_ha", 
      "dur_scale_ari", "dur_scale_dia", "dur_shape_ari", "dur_shape_dia",
      "equi", "f_bw_fm", "f_hall_kids", "inc_scale_ari", 
      "inc_scale_dia", "intake", "kcal_ratio_ari", "kcal_ratio_dia", 
      "mam_tx_daily", "n_kids", "req_intake", "sacrifice", "sam_tx_daily")
    obj_sim <- vector("list", length = length(x))
    names(obj_sim) <- x
    for (i in x) {if (exists(i)) {obj_sim[[i]] <- get(i)} }

  #...................................      
  ## Initiate simulation output

    # Output of all retrospective analysis runs
    if (retro == "yes") {
      out_sim <- expand.grid(area = areas, 
        date = as.Date(date_crisis : (date_start - 1)), run = 1:n_runs)
      out_sim$scenario <- "retrospective"
      out_sim[, outcomes] <- NA
    }
    
    # Initialise scenario projection timeline by area and outputs of scenarios
    if (project == "yes") {
      
      # initialise timeline for each scenario
      stl <- expand.grid(date = as.Date(date_start : date_end), area = areas)
      stl[, c("scenario", "intake", 
        paste("tx_mam", c("coverage", "admissions"), sep = "_"),
        paste("tx_sam", c("coverage", "admissions"), sep = "_"),
        paste("ari", c("period_prev", "point_prev", "incidence"), sep = "_"),
        paste("dia", c("period_prev", "point_prev", "incidence"), sep = "_"),
        "prop_formula", "req_ch_protected", "req_ad_protected", outcomes)] <- NA     

      # add population data if they exist
      if (exists("df_po")) {
        stl <- merge(stl, df_po, by = c("date", "area"), all.x = T)
      }
      stl <- stl[order(stl$area, stl$date), ]
      
      # initialise  / add output of all runs of each scenario projection
      x <- expand.grid(scenario = names(scens), area = areas, 
        date = as.Date(date_start:date_end), run = 1:n_runs)
      x[, outcomes] <- NA
      if (retro == "yes") {out_sim <- rbind(out_sim, x[, colnames(out_sim)])}
      if (retro == "no") {out_sim <- x}
    }



#...............................................................................  
### Running a simulation (many runs) in parallel
      # retrospective estimates + 1 projection estimate per scenario)
#...............................................................................

# clusters <- makeCluster(detectCores() - 1)
# registerDoParallel(clusters)
# x <- c("obj_co", "obj_sim", "retro", "project", "n_runs", "out_sim", "f_co",
#   "f_sim", "df_ha")
# if (retro == "yes") {x <- c(x, "rtl")}
# if (project == "yes") {x <- c(x, "scens", "stl")}
# foreach (run_i = 1:n_runs, .export = x, .packages = c("gamlss")) %dopar% {

for (run_i in 1:n_runs) {    
  #...................................      
  ## Preparatory steps
  
    # Progress
    print(paste0("run ", run_i, " of ", n_runs))

    # Initialise a new cohort and run the cohort to the start of the crisis
    print(" initialising a new cohort")
    co_i <- f_co()  # 100 kids x area ~ 1 min
    
  #...................................      
  ## Implement retrospective estimation, if desired
  if (retro == "yes") {
    
    # Initialise fresh timeline
    rtl_i <- rtl
    
    # If absolute daily number of SAM/MAM treatment admissions is provided,
          # draw a Poisson daily random count
    if (! all(is.na(rtl_i$tx_mam_admissions))) {
      rtl_i$tx_mam_admissions <- rpois(1:nrow(rtl_i), rtl_i$tx_mam_admissions)
    }
    if (! all(is.na(rtl_i$tx_sam_admissions))) {
      rtl_i$tx_sam_admissions <- rpois(1:nrow(rtl_i), rtl_i$tx_sam_admissions)        
    }

    # Take the cohort through the crisis period
    print(" retrospective estimation")
    out_rtl <- f_sim(co_f = co_i, tl_f = rtl_i) # 100 kids x area ~ 3 min
    
    # Manage output
    rtl_i <- out_rtl$tl
    rtl_i$run <- run_i

    # Update output
    out_sim[which(out_sim$run==run_i & out_sim$scenario=="retrospective"), ] <- 
      rtl_i[order(rtl_i$area, rtl_i$date), colnames(out_sim)]
  }
    
  #...................................      
  ## Implement scenario projection, if desired    
  if (project == "yes") {
    
    # For each scenario...
    for (i in names(scens)) {    
      # progress
      print(paste0(" projection for scenario ", i))
      
      # cohort at the start of the projection period
      if (retro == "yes") {co_start <- out_rtl$co}
      if (retro == "no") {co_start <- co_i}
      
      # select scenario and initialise its timeline
      scen_i <- scens[[i]]
      stl_i <- stl

      # compute value of each factor over the scenario period
      for (j in 1:nrow(scen_i)) {
        # for each row of the scenario...
        x <- which(stl_i$area == scen_i[j, "area"] & 
          stl_i$date %in% as.Date(scen_i[j, "date_1"] : scen_i[j, "date_2"]))
        
        # if absolute value is provided...
        if (! is.na(scen_i[j, "abs_value"])) {
          stl_i[x, scen_i[j, "variable"]] <-  scen_i[j, "abs_value"]
        }
        # if relative value is provided (which implies retro estimate!)...
        if (! is.na(scen_i[j, "rel_risk"])) {
          stl_i[x, scen_i[j, "variable"]] <-  scen_i[j, "rel_risk"] *
            mean(rtl_i[which(rtl_i$area == scen_i[j, "area"] & rtl_i$date %in% 
             as.Date((max(rtl_i$date)-scen_i[j,"days_ref"]+1):max(rtl_i$date))), 
             scen_i[j, "variable"]])
        }
      }
      
      # make sure each value is bounded by 0 and 1 if it is a proportion
      x <- c("tx_mam_coverage", "tx_sam_coverage", "ari_point_prev",
        "dia_point_prev", "prop_formula", "req_ch_protected","req_ad_protected")
      for(j in x) {
        stl_i[, j] <- pmin(stl_i[, j], 1)
        stl_i[, j] <- pmax(stl_i[, j], 0)
      }

      # if absolute number of SAM/MAM treatment admissions are provided,
            # scale them to the cohort size, draw a Poisson daily random count
      if (! all(is.na(stl_i$tx_mam_admissions))) {
        stl_i$tx_mam_admissions<-stl_i$tx_mam_admissions*n_kids/stl_i$pop_0059
        stl_i$tx_mam_admissions <- rpois(1:nrow(stl_i), stl_i$tx_mam_admissions)
      }
      if (! all(is.na(stl_i$tx_sam_admissions))) {
        stl_i$tx_sam_admissions<-stl_i$tx_sam_admissions*n_kids/stl_i$pop_0059
        stl_i$tx_sam_admissions <- rpois(1:nrow(stl_i), stl_i$tx_sam_admissions)        
      }
      
      # establish equivalence of prevalence and incidence of ARI and diarrhoea
        # specify indicators
        indicators <- c("period_prev", "point_prev", "incidence")
        equi$incidence <- equi$inc_shape
      
      for (k in c("ari", "dia")) {
        
        # select equivalence table for this disease
        equi_i <- equi[which(equi$disease == k), ]
        
        # establish which indicator is non-missing during the timeline
        which_ind <- NA
        for (j in indicators) {
          if (! all(is.na(stl_i[, paste(k, j, sep = "_")]))) {which_ind <- j} }
        
        # fill in the others based on equivalence  
        for (j in indicators[indicators != which_ind]) {
          var <- paste(k, j, sep = "_")
          stl_i[, var] <- na.replace(approx(equi_i[, which_ind], equi_i[, j], 
            xout = stl_i[, paste(k, which_ind, sep = "_")])$y, 0)
        }
      }
      
      # simulate scenario
      out_stl <- f_sim(co_f = co_start, tl_f = stl_i)
      
      # manage output
      stl_i <- out_stl$tl
      stl_i$run <- run_i
      stl_i$scenario <- i

      # update output
      out_sim[which(out_sim$run == run_i & out_sim$scenario == i), ] <- 
        stl_i[order(stl_i$area, stl_i$date), colnames(out_sim)]
    }
  }
} # close runs loop    
# stopCluster(cl = clusters) # stop cluster

    
#...............................................................................  
### Visualising model output
#...............................................................................
    
  #...................................      
  ## Aggregate and compute uncertainty intervals

    # Compute SAM, MAM and GAM prevalence
    out_sim$pmam_0659 <- out_sim$mam_0659 / out_sim$n_0659
    out_sim$psam_0659 <- out_sim$sam_0659 / out_sim$n_0659
    out_sim$pgam_0659 <- out_sim$psam_0659 + out_sim$pmam_0659
    out_sim$pmam_0623 <- out_sim$mam_0623 / out_sim$n_0623
    out_sim$psam_0623 <- out_sim$sam_0623 / out_sim$n_0623
    out_sim$pgam_0623 <- out_sim$psam_0623 + out_sim$pmam_0623
    outcomes_new <- c(outcomes, "pmam_0659", "psam_0659", "pgam_0659",
      "pmam_0623", "psam_0623", "pgam_0623")
    
    # Compute mean, medians and 80%, 95% CIs of run outputs (and reshape long)
    out_agg <- aggregate(out_sim[, outcomes_new], by = out_sim[, c("area", 
      "scenario", "date")], 
      function(x) {c(mean(x, na.rm = T), 
        quantile(x, c(0.50, 0.10, 0.90, 0.025, 0.975), na.rm = T))})
    x <- as.data.frame(do.call(rbind, out_agg[, outcomes_new]))
    x$outcome <- rep(outcomes_new, each = nrow(out_agg))
    x[,c("area", "scenario", "date")] <- out_agg[,c("area", "scenario", "date")]
    out_agg <- x
    colnames(out_agg) <- c("mean", "median", "lci80", "uci80", "lci95", "uci95",
      "outcome", "area", "scenario", "date")
    x <- c()
    if (retro == "yes") {x <- c("retrospective")}
    if (project == "yes") {x <- c(x, names(scens))}
    out_agg$scenario <- factor(out_agg$scenario, levels = x)
    out_agg$outcome <- factor(out_agg$outcome, levels = outcomes_new,
      labels = c("SAM cases", "MAM cases", "SAM incidence", "MAM incidence",
        "children 6-59mo", "SAM cases (6-59mo)", "MAM cases (6-59mo)", 
        "children 6-23mo","SAM cases (6-23mo)","MAM cases (6-23mo)", "mean WHZ", 
        "SAM on treatment", "MAM on treatment", 
        "SAM treatment coverage", "MAM treatment coverage",
        "diarrhoea cases", "ARI cases",
        "proportion formula-fed", "infants <6mo", "Kcal intake",
        "MAM prevalence (6-59mo)",
        "SAM prevalence (6-59mo)", "GAM prevalence (6-59mo)", 
        "MAM prevalence (6-23mo)", "SAM prevalence (6-23mo)", 
        "GAM prevalence (6-23mo)"))

    # Save simulation output
    saveRDS(out_sim, paste0(dir_path, "out/22_out_sim.rds"))
    saveRDS(out_agg, paste0(dir_path, "out/22_out_agg.rds"))
    
    
  #...................................      
  ## Produce generic visualisations
    
    # Graph trends in all the outcomes (point estimate only), by area
    pl_list <- vector("list", length(areas))
    for (i in 1:length(areas)) {
      pl_list[[i]] <- ggplot(subset(out_agg, area == areas[i]), aes(x = date, 
        y = median, colour = scenario, linetype = scenario)) +
        geom_step() +
        theme_bw() +
        scale_x_date("date", date_labels = "%b-%Y", breaks = "1 month", 
          expand = c(0,0)) +
        scale_y_continuous("value") +
        scale_linetype_manual("scenario", 
          values = c("solid","11","31","22","12","33")) +
        scale_colour_manual("scenario", values = palette_gen[c(1,5,9,15,12,3)])+
        facet_wrap(. ~ outcome, ncol = 5, scales = "free_y") +
        theme(legend.position = "top", axis.text.x = element_text(angle = 45,
          hjust = 1, vjust = 1), axis.title.x = element_blank(),
          plot.margin = margin(2,0,0,0, unit = "cm"))
    }
    ggarrange(plotlist = pl_list, ncol = 1, nrow = length(areas), 
      labels = areas, common.legend = T)
    ggsave(paste0(dir_path, "out/22_results_all_outcomes.png"), units = "cm",
      dpi = "print", width = 40, height = 30 * length(areas))
    
    # Graph trends in prevalence of malnutrition (with CIs), by area
    ncols <- length(unique(grep("SAM prevalence|GAM prevalence",
      out_agg$outcome, value = T)))
    df <- out_agg[grepl("SAM prevalence|GAM prevalence",out_agg$outcome),]
    ggplot(df, aes(x = date, 
      y = median, colour = scenario, fill = scenario, linetype = scenario)) +
      geom_step() +
      geom_stepribbon(aes(ymin = lci95, ymax = uci95), alpha = 0.1, colour=NA) +
      geom_stepribbon(aes(ymin = lci80, ymax = uci80), alpha = 0.2, colour=NA) +
      theme_bw() +
      scale_x_date("date", date_labels = "%b-%Y", breaks = "1 month", 
        expand = c(0,0)) +
      scale_y_continuous("prevalence", labels = percent, limits = c(0, NA)) +
      scale_linetype_manual("scenario", 
        values = c("solid","11","31","22","12","33")) +
      scale_colour_manual("scenario", values = palette_gen[c(1,5,9,15,12,3)]) +
      scale_fill_manual("scenario", values = palette_gen[c(1,5,9,15,12,3)]) +
      facet_wrap(area ~ outcome, scales = "free_y", ncol = ncols) +
      theme(legend.position = "top", axis.text.x = element_text(angle = 45,
        hjust = 1, vjust = 1), axis.title.x = element_blank())
    ggsave(paste0(dir_path, "out/22_results_nut_prevalence.png"), units="cm",
      dpi = "print", width = 10 * ncols, height = 10 * length(areas))
         
    # Graph trends in incidence of malnutrition (with CIs), by area
    ncols <- length(unique(grep("SAM incidence|MAM incidence",
      out_agg$outcome, value = T)))
    df <- out_agg[grepl("SAM incidence|MAM incidence",out_agg$outcome),]
    ggplot(df, aes(x = date, 
      y = median, colour = scenario, fill = scenario, linetype = scenario)) +
      geom_step() +
      geom_stepribbon(aes(ymin = lci95, ymax = uci95), alpha = 0.1, colour=NA) +
      geom_stepribbon(aes(ymin = lci80, ymax = uci80), alpha = 0.2, colour=NA) +
      theme_bw() +
      scale_x_date("date", date_labels = "%b-%Y", breaks = "1 month", 
        expand = c(0,0)) +
      scale_y_continuous("incidence", limits = c(0, NA)) +
      scale_linetype_manual("scenario", 
        values = c("solid","11","31","22","12","33")) +
      scale_colour_manual("scenario", values = palette_gen[c(1,5,9,15,12,3)]) +
      scale_fill_manual("scenario", values = palette_gen[c(1,5,9,15,12,3)]) +
      facet_wrap(area ~ outcome, scales = "free_y", ncol = ncols) +
      theme(legend.position = "top", axis.text.x = element_text(angle = 45,
        hjust = 1, vjust = 1), axis.title.x = element_blank())
    ggsave(paste0(dir_path, "out/22_results_nut_incidence.png"), units="cm",
      dpi = "print", width = 10 * ncols, height = 10 * length(areas))
    
    # Graph trends in weight-for-height (with CIs), by area
    df <- subset(out_agg, outcome == "mean WHZ")
    pl <- ggplot(df, aes(x = date, 
      y = median, colour = scenario, fill = scenario, linetype = scenario)) +
      geom_step() +
      geom_stepribbon(aes(ymin = lci95, ymax = uci95), alpha = 0.1, colour=NA) +
      geom_stepribbon(aes(ymin = lci80, ymax = uci80), alpha = 0.2, colour=NA) +
      theme_bw() +
      scale_x_date("date", date_labels = "%b-%Y", breaks = "1 month", 
        expand = c(0,0)) +
      scale_y_continuous("mean weight-for-height Z-score") +
      scale_linetype_manual("scenario", 
        values = c("solid","11","31","22","12","33")) +
      scale_colour_manual("scenario", values = palette_gen[c(1,5,9,15,12,3)]) +
      scale_fill_manual("scenario", values = palette_gen[c(1,5,9,15,12,3)]) +
      theme(legend.position = "top", axis.text.x = element_text(angle = 45,
        hjust = 1, vjust = 1), axis.title.x = element_blank())
    if (length(areas) > 1) {pl <- pl + facet_wrap(. ~ area, ncol = 2, 
      scales = "free_y")}
    ggsave(paste0(dir_path, "out/22_results_whz_only.png"), units = "cm",
      dpi = "print", width = 15 * length(areas), height = 15)

    
  #...................................      
  ## Produce visualisations for publication
    
    # Graph trends in select nutritional outcomes, by area (with CIs)
    x <- c("SAM prevalence (6-59mo)", "GAM prevalence (6-59mo)",
      "GAM prevalence (6-23mo)")
    ncols <- length(x)
    df <- out_agg[which(out_agg$outcome %in% x), ]
    pl_a <- ggplot(df, aes(x = date, 
      y = median, colour = scenario, fill = scenario, linetype = scenario)) +
      geom_step() +
      geom_stepribbon(aes(ymin = lci95, ymax = uci95), alpha = 0.1, colour=NA) +
      geom_stepribbon(aes(ymin = lci80, ymax = uci80), alpha = 0.2, colour=NA) +
      theme_bw() +
      scale_x_date("date", date_labels = "%b-%Y", breaks = "1 month", 
        expand = c(0,0)) +
      scale_y_continuous("prevalence", labels = percent, limits = c(0, NA),
        breaks = seq(0, 1, 0.2)) +
      scale_linetype_manual("scenario", 
        values = c("solid","11","31","22","12","33")) +
      scale_colour_manual("scenario", values = palette_gen[c(1,5,9,15,12,3)]) +
      scale_fill_manual("scenario", values = palette_gen[c(1,5,9,15,12,3)]) +
      facet_grid(outcome ~ area, scales = "free_y") +
      theme(legend.position = "top", axis.text.x = element_blank(), 
        axis.title.x = element_blank(), axis.ticks.x = element_blank())
    df <- out_agg[which(out_agg$outcome == "mean WHZ"), ]
    pl_b <- ggplot(df, aes(x = date, 
      y = median, colour = scenario, fill = scenario, linetype = scenario)) +
      geom_step() +
      geom_stepribbon(aes(ymin = lci95, ymax = uci95), alpha = 0.1, colour=NA) +
      geom_stepribbon(aes(ymin = lci80, ymax = uci80), alpha = 0.2, colour=NA) +
      theme_bw() +
      scale_x_date("date", date_labels = "%b-%Y", breaks = "1 month", 
        expand = c(0,0)) +
      scale_y_continuous("mean WHZ", limits = c(NA, NA)) +
      scale_linetype_manual("scenario", 
        values = c("solid","11","31","22","12","33")) +
      scale_colour_manual("scenario", values = palette_gen[c(1,5,9,15,12,3)]) +
      scale_fill_manual("scenario", values = palette_gen[c(1,5,9,15,12,3)]) +
      facet_grid(outcome ~ area, scales = "free_y") +
      theme(legend.position = "top", axis.text.x = element_text(angle = 30,
        hjust = 0.7, vjust = 0.7), axis.title.x = element_blank(),
        strip.background.x = element_blank(), strip.text.x = element_blank())
    ggarrange(pl_a, pl_b, labels = NA, common.legend = T, heights = c(3,1),
      ncol = 1, nrow = 2, align = "v")
    ggsave(paste0(dir_path, "out/22_results_nut_select_wide.png"), units="cm",
      dpi = "print", width = 8 * (ncols + 1), height = 10 * length(areas), 
      bg="white")
    ggsave(paste0(dir_path, "out/22_results_nut_select_long.png"), units="cm",
      dpi = "print", height = 6 * (ncols + 1), width = 10 * length(areas),
       bg="white")

    # Graph trends in select nutritional outcomes, by area (simpler without CIs)
    x <- c("SAM prevalence (6-59mo)", "GAM prevalence (6-59mo)")
    ncols <- length(x)
    df <- out_agg[which(out_agg$outcome %in% x), ]
    pl_a <- ggplot(df, aes(x = date, 
      y = median, colour = scenario, fill = scenario, linetype = scenario)) +
      geom_step(linewidth = 1) +
      # geom_stepribbon(aes(ymin = lci95, ymax = uci95), alpha = 0.1, colour=NA) +
      # geom_stepribbon(aes(ymin = lci80, ymax = uci80), alpha = 0.2, colour=NA) +
      theme_bw() +
      scale_x_date("date", date_labels = "%b-%Y", breaks = "1 month", 
        expand = c(0,0)) +
      scale_y_continuous("prevalence", labels = percent, limits = c(0, NA),
        breaks = seq(0, 1, 0.05), expand = expansion(mult=0, add=c(0,0.02))) +
      scale_linetype_manual("scenario", 
        values = c("solid","11","31","22","12","33")) +
      scale_colour_manual("scenario", values = palette_gen[c(1,5,9,15,12,3)]) +
      scale_fill_manual("scenario", values = palette_gen[c(1,5,9,15,12,3)]) +
      facet_grid(outcome ~ area, scales = "free_y") +
      theme(legend.position = "top", axis.text.x = element_blank(), 
        axis.title.x = element_blank(), axis.ticks.x = element_blank())
    df <- out_agg[which(out_agg$outcome == "mean WHZ"), ]
    pl_b <- ggplot(df, aes(x = date, 
      y = median, colour = scenario, fill = scenario, linetype = scenario)) +
      geom_step(linewidth = 1) +
      # geom_stepribbon(aes(ymin = lci95, ymax = uci95), alpha = 0.1, colour=NA) +
      # geom_stepribbon(aes(ymin = lci80, ymax = uci80), alpha = 0.2, colour=NA) +
      theme_bw() +
      scale_x_date("date", date_labels = "%b-%Y", breaks = "1 month", 
        expand = c(0,0)) +
      scale_y_continuous("mean WHZ", limits = c(NA, NA)) +
      scale_linetype_manual("scenario", 
        values = c("solid","11","31","22","12","33")) +
      scale_colour_manual("scenario", values = palette_gen[c(1,5,9,15,12,3)]) +
      scale_fill_manual("scenario", values = palette_gen[c(1,5,9,15,12,3)]) +
      facet_grid(outcome ~ area, scales = "free_y") +
      theme(legend.position = "top", axis.text.x = element_text(angle = 30,
        hjust = 0.7, vjust = 0.7), axis.title.x = element_blank(),
        strip.background.x = element_blank(), strip.text.x = element_blank())
    ggarrange(pl_a, pl_b, labels = NA, common.legend = T, heights = c(3,1),
      ncol = 1, nrow = 2, align = "v")
    ggsave(paste0(dir_path, "out/22_results_nut_select_wide_noci.png"), 
      units = "cm", dpi = "print", width = 10 * (ncols + 1), 
      height = 8 * length(areas), bg="white")
    ggsave(paste0(dir_path, "out/22_results_nut_select_long_noci.png"), 
      units = "cm", dpi = "print", width = 10 * length(areas),
       height = 8 * (ncols + 1), bg="white")

           
#...............................................................................  
### ENDS
#...............................................................................
    
    