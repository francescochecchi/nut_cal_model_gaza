#...............................................................................
### + Evolution of child acute malnutrition during war in the Gaza Strip, ++ ###
### + 2023-2024: retrospective estimates and scenario-based projections   ++ ###
#...............................................................................

#...............................................................................
## -- R SCRIPT TO SET KCAL INTAKE SCENARIOS AND PRODUCE GRAPHS FOR PAPER  --- ##
#...............................................................................



#...............................................................................
### Preparatory steps
#...............................................................................

  #...................................      
  ## Install or load required R packages
  if (!"pacman" %in% rownames(installed.packages())){install.packages("pacman")}
  
  pacman::p_load(
    ggplot2,       # Data visualisation
    ggpubr,        # Arrange multiple plots into a single plot
    ggrepel,       # Improve labelling of plots
    readxl,        # Read Excel files
    scales,        # Scale and format data
    tidyverse,     # Tidyverse suite of packages
    viridis,       # Colour palettes
    zoo            # Compute running means
  )       

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
### Setting Kcal per capita by scenario
#...............................................................................

  #...................................      
  ## Prepare inputs and outputs

    # Read and bind scenario data
    best <- readRDS(paste0(dir_path, "kcal_capita_by_area_044.rds"))
    central <- readRDS(paste0(dir_path, "kcal_capita_by_area_033.rds"))
    worst <- readRDS(paste0(dir_path, "kcal_capita_by_area_022.rds"))
    central <- subset(central, week >= as.Date("2024-05-07"))
    worst <- subset(worst, week >= as.Date("2024-05-07"))
    all <- do.call(rbind, list(best, central, worst))
    all$scenario <- c(
      rep("reasonable best", nrow(best)),
      rep("central", nrow(central)), 
      rep("reasonable worst", nrow(worst))
    )
    all[which(all$week < as.Date("2024-05-07") & 
      all$scenario == "reasonable best"), "scenario"] <- "retrospective"
    
    # Add weeks till end of year
    all <- all[, c("week", "area", "est", "scenario")]
    x <- expand.grid(
      week = seq(max(all$week) + 7, as.Date("2024-12-31"), by = 7), 
      scenario = c("reasonable best", "central", "reasonable worst"),
      area = c("north", "south-central")
    )
    x$est <- NA
    x <- x[, colnames(all)]
    all <- rbind(all, x)

  #...................................      
  ## Specify values
        
    # Specify south-central scenario values
    date_ref_start <- as.Date("2024-06-05")
    date_ref_end <- as.Date("2024-09-28")
    dates_ref <- as.Date(date_ref_start : date_ref_end)
    
    x <- which(all$scenario == "reasonable best" & all$week %in% dates_ref &
      all$area == "south-central")
    value <- mean(all[x, "est"])
    x <- which(all$scenario == "reasonable best" & all$week > date_ref_end &
      all$area == "south-central")
    all[x, "est"] <- value * 0.80
    
    x <- which(all$scenario == "central" & all$week %in% dates_ref &
      all$area == "south-central")
    value <- mean(all[x, "est"])
    x <- which(all$scenario == "central" & all$week > date_ref_end &
      all$area == "south-central")
    all[x, "est"] <- value * 0.70
    
    x <- which(all$scenario == "reasonable worst" & all$week %in% dates_ref &
      all$area == "south-central")
    value <- mean(all[x, "est"])
    x <- which(all$scenario == "reasonable worst" & all$week > date_ref_end &
      all$area == "south-central")
    all[x, "est"] <- value * 0.60
    
    # Specify north scenario values
    date_ref_start <- as.Date("2024-04-06")
    date_ref_end <- as.Date("2024-05-04")
    dates_ref <- as.Date(date_ref_start : date_ref_end)
    
    x <- which(all$scenario == "retrospective" & all$week %in% dates_ref &
      all$area == "north")
    value <- mean(all[x, "est"])
    
    x <- which(all$scenario == "reasonable best" & all$week > date_ref_end &
      all$area == "north")
    all[x, "est"] <- value * 1.00
    x <- which(all$scenario == "reasonable best" & 
      all$week > as.Date("2024-10-01") & all$area == "north")
    all[x, "est"] <- value * 0.60
    
    x <- which(all$scenario == "central" & all$week > date_ref_end &
      all$area == "north")
    all[x, "est"] <- value * 0.90
    x <- which(all$scenario == "central" &
      all$week > as.Date("2024-10-01") & all$area == "north")
    all[x, "est"] <- value * 0.50

    x <- which(all$scenario == "reasonable worst" & all$week > date_ref_end &
      all$area == "north")
    all[x, "est"] <- value * 0.80
    x <- which(all$scenario == "reasonable worst" & 
      all$week > as.Date("2024-10-01") & all$area == "north")
    all[x, "est"] <- value * 0.40
    all$scenario <- factor(all$scenario, levels = c("retrospective",
      "reasonable best", "central", "reasonable worst"))
    all <- all[order(all$scenario, all$area, all$week), ]

  #...................................      
  ## Save and visualise output
        
    # Save (and graph)can be directly pasted onto 'scenarios' tab of 
        # 'crisis_specs.xlsx'
    write.csv(all, paste0(dir_path, "31_kcal_capita_by_week.csv"), row.names=F)
    
    # Graph
    ggplot(all, aes(x = week, y = est, colour = scenario, linetype = scenario))+
      geom_step(linewidth = 1) +
      scale_colour_manual("scenario", values = palette_gen[c(1,5,9,15,12,3)]) +
      scale_linetype_manual("scenario", 
        values = c("solid","11","31","22","12","33")) +
      facet_grid(area ~ .) +
      theme_bw() +
      scale_x_date("date", date_labels = "%b-%Y", breaks = "1 month", 
        expand = c(0,0)) +
      scale_y_continuous("Kcal per capita") +
      theme(legend.position = "top", axis.text.x = element_text(angle = 30,
        hjust = 0.25, vjust = 0.25))
    ggsave(paste0(dir_path, "31_kcal_scenarios.png"), units = "cm",
      dpi = "print", width = 23, height = 17)
      

    
#...............................................................................
### Visualising GAM ground prevalence estimates
#...............................................................................

  #...................................      
  ## Read ground metadata and manage variables
    
    # Read dataset and manage variables
    df <- read_excel(paste0(dir_path, "gaza_gam_ground.xlsx"))
    df[which(df$governorate == "Gaza"), "governorate"] <- "Gaza City"
    df$date_1 <- as.Date(df$date_1)
    df$date_2 <- as.Date(df$date_2)
    df$date_mid <- as.Date(df$date_1 + as.numeric(df$date_2 - df$date_1) / 2)
    df[which(df$governorate == "Deir el Balah"), "governorate"] <- 
      "south-central (Deir el Balah)"
    df[which(df$governorate == "Gaza City"), "governorate"] <- 
      "north (Gaza City)"
    df[which(df$governorate == "Khan Younis"), "governorate"] <- 
      "south-central (Khan Younis)"
    df[which(df$governorate == "North Gaza"), "governorate"] <- 
      "north (North Gaza)"
    df[which(df$governorate == "Rafah"), "governorate"] <- 
      "south-central (Rafah)"
    
    df$governorate <- factor(df$governorate, 
      levels = sort(unique(df$governorate)))
    df <- df[order(df$governorate), ]
    df <- subset(df, status_kids != "sick")
 
    
  #...................................      
  ## Smooth trends
    
    # Set up output
    sm <- expand.grid(
      date_mid = as.Date(as.Date("2024-01-01") : as.Date("2024-11-30")),
      governorate = unique(df$governorate)
    )
    x <- length(unique(df$governorate))
    sm$t_0623 <- rep(1:(nrow(sm) / x), x)
    sm$t_0659 <- rep(1:(nrow(sm) / x), x)
    sm <- merge(sm, df, by = c("governorate", "date_mid"), all.x = T)
    
    for (i in unique(sm$governorate)) {
      x <- df[which(df$governorate == i & ! is.na(df$gam_muacz_0623_est)), ]
      sm[which(sm$governorate == i & (sm$date_mid < min(x$date_mid) | 
        sm$date_mid > max(x$date_mid))), "no_0623"] <- 1
      x <- df[which(df$governorate == i & ! is.na(df$gam_muacz_0659_est)), ]
      sm[which(sm$governorate == i & (sm$date_mid < min(x$date_mid) | 
        sm$date_mid > max(x$date_mid))), "no_0659"] <- 1
    }
    sm[which(sm$no_0623 == 1 & is.na(sm$gam_muacz_0623_est)), "t_0623"] <- NA
    sm[which(sm$no_0659 == 1 & is.na(sm$gam_muacz_0659_est)), "t_0659"] <- NA
            
    # Smooth
    x <-by(sm, sm$governorate, function(xx) {
      x1 <- xx[complete.cases(xx$date_mid, xx$gam_muacz_0623_est), ]
      x1 <- x1[order(x1$date_mid), ]
      fit <- smooth.spline(x = x1$t_0623, y = x1$gam_muacz_0623_est, 
        w = x1$sampsi, spar = 0.7)
      xx$gam_muacz_0623_smooth <- NA
      xx[which(! is.na(xx$t_0623)), "gam_muacz_0623_smooth"] <- 
        predict(fit, xx[which(! is.na(xx$t_0623)), "t_0623"])$y
      x1 <- xx[complete.cases(xx$date_mid, xx$gam_muacz_0659_est), ]
      x1 <- x1[order(x1$date_mid), ]
      fit <- smooth.spline(x = x1$t_0659, y = x1$gam_muacz_0659_est, 
        w = x1$sampsi, spar = 0.7)
      xx$gam_muacz_0659_smooth <- NA
      xx[which(! is.na(xx$t_0659)), "gam_muacz_0659_smooth"] <- 
        predict(fit, xx[which(! is.na(xx$t_0659)), "t_0659"])$y
      return(xx)      
    })
    sm <- do.call(rbind, x)

    
  #...................................      
  ## Reshape and graph
    
    # Reshape long
    sm <- sm[, c("governorate", "date_mid", "sampsi", "gam_muacz_0623_est",
      "gam_muacz_0623_smooth", "gam_muacz_0659_est", "gam_muacz_0659_smooth")]  
    x1 <- sm[, c("governorate", "date_mid", "sampsi", "gam_muacz_0623_est",
      "gam_muacz_0623_smooth")]
    colnames(x1) <- c("governorate", "date_mid", "sampsi", "est", "smooth")
    x1$age <- "GAM prevalence (MUACZ < 2SD), 6-23mo"
    x2 <- sm[, c("governorate", "date_mid", "sampsi", "gam_muacz_0659_est",
      "gam_muacz_0659_smooth")]
    colnames(x2) <- c("governorate", "date_mid", "sampsi", "est", "smooth")
    x2$age <- "GAM prevalence (MUACZ < 2SD), 6-59mo"      
    sm <- rbind(x1, x2)
    sm <- sm[order(as.character(sm$governorate)),]
    
    # Plot non-sick estimates by governorate (long)
    ggplot(sm, aes(y = est,
      x = date_mid, colour = governorate, fill = governorate, size = sampsi,
      group = governorate)) +
      geom_point(alpha = 0.5) +
      geom_line(aes(y = smooth, x = date_mid,
        colour = governorate, group = governorate), 
        linewidth = 0.5, linetype = "11") +
      scale_y_continuous("prevalence", limits = c(0, 0.35), labels = percent) +
      scale_x_date("date", date_labels = "%b-%Y", breaks = "1 month", 
        expand = c(0,0)) +
      scale_size_continuous("sample size") +
      facet_grid(governorate~age) +
      theme_bw() +
      scale_colour_manual("governorate", values = palette_gen[c(1,4,7,11,14)]) +
      scale_fill_manual("governorate", values = palette_gen[c(1,4,7,11,14)]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        legend.position = "top") +
      guides(colour = "none", fill = "none")
    ggsave(paste0(dir_path, "31_trends_gam_muacz_long.png"), units = "cm",
      dpi = "print", width = 21, height = 26)
    
    # Plot non-sick estimates by governorate (wide)
    ggplot(sm, aes(y = est,
      x = date_mid, colour = governorate, fill = governorate, size = sampsi,
      group = governorate)) +
      geom_point(alpha = 0.5) +
      geom_line(aes(y = smooth, x = date_mid,
        colour = governorate, group = governorate), 
        linewidth = 0.5, linetype = "11") +
      scale_y_continuous("prevalence", limits = c(0, 0.35), labels = percent) +
      scale_x_date("date", date_labels = "%b-%Y", breaks = "1 month", 
        expand = c(0,0)) +
      scale_size_continuous("sample size") +
      facet_grid(age~governorate) +
      theme_bw() +
      scale_colour_manual("governorate", values = palette_gen[c(1,4,7,11,14)]) +
      scale_fill_manual("governorate", values = palette_gen[c(1,4,7,11,14)]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        legend.position = "top") +
      guides(colour = "none", fill = "none")
    ggsave(paste0(dir_path, "31_trends_gam_muacz_wide.png"), units = "cm",
      dpi = "print", width = 32, height = 17)
    
  
#...............................................................................
### Visualising test model output
      # first need to run model with `crisis_specs_test.xlsx`;
      # for convenience, output has been renamed 'out_agg_test_a.rds' 
      # and 'out_agg_test_b.rds', and is included within this sub-folder
#...............................................................................

  #...................................      
  ## Test run plots 

    # Plot A
    out_agg <- readRDS(paste0(dir_path, "22_out_agg_test_a.rds"))
    df <- subset(out_agg, outcome == "GAM prevalence (6-59mo)")
    df$day <- as.integer(df$date - min(df$date))
    pl_a <- ggplot(df, aes(x = day,
      y = median, colour = scenario, fill = scenario, linetype = scenario)) +
      geom_line() +
      geom_ribbon(aes(ymin = lci95, ymax = uci95), alpha = 0.2, colour=NA) +
      theme_bw() +
      scale_x_continuous("day", breaks = seq(0,100, 14), expand = c(0,0)) +
      scale_y_continuous(
        "prevalence of global acute malnutrition (children 6 to 59mo)", 
        limits=c(0,1), breaks = seq(0,1,0.1), labels=percent, expand = c(0,0)) +
      scale_linetype_manual("scenario",
        values = c("solid","11","31","22","12","33")) +
      scale_colour_manual("scenario", values = palette_gen[c(1,5,9,15,12,3)]) +
      scale_fill_manual("scenario", values = palette_gen[c(1,5,9,15,12,3)]) +
      theme(legend.position = "top") +
      annotate("text", x = 1, y = 0.9, size = 3, hjust = 0, lineheight = 1,
        label = " Scenarios:\n 1) 2028 Kcal/day (recommended), other factors at pre-crisis levels\n 2) 1700 Kcal/day\n 3) 1700 Kcal/day + adult caloric sacrifice = 20%\n 4) 1700 Kcal/day + adult caloric sacrifice = 20% + ARI/diarrhoea excess point prevalence = 30%\n 5) 1700 Kcal/day + adult caloric sacrifice = 20% + ARI/diarrhoea excess point prevalence = 30% + SAM/MAM treatment coverage = 100%"
      )
    ggsave(paste0(dir_path, "31_plot_test_a.png"), units = "cm",
      dpi = "print", width = 23, height = 17)
    saveRDS(pl_a, paste0(dir_path, "31_plot_test_a.rds"))
    
    # Plot B
    out_agg <- readRDS(paste0(dir_path, "22_out_agg_test_b.rds"))
    df <- subset(out_agg, outcome == "mean WHZ")
    df$day <- as.integer(df$date - min(df$date))
    df$sacrifice <- (as.numeric(gsub("%","",as.character(df$scenario))))/100
    df$sacrifice <- factor(df$sacrifice, levels = c(0.10,0.20,0.30,0.40),
      labels = c("10%","20%","30%","40%"))
    pl_b <- ggplot(df, aes(x = day,
      y = median, colour = sacrifice, fill = sacrifice, linetype = sacrifice)) +
      geom_line(alpha = 0.5) +
      geom_ribbon(aes(ymin = lci95, ymax = uci95), alpha = 0.2, colour=NA) +
      theme_bw() +
      scale_x_continuous("day", breaks = seq(0,100, 14), expand = c(0,0)) +
      scale_y_continuous("mean weight-for-height Z-score", limits=c(-3.5,0.1),
        breaks = seq(-3.5,0.5,0.5)) +
      scale_linetype_manual("percent sacrificed",
        values = c("solid","11","31","22","12","33")) +
      scale_colour_manual("percent sacrificed", 
        values = palette_gen[c(1,5,9,15,12,3)]) +
      scale_fill_manual("percent sacrificed", 
        values = palette_gen[c(1,5,9,15,12,3)]) +
      theme(legend.position = "top") +
      annotate("text", x = 45, y = 0.1, size = 3, hjust = 0.5, lineheight = 1,
        label = " All scenarios assume 1700 Kcal/day")
    ggsave(paste0(dir_path, "31_plot_test_b.png"), units = "cm",
      dpi = "print", width = 23, height = 17)
    saveRDS(pl_b, paste0(dir_path, "31_plot_test_b.rds"))
    
    # Combination plot
    pl_a <- readRDS(paste0(dir_path, "31_plot_test_a.rds"))
    pl_b <- readRDS(paste0(dir_path, "31_plot_test_b.rds"))
    ggarrange(pl_a + theme(axis.text.x = element_blank(), 
      axis.title.x = element_blank()), pl_b, nrow=2, ncol=1, labels=c("A","B"))
    ggsave(paste0(dir_path, "31_plot_test_combi.png"), units = "cm",
      dpi = "print", width = 23, height = 28)
        
    
#...............................................................................
### Preparing other graphs for model presentations
#...............................................................................
    
  #...................................      
  ## Graph of weight-for-height distribution before and during crisis
      
  mean_b <- -0.3
  sd_b = 1.1
  mean_c <- -1.3
  sd_c = 1.1
      
  df <- data.frame(whz = rep(seq(-5, 5, 0.001), 2))
  df$period <- c(rep("baseline", nrow(df)/2), rep("crisis", nrow(df)/2))
  df$density <- NA
  x <- which(df$period == "baseline")
  df[x, "density"] <- dnorm(df[x, "whz"], mean = mean_b, sd = sd_b)
  x <- which(df$period == "crisis")
  df[x, "density"] <- dnorm(df[x, "whz"], mean = mean_c, sd = sd_c)
  df$mean_whz <- ifelse(df$period == "crisis", mean_c, mean_b)

  ggplot(df, aes(x = whz, y = density, colour = period, fill = period)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymax = density), ymin = 0,
      alpha = 0) +
    geom_ribbon(data = df[which(df$whz < -2), ], aes(ymax = density), ymin = 0,
      alpha = 0.5) +
    geom_ribbon(data = df[which(df$whz < -3), ], aes(ymax = density), ymin = 0,
      alpha = 0.8) +
    scale_x_continuous("weight-for-height Z-score", breaks = seq(-5, 5, 1)) +
    scale_colour_manual(values = palette_gen[c(5,10)]) +
    scale_fill_manual(values = palette_gen[c(5,10)]) +
    geom_vline(aes(xintercept = 0), linetype = "22", colour = palette_gen[16],
      linewidth = 1) +
    geom_segment(aes(x = mean_whz, xend = mean_whz, 
      y = 0, yend = max(density)), linetype = "22", linewidth = 1) +
    # geom_segment(aes(x = mean_c, xend = mean_c, 
    #   y = 0, yend = max(df[which(period == "crisis"), "density"])), 
    #   colour = palette_gen[10], linetype = "22", linewidth = 1) +
    facet_wrap(period ~., nrow = 2) +
    theme_bw() +
    theme(legend.position = "none")
  
  ggsave(paste0(dir_path, "31_whz_distributions.png"), dpi = "print", 
    units = "cm", height = 15, width = 20)
  

  #...................................      
  ## Graph of individual calibration cumulative log-normal function
  
    # Prepare set of exemplar age-dependent parameter values (girls)
    pars <- data.frame(
      centile = paste0(c(0,5,20,50,80,95,99), "%"),
      centile_num = c(0,5,20,50,80,95,99),
      par0 = c(0.081, 0.039, 0.018, 0, 0.024, 0.051, 0.075),
      par1 = c(4.2, 3.9, 3.6, 0, 5.1, 5.4, 5.1)
    )    
    
    # Output dataset with all ages
    df <- expand.grid(age = 0:1826, centile = paste0(c(0,5,20,50,80,95,99),"%"))
    df <- merge(df, pars, by = "centile", all.x = TRUE)
    df <- df[order(df$centile, df$age), ]
    
    # Cumulative log-normal parameters
    df$par1b <- dLOGNO(df$age / 100, df$par1)
    df$theta0 <- df$par0
    x <- by(df, df$centile, function(x) {x$par0 * x$par1b / max(x$par1b)})
    df$theta1 <- df$par0 * as.vector(unlist(x))
    df[df$centile_num >= 50, c("theta0", "theta1")] <- 
      - df[df$centile_num >= 50, c("theta0", "theta1")]
    
    # Adjustment
    df$theta <- df$theta0 + df$theta1
    
    # Plot  
    ggplot(df, aes(x = age / 30.44, y = theta, colour = centile, 
      linetype = centile)) +
      geom_line(linewidth = 1) +
      scale_colour_manual("growth percentile",
        values = palette_gen[c(1,3,6,9,11,13,16)]) +
      scale_linetype_discrete("growth percentile") +
      theme_bw() +
      theme(legend.position = "top") +
      scale_y_continuous("relative adjustment to energy need") +
      scale_x_continuous("age (months)", breaks = seq(0,60, 6))
    ggsave(paste0(dir_path, "31_calibration_adjustments.png"), 
      dpi = "print", units = "cm", height = 12, width = 20)

    
    
#...............................................................................
### ENDS
#...............................................................................
    
    