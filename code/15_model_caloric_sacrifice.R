#...............................................................................
### + GENERIC MODEL OF ACUTE MALNUTRITION AS A FUNCTION OF CALORIC INTAKE ++ ###
#...............................................................................

#...............................................................................
## ------ R SCRIPT TO MODEL CALORIC SACRIFICE BY ADULTS FOR CHILDREN  ------- ##
#...............................................................................


if (! file.exists(paste0(dir_path, "out/15_sacrifice.rds"))) {
#...............................................................................  
### Computing adjustments to child and adult Kcal intake, based on 
    # varying levels of Kcal available and adult Kcal sacrifice
#...............................................................................

  #...................................      
  ## Initialise needed objects

    # Initialise dataframe of 1000 households to be simulated
    hh <- data.frame("hh" = 1:1000)
    for (i in c("m", "age", "sex", "u6mo", "req", "act")) 
      {hh[, paste0(i, 1:max(max_size))] <- NA}
    hh[, c("size", "n_preg", "n_lact", "req_ch", "act_ch", "act_ad", "deficit",
      "act_spare", "act_ch_adj", "act_ad_adj")] <- NA
    cols_hh <- colnames(hh)[2:ncol(hh)]

    # Initialise dataframe of simulated input values and corresponding outputs
    sim <- expand.grid(kcal_capita = seq(100, 3000, 100), 
      req_ch_protected = seq(from = 0.8, to = 1.0, 0.05),
      req_ad_protected = seq(from = 0.3, to = 1.0, 0.05))
    sim[, c("kcal_ch_adj", "kcal_ad_adj")] <- NA
    
    # Loop progress bar   
    pb <- txtProgressBar(min = 1, max = nrow(sim), style = 3)
    
  #...................................
  ## Run simulation for each set of input values
  for (sim_i in 1:nrow(sim)) {

    # Update progress bar
    setTxtProgressBar(pb, sim_i)

    # Specify input values
    req_ch_protected <- sim[sim_i, "req_ch_protected"]
    req_ad_protected <- sim[sim_i, "req_ad_protected"]
    kcal_capita_sim <- sim[sim_i, "kcal_capita"]
    hh[, cols_hh] <- NA

    # Run simulation and collect results
    sim[sim_i, c("kcal_ch_adj", "kcal_ad_adj")] <- f_sacr()
  }
  close(pb)
  
  #...................................      
  ## Save and visualise output
  
    # Save output
    write_rds(sim, paste0(dir_path, "out/15_sacrifice.rds"))
    
    # Graph relationships
    sim$req_ad_sacrificed <- 1 - sim$req_ad_protected
    sim$group <- paste(sim$req_ad_sacrificed, sim$req_ch_protected, sep = "_")
    df <- subset(sim, req_ch_protected >= 0.9 & req_ad_sacrificed >= 0.3 & 
      req_ad_sacrificed <= 0.7)
    ggplot(df, aes(x = kcal_capita, y = kcal_ch_adj, colour = group)) +
      geom_line(linewidth = 1, alpha = 0.75) +
      facet_grid(req_ad_sacrificed ~ req_ch_protected) +
      theme_bw() +
      scale_y_continuous("resulting adjustment to child's intake (multiplier)",
        sec.axis = sec_axis(~ . , 
          name = "proportion of intake that adults are willing to sacrifice", 
          breaks = NULL, labels = NULL)) +
      scale_x_continuous("mean Kcal per person-day at the population level",
        breaks = seq(0, 3000, 200),
        sec.axis = sec_axis(~ . , name = 
        "proportion of child's recommended intake that adults try to safeguard", 
          breaks = NULL, labels = NULL)) +
      scale_color_viridis_d() +
      theme_bw() +
      theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0(dir_path, "out/15_sacrifice.png"),
        dpi = "print", units = "cm", height = 15, width = 22) 
}


#...............................................................................  
### ENDS
#...............................................................................

