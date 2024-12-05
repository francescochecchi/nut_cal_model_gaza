#...............................................................................
### + GENERIC MODEL OF ACUTE MALNUTRITION AS A FUNCTION OF CALORIC INTAKE ++ ###
#...............................................................................

#...............................................................................
## - R SCRIPT TO DEVELOP SUB-MODEL OF WEIGHT LOSS DURING INFECTIOUS DISEASE - ##
#...............................................................................

    
#...............................................................................  
### Preparing inputs for modelling weight loss due to disease
#...............................................................................

  #...................................      
  ## Work out generic correspondence between period prevalence and longitudinal
      # prevalence, given different ARI and diarrhoea incidence levels
  if (! file.exists(paste0(dir_path,"out/14_disease_equi.rds"))) {

    # Set incidence values to explore (episodes per child-year)
    inc_set <- seq(0.1, 10, by = 0.1)

    # Initialise output
    equi <- expand.grid(disease = c("ari", "dia"), inc_shape = inc_set)
    equi[, c("period_prev", "point_prev")] <- NA
    
    # Loop progress bar
    pb <- txtProgressBar(min = 1, max = nrow(equi), style = 3)

    # Work out correspondence
    for (i in 1:nrow(equi)) {
      setTxtProgressBar(pb, i)
      equi[i, c("period_prev", "point_prev")] <- 
        f_prev(inc_shape_f = equi[i, "inc_shape"], disease_f=equi[i, "disease"])
    }

    # Save output
    write_rds(equi, paste0(dir_path,"out/14_disease_equi.rds"))

    # Visualise output
    df <- reshape(equi, direction = "long", 
      varying = c("period_prev", "point_prev"),
      idvar = c("disease", "inc_shape"), timevar = "indicator", 
      times = c("period prevalence", "point prevalence"), v.names= "prevalence")
    df$disease <- ifelse(df$disease == "ari", "ARI", "diarrhoea")
    ggplot(df, aes(x = inc_shape, y = prevalence, colour = indicator,
      fill = indicator)) +
      geom_point(alpha = 0.75) +
      geom_line(linewidth = 1) +
      facet_wrap(disease ~., ncol = 2) +
      theme_bw() +
      theme(legend.position = "top", panel.spacing.x = unit(10, "pt")) +
      scale_colour_manual(values = palette_gen[c(5,14)]) +
      scale_fill_manual(values = palette_gen[c(5,14)]) +
      scale_x_continuous("incidence (episodes per child-year)", 
        breaks = seq(0, 10, 1), expand = c(0,0), limits = c(0, NA),) +
      scale_y_continuous("prevalence", labels = scales::percent_format(),
        breaks = seq(0, 1, 0.1), expand = c(0,0), limits = c(0, 1))
    ggsave(paste0(dir_path, "out/14_relationship_inc_prev.png"),
      dpi = "print", units = "cm", height = 12, width = 25)
  }    
  

    
#...............................................................................  
### ENDS
#...............................................................................
    
    
    