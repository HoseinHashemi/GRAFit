GRAFitAddhdr <- function(out_catalog) {

    header <- c("CATAID", "G10ID", "z", "LookBackT", "DEVILS_RA", "DEVILS_DEC",
                "SEMIMAJ", "Ymag", "stellar_mass",
                "GRAFitRA", "GRAFitDEC",
                "ProFound_R50[asec]", "ProFound_R90[asec]", "ProFound_R100[asec]", 
                "ProFound_semimaj[asec]", "ProFound_semimin[asec]",
                "ProFound_mag", "ProFound_magErr", "edge_frac", "asym", "sep", "con",
                "nComp", "BTflux",
                "model_xcen1", "SDxcen1", "MCSExcen1", "model_xcen2", "SDxcen2", "MCSExcen2",
                "model_ycen1", "SDycen1", "MCSEycen1", "model_ycen2", "SDycen2", "MCSEycen2",
                "model_mag1", "SDmag1", "MCSEmag1", "model_mag2", "SDmag2", "MCSEmag2", 
                "model_re1", "SDre1", "MCSEre1", "model_re2", "SDre2", "MCSEre2", 
                "model_n1", "SDnser1", "MCSEnser1", "model_n2", "SDnser2", "MCSEnser2",
                "model_ang1", "SDang1", "MCSEang1", "model_ang2", "SDang2", "MCSEang2",
                "model_axrat1", "SDaxrat1", "MCSEaxrat1", "model_axrat2", "SDaxrat2", "MCSEaxrat2", 
                "logLike", "dof", "AIC", "DIC1", "DIC2", "LML", "fitClass", "elapsed_time", "model_path")
    
  cat <- read.csv(file = out_catalog, header = FALSE, sep = ",")
  
  write.table(cat, file = out_catalog, append = FALSE,
              col.names = header, row.names = FALSE, sep = ",")
  
}
