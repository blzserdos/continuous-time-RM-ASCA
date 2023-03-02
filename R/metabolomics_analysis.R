source("funcs_metab.R")
lipo_dat_ <- readRDS("lipo_dat_VLDL_LDL.Rds") # load data
lipo_dat_$treatment <- recode(lipo_dat_$treatment, "1" = "A", "2" = "B") # recode treatment
lipo_dat_$concentration_raw <- lipo_dat_$concentration # create var with original outcome (for bootstrapping)

lipo_dat <- SDscale(lipo_dat_) # create scaled outcome var

fig_data(lipo_dat)

# specify design properties of data
Design_Info <- data_design(
  n_treatments       = 2,
  n_timepoints       = 11,
  timepoints         = c(
    "0",
    "30",
    "60",
    "90",
    "120",
    "150",
    "180",
    "210",
    "240",
    "270",
    "300"
  ) %>%
    factor(x      = .,
           levels = .),
  treatments         = c("A", "B") %>%
    factor(x      = .,
           levels = .),
  n_PC               = 3 ,
  n_variables        = 36,
  n_subjects         = 20
) 

# create design matrix
design_matrix <- gen_design(Design_Info) %>%
  arrange(subject, treatment, t)

Design_Info$design <- design_matrix

model_mat <- merge(design_matrix, 
                   lipo_dat, 
                   by.x = c("t", "treatment", "subject"), 
                   by.y = c("Time", "treatment", "subject")) %>% 
  arrange(metabolite, subject, treatment, t) 


start_time <- Sys.time()
models <- model_mat %>% 
  group_by(metabolite) %>%
  group_nest() %>% 
  fit_lmm() %>%  
  extract_coefs(model_df = ., include_residuals = T)
print(paste0("finished in ",format(Sys.time()-start_time, digits=2)))


start_time <- Sys.time()
bootstrap_models <- fit_bootstrap_models(models,
                                         Design_Info)
print(paste0("finished in ",format(Sys.time()-start_time, digits=2)))

bootstrap_results <- process_bootstrap_models(bootstrap_models,
                                              Design_Info,
                                              c("time","interaction"))

titles <- c('melc_full', 'melc_time','melc_int') 
Scores <- list(bootstrap_results$Tfull, bootstrap_results$Ta,  bootstrap_results$Tab) 
Loads <- list(bootstrap_results$Pfull, bootstrap_results$Pa, bootstrap_results$Pab) 
for(i in 1:3){
  s <- str_split(Loads[[i]]$var,'_')
  Loads[[i]]$frac <- sapply(s, '[', 3)
  Loads[[i]]$size <- sapply(s, '[', 1)
  Loads[[i]]$lipoclass <- sapply(s, '[', 2)
  Loads[[i]] <- Loads[[i]] %>% arrange(size, frac, lipoclass)
  
  Loads[[i]]$frac <- factor(Loads[[i]]$frac, levels = c("CE","FC","PL","TG"))
  Loads[[i]]$size <- factor(Loads[[i]]$size, levels = c("XXL","XL","L","M","S","XS"))
  Loads[[i]]$lipoclass <- factor(Loads[[i]]$lipoclass, levels = c("VLDL","LDL"))

}
Variances <- list(bootstrap_results$perc_full, bootstrap_results$perc_a, bootstrap_results$perc_ab) 

fig_est_scores_loadings_JK_re2(titles, Scores, Loads, Variances)

fig_random_effects(models, Design_Info$design)

