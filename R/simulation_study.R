# > sessionInfo() # for reproducibility
# R version 4.2.1 (2022-06-23 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=English_United States.utf8 
# [2] LC_CTYPE=English_United States.utf8   
# [3] LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# attached base packages:
# [1] grid      parallel  splines   stats4    stats     graphics  grDevices
# [8] utils     datasets  methods   base     
# 
# other attached packages:
# [1] optimx_2022-4.30   dfoptim_2020.10-1  lme4_1.1-30       
# [4] Matrix_1.4-1       cowplot_1.1.1      factoextra_1.0.7  
# [7] AlgDesign_1.2.1    DoE.base_1.2-1     conf.design_2.0.0 
# [10] pracma_2.3.8       gridExtra_2.3      cds_1.0.3         
# [13] VGAM_1.1-7         furrr_0.3.1        DHARMa_0.4.5      
# [16] future.callr_0.8.0 future_1.28.0      vegan_2.6-2       
# [19] lattice_0.20-45    permute_0.9-7      ggpubr_0.4.0      
# [22] forcats_0.5.2      stringr_1.4.1      dplyr_1.0.10      
# [25] purrr_0.3.4        readr_2.1.2        tidyr_1.2.0       
# [28] tibble_3.1.8       ggplot2_3.3.6      tidyverse_1.3.2 


source("funcs_sim.R")
set.seed(1)

n_metabolites <- 25 # number of metabolites
n_timepoints <- 11 # number of time-points
n_subjects <- 40 # number of subjects
n_treatments <- 2 # number of treatments
n_obs <- n_subjects*n_treatments*n_timepoints # number of observations
x <- seq(0,10,1) # measurement times

# create random effects
n_re <- 5  # random effects dimension
A <- matrix(runif(n_re^2)-0.5, ncol=n_re) 
Sigma <- (t(A) %*% A)
random.effects <- MASS::mvrnorm(n = n_metabolites*n_subjects, mu = c(0, 0, 0, 0, 0), Sigma = Sigma)

# create random effects
fixed.effects <- readRDS("fixed_effects.Rds")

# create 
df_metab <- data.frame(metabolite = as.factor(c(1:n_metabolites))) %>%
  cbind(fixed.effects)

df_ind <- data.frame(subject = rep(as.factor(c(1:n_subjects)),times=n_metabolites)) %>% 
  mutate(b0 = random.effects[,1]*0.2,
         b1 = random.effects[,2]*0.05,
         b2 = random.effects[,3]*0.005,
         b3 = random.effects[,4]*0.05,
         b4 = random.effects[,5]*0.005)

df_within <- data.frame(subject = as.factor(sort(rep(c(1:n_subjects), n_treatments*n_timepoints))),
                  treatment = rep(c(1,-1), n_subjects*n_timepoints), 
                  tp_L = rep(x, n_treatments*n_subjects)) %>% 
  mutate(tp_Q = tp_L^2) %>%
  arrange(subject,treatment,tp_L)

DF <- merge(df_ind, df_within) %>%
  mutate(metabolite=rep(as.factor(c(1:n_metabolites)),times=n_subjects,each=n_treatments*n_timepoints))
response_matrix <- merge(df_metab,DF)

# compose response matrix
response_matrix <- response_matrix %>%
group_by(metabolite) %>%
mutate(noise = 0.1 * rnorm(n = n_obs),
       concentration = (B0+b0)+(B1+b1)*tp_L+(B2+b2)*tp_Q+(B3+b3)*tp_L*treatment+(B4+b4)*tp_Q*treatment+noise,# simulated data
       Mfull = B1*tp_L+B2*tp_Q+B3*tp_L*treatment+B4*tp_Q*treatment,
       Ma =  B1*tp_L+B2*tp_Q,
       Mab = B3*tp_L*treatment+B4*tp_Q*treatment,
       fxd_model = B0+B1*tp_L+B2*tp_Q+B3*tp_L*treatment+B4*tp_Q*treatment,
       rand_eff = b0+b1*tp_L+b2*tp_Q+b3*tp_L*treatment+b4*tp_Q*treatment,
       mxd_model = (B0+b0)+(B1+b1)*tp_L+(B2+b2)*tp_Q+(B3+b3)*tp_L*treatment+(B4+b4)*tp_Q*treatment,
       mxd_model_noise = (B0+b0)+(B1+b1)*tp_L+(B2+b2)*tp_Q+(B3+b3)*tp_L*treatment+(B4+b4)*tp_Q*treatment+noise) %>% 
  ungroup() %>%
  arrange(metabolite,subject,treatment,tp_L)

simulated_data <- response_matrix %>%
  select(tp_L,tp_Q,treatment,subject,metabolite,concentration)

data_gen_process(response_matrix)

fig_sim_data(response_matrix)

# Compose the encoded effect matrices
# generate the design matrix
design_matrix_enc <- response_matrix %>% 
  select(metabolite,subject,treatment,tp_L,Mfull) %>%
  pivot_wider(names_from = metabolite, values_from = Mfull) %>%
  select(c(subject,treatment,tp_L)) %>%
  `colnames<-` (c("subject","treatment","tp_L"))

# generate efffect matrices
effect_mat <- list()

effect_mat[['Ma']] <- response_matrix %>% 
  select(Ma,metabolite,subject,treatment,tp_L) %>%
  pivot_wider(names_from = metabolite, values_from = Ma) %>% 
  select(!c(subject,treatment,tp_L))

effect_mat[['Mab']] <- response_matrix %>% 
  select(Mab,metabolite,subject,treatment,tp_L) %>%
  pivot_wider(names_from = metabolite, values_from = Mab) %>% 
  select(!c(subject,treatment,tp_L))

# PCA analyses
PCA_res <- effect_mat %>%
  pca_effect_matrices()

# STORE MULTIVARIATE RESULTS OF ENCODED EFFECTS FOR PLOTTING
titles <- c("time", "time_treatment")

scores <- list(PCA_res$Ma$scores,  
               PCA_res$Mab$scores)
loads <- list(PCA_res$Ma$loadings,
              PCA_res$Mab$loadings)
vars <- list(expl_var_from_svd(PCA_res$Ma$singular_values),
             expl_var_from_svd(PCA_res$Mab$singular_values))

design_matrix <- response_matrix %>% 
  filter(metabolite==1) %>%
  select(tp_L,tp_Q,treatment,subject) %>%
  as_tibble()

# create model matrix
X <-
  model.matrix(~ tp_L +
                 tp_Q +
                 tp_L:treatment +
                 tp_Q:treatment,
               design_matrix)

# fit reference models
sim_models <- simulated_data %>% 
  group_by(metabolite) %>%
  group_nest() %>%
  fit_lmm_simple() %>%
  extract_coefs(model_df = ., include_residuals = TRUE)

# Plot individual-specific results
fig_random_effects(sim_models)
fig_random_effects2(sim_models,design_matrix)

# fit bootstrapped models
bootstrap_models <- fit_bootstrap_models(sim_models)

bootstrap_results <- process_bootstrap_models(bootstrap_models,
                                              design_matrix,
                                              c("time","time-treatment"))


# STORE MULTIVARIATE RESULTS OF ESTIMATED EFFECTS FOR PLOTTING
titles_ <- c('res_sim_time_new','res_sim_int_new') 
scores_ <- list(bootstrap_results$Ta,  bootstrap_results$Tab) 
loads_ <- list(bootstrap_results$Pa, bootstrap_results$Pab) 
varis_ <- list(bootstrap_results$perc_a, bootstrap_results$perc_ab) 

# plot population level results
fig_enc_sim_full(titles, scores, loads, vars, design_matrix_enc,
                 titles_, scores_, loads_, varis_)
fig_enc_est()
# END