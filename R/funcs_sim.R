# load packages
library("utils")
library("stats")
library("tidyverse")
library("ggpubr")
library("vegan")
library("future")
library("future.callr")
library("DHARMa")
library("furrr")
library("VGAM")
library("cds")
library("gridExtra")
library("pracma")
library("DoE.base")
library("AlgDesign")
library("factoextra")
library("cowplot")
library("lme4")
library("dfoptim")
library("optimx")

fit_lmm_simple <- function(long_response_df) {
  # fit LMM
  
  cl <- parallel::makeCluster(6)
  plan(cluster, workers = cl)
  
  res <- long_response_df %>% 
    mutate(M = future_map(data, function(x)
        lme4::lmer(
          concentration ~ 
            tp_L + 
            tp_Q + 
            tp_L:treatment +
            tp_Q:treatment + 
            (
              tp_L +
              tp_Q +
              tp_L:treatment +
              tp_Q:treatment |
              subject
            ),
          data = x,
          na.action=na.exclude,
          REML = TRUE)))
  parallel::stopCluster(cl)
  future:::ClusterRegistry("stop")
  plan(sequential)
  return(res)
}

extract_coefs <- function(model_df,
                           include_residuals = FALSE){
  # extract LMM coefficients
  
  output<- model_df %>%
    mutate(
      # Fixed Effects
      fixed = map(M, function(x) lme4::fixef(x)),
      # Random effects
      random = map(M, function(x) lme4::ranef(x)))

  if(include_residuals){
    output <- output %>%
      mutate(
        # Residuals
        resids = map(M,function(x)  residuals(x)))}
  return(output)
}

perform_pca <- function(input_matrix){
  # perform PCA via SVD
  
  USV <- svd(input_matrix)
  output <- list()
  output$scores <- (USV$u*sqrt(dim(input_matrix)[1]-1))[,1:3]
  output$loadings <- (USV$v %*% diag(USV$d) / sqrt(dim(input_matrix)[1]-1))[,1:3]
  output$singular_values = USV$d[1:3]
  return(output)
}


pca_effect_matrices <- function(input_matrix) {
  # Center data and perform PCA
  
  return(
    lapply(input_matrix,
           FUN = "col_center") %>%
      lapply(FUN = "perform_pca")
  )
}

expl_var_from_svd <- function(singular_values,
                              decimals = 1){
  # Calculate explained variance
  
  return(round(singular_values^2/sum(singular_values^2)*100,decimals))
}

col_center <- function(data_matrix){
  # center data
  
  return(sweep(x      = data_matrix,
               MARGIN = 2,
               STATS  = colMeans(data_matrix),
               FUN    = "-"))}

# Orthogonal procrustes rotation
orthprocr_pca <- function(target,query){
  sol <- query
  sol$loadings <- cds::orthprocr(Z = target$loadings,
                                 X = query$loadings)$XQ
  sol$scores <- sol$scores %*% (cds::orthprocr(Z = target$loadings,
                                               X = query$loadings)$Q %>% solve())
  return(sol)
}

fit_bootstrap_models <- function(fitted_models){
  
  # Recreate data from fitted models
  original_data <- do.call("rbind", fitted_models$data) %>%
    mutate(metabolite = rep(fitted_models$metabolite, each = nrow(fitted_models$data[[1]])),) %>%
    arrange(tp_L, treatment, subject)
 
  bootstrap_models <- list()
  bootstrap_metabolites  <- list()
  bootstrap_PCA <- list()
  bootstrap_res <- list()

  # set up bootstrap details
  no_samples <- 100
  bootstraps <- rsample::group_bootstraps(original_data,
                            group=subject,
                            times=no_samples)

  # Iterate over bootstrap folds
  for(i in 1:no_samples
  ){
    print(paste0("starting fold ",i))
    fold_start_time <- Sys.time()
    bootstrap_models[[i]]<-
      as.data.frame(bootstraps$splits[[i]]) %>% # BS
      group_by(subject,treatment,tp_L) %>%
      mutate(dupeID = row_number()) %>%
      ungroup() %>%
      group_by(metabolite) %>%
      arrange(dupeID,subject,treatment,tp_L) %>%
      group_nest() %>%
      fit_lmm_simple(long_response_df = .) %>%
      extract_coefs(model_df = .,
                    include_residuals = T) %>% 
      filter(!is.na(M)) %>%
      select(-c(data))
    bootstrap_metabolites[[i]]<- bootstrap_models[[i]]$metabolite
    print(paste0("fold ",i," finished in ",format(Sys.time()-fold_start_time, digits=2)))
  }
  
  shared_fitted_metabolites <- Reduce(intersect,bootstrap_metabolites)
  
  output <- list()
  output$metabolites <- shared_fitted_metabolites
  output$models<- bootstrap_models
  output$full_model <- fitted_models
  output$bootstraps <- bootstraps
  return(output)}


process_bootstrap_models <- function(bootstrap_models_list,
                                     design_matrix,
                                     effects){
  shared_fitted_metabolites <- bootstrap_models_list$metabolites
  bootstrap_models <- bootstrap_models_list$models
  fitted_models <- bootstrap_models_list$full_model 
  subject_group_df <- bootstrap_models_list$subject_groups
  bootstraps <- bootstrap_models_list$bootstraps
  bootstrap_res <- list()

  effect_matrices <- compose_effmat(fitted_models,design_matrix)
  
  reference_pca <- effect_matrices %>%
    pca_effect_matrices()

  for(i in 1:length(bootstrap_models)){

    fold_coef_mat <-compose_effmat(bootstrap_models[[i]],design_matrix)
    
    print('clear2')

    fold_pca_res <- fold_coef_mat %>%
      pca_effect_matrices()
    
    # perform procrustes rotation
    fold_pca_res$Ma <- orthprocr_pca(reference_pca$Ma,
                                     fold_pca_res$Ma)
    fold_pca_res$Mab <- orthprocr_pca(reference_pca$Mab,
                                      fold_pca_res$Mab)
    fold_pca_res$Mfull <- orthprocr_pca(reference_pca$Mfull,
                                        fold_pca_res$Mfull)
    # Calculate explained variance
    fold_expl_a <- tibble(fold = i,
                          PC = seq(1:length(fold_pca_res$Ma$singular_values)),
                          perc_expl = fold_pca_res$Ma$singular_values %>% expl_var_from_svd())%>%
      filter(PC %in% seq(1:3))
 
    fold_expl_ab <- tibble(fold = i,
                           PC = seq(1:length(fold_pca_res$Mab$singular_values)),
                           perc_expl = fold_pca_res$Mab$singular_values %>% expl_var_from_svd())%>%
      filter(PC %in% seq(1:3))
    
    fold_expl_full <- tibble(fold = i,
                             PC = seq(1:length(fold_pca_res$Mfull$singular_values)),
                             perc_expl = fold_pca_res$Mfull$singular_values %>% expl_var_from_svd())%>%
      filter(PC %in% seq(1:3))

    fold_Ta <- fold_pca_res$Ma$scores %>%
      as_tibble(.name_repair) %>%
      cbind(design_matrix) %>%
      pivot_longer(cols=-colnames(design_matrix),
                   names_to = "PC",
                   values_to = "score",
                   names_prefix = "V") %>%
      mutate(PC = as.factor(PC),
             treatment = as.factor(treatment)) %>%
      filter(PC %in% seq(1:3))

    fold_Tab <- fold_pca_res$Mab$scores %>%
      as_tibble(.name_repair)%>%
      cbind(design_matrix) %>%
      pivot_longer(cols=-colnames(design_matrix),
                   names_to = "PC",
                   values_to = "score",
                   names_prefix = "V")%>%
      mutate(PC = as.factor(PC),
             treatment = as.factor(treatment)) %>%
      filter(PC %in% seq(1:3))
   
    fold_Tfull <- fold_pca_res$Mfull$scores %>%
      as_tibble(.name_repair)%>%
      cbind(design_matrix) %>%
      pivot_longer(cols=-colnames(design_matrix),
                   names_to = "PC",
                   values_to = "score",
                   names_prefix = "V")%>%
      mutate(PC = as.factor(PC),
             treatment = as.factor(treatment)) %>%
      filter(PC %in% seq(1:3))

    fold_Pa <- fold_pca_res$Ma$loadings %>%
      as_tibble(.name_repair)%>%
      mutate(var = shared_fitted_metabolites,
             fold = i)%>%
      pivot_longer(col=-var,
                   names_to = "PC",
                   values_to = "loading",
                   names_prefix = "V")%>%
      filter(PC %in% seq(1:3))

    fold_Pab <- fold_pca_res$Mab$loadings %>%
      as_tibble(.name_repair)%>%
      mutate(var = shared_fitted_metabolites,
             fold = i)%>%
      pivot_longer(col=-var,
                   names_to = "PC",
                   values_to = "loading",
                   names_prefix = "V") %>%
      filter(PC %in% seq(1:3))
   
    fold_Pfull <- fold_pca_res$Mfull$loadings %>%
      as_tibble(.name_repair)%>%
      mutate(var = shared_fitted_metabolites,
             fold = i)%>%
      pivot_longer(col=-var,
                   names_to = "PC",
                   values_to = "loading",
                   names_prefix = "V") %>%
      filter(PC %in% seq(1:3))

    if(i == 1){
      bootstrap_res$Ta <- fold_Ta
      bootstrap_res$Pa <- fold_Pa
      bootstrap_res$Tab <- fold_Tab
      bootstrap_res$Pab <- fold_Pab
      bootstrap_res$Tfull <- fold_Tfull
      bootstrap_res$Pfull <- fold_Pfull
      bootstrap_res$perc_a <- fold_expl_a
      bootstrap_res$perc_ab <- fold_expl_ab
      bootstrap_res$perc_full <- fold_expl_full
    }else{
      bootstrap_res$Ta <- rbind(bootstrap_res$Ta,fold_Ta)
      bootstrap_res$Pa <- rbind(bootstrap_res$Pa,fold_Pa)
      bootstrap_res$Tab <- rbind(bootstrap_res$Tab,fold_Tab)
      bootstrap_res$Pab <- rbind(bootstrap_res$Pab,fold_Pab)
      bootstrap_res$Tfull <- rbind(bootstrap_res$Tfull,fold_Tfull)
      bootstrap_res$Pfull <- rbind(bootstrap_res$Pfull,fold_Pfull)
      bootstrap_res$perc_a <- rbind(bootstrap_res$perc_a,fold_expl_a)
      bootstrap_res$perc_ab<- rbind(bootstrap_res$perc_ab,fold_expl_ab)
      bootstrap_res$perc_full<- rbind(bootstrap_res$perc_full,fold_expl_full)
      }
  }

  bootstrap_res$reference <- reference_pca
  bootstrap_res$effect_matrices <- effect_matrices
  return(bootstrap_res)}

compose_effmat <-
  function(ref_mod, design_matrix,
           include_residuals = F) {
    output_ <- list()
    for (i in 1:length(ref_mod$metabolite)) {
        design_mat <- design_matrix %>% 
          model.matrix(
            ~ tp_L +
              tp_Q +
              tp_L:treatment +
              tp_Q:treatment
            ,
            .
          )

      # Intercept
      output_$M0[[paste0('metab_', i)]] <-
        as.matrix(select(as.data.frame(design_mat), contains("Intercept"))) %*% as.matrix(select(as.data.frame(t(ref_mod$fixed[[i]])), contains('Intercept')))
      # Time effect
      output_$Ma[[paste0('metab_', i)]] <-
        as.matrix(select(as.data.frame(design_mat), ends_with(c(
          "L", "Q"
        )))) %*% t(as.matrix(select(as.data.frame(
          t(ref_mod$fixed[[i]])
        ), ends_with(
          c("L", "Q")
        ))))
      # Time-Treatment interaction effect
      output_$Mab[[paste0('metab_', i)]] <-
        as.matrix(select(as.data.frame(design_mat), contains("treatment"))) %*% t(as.matrix(select(
          as.data.frame(t(ref_mod$fixed[[i]])), contains('treatment')
        )))
      output_$Mfull[[paste0('metab_', i)]] <-
        as.matrix(select(as.data.frame(design_mat), !contains("Intercept"))) %*% t(as.matrix(select(
          as.data.frame(t(ref_mod$fixed[[i]])), !contains('Intercept')
        )))
    }
    
    output <- list()
    output$M0 <- do.call("cbind", output_$M0)
    output$Ma <- do.call("cbind", output_$Ma)
    output$Mab <- do.call("cbind", output_$Mab)
    output$Mfull <- do.call("cbind",output_$Mfull)
    if (include_residuals) {
      output$Mres <-
        do.call("cbind", ref_mod$resids)
    }
    return(output)
  }

create_re_effect_matrices <- function(design_matrix, sim_models, coeffs_, i) {
  
  des <- list()
  
  X <- design_matrix %>%
    filter(subject==1) %>%
    model.matrix(
      ~ tp_L +
        tp_Q +
        tp_L:treatment +
        tp_Q:treatment     ,
      .
    ) %>% 
    as.data.frame() %>%
    select(coeffs_) %>%
    as.matrix()
  
  gamma <- sim_models$random[[i]]$subject %>% 
    select(coeffs_) %>%
    t()
  
  des$X <- X
  des$gamma <- gamma
  return(des)
}


score_plot_enc_avg <- function(design_matrix, Scores, Variances,label,j,i){
  
  df <- design_matrix
  n_PC <- 3
  df$score <- Scores[[j]][,i]
  
  S <- df %>%
    group_by(tp_L) %>%
    ggplot(aes(x = tp_L,
               y = score)) +
    geom_hline(yintercept = 0, linetype='dotted', color="grey") +
    geom_line(aes(x = tp_L,
                  y = score)) +

    scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
    scale_x_continuous(breaks = c(0,2,4,6,8,10))+
    
    theme_classic() +
    theme(
      plot.margin = margin(0, 0, 0, 0, unit = "cm"),
      axis.text.x = element_blank(),
      axis.text.y = element_text(color = "black"),
      axis.title.y = element_text(vjust = -1),
      axis.ticks = element_line(color = "black"),
      legend.background = element_rect(color = NA, fill = "transparent"),
      legend.box.background = element_rect(color = NA, fill = "transparent"),
      legend.key = element_rect(color = NA, fill = "transparent"),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    labs(x = "", y = label[i])
}

score_plot_enc <- function(design_matrix, Scores, Variances,label,j,i){
  
  df <- design_matrix
  n_PC <- 3
  df$score <- Scores[[j]][,i]
  
  S <- df %>%
    #filter(PC == i) %>%
    mutate(treatment = recode(as.factor(treatment), `1` = "B", `-1` = "A")) %>%
    group_by(tp_L, treatment) %>%
    ggplot(aes(x = tp_L,
               y = score)) +
    geom_hline(yintercept = 0, linetype='dotted', color="grey") +
    geom_line(aes(x = tp_L,
                  y = score,
                  linetype = treatment)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
    scale_x_continuous(breaks = c(0,2,4,6,8,10))+
    theme_classic() +
    theme(
      plot.margin = margin(0, 0, 0, 0, unit = "cm"),
      axis.text.x = element_blank(),
      axis.text.y = element_text(color = "black"),
      axis.title.y = element_text(vjust = -1),
      axis.ticks = element_line(color = "black"),
      legend.background = element_rect(color = NA, fill = "transparent"),
      legend.box.background = element_rect(color = NA, fill = "transparent"),
      legend.key = element_rect(color = NA, fill = "transparent"),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    labs(x = "", y = label[i])
}

loads_plot_enc <- function(Loads, Variances, j, i){
  
  pos <- as.factor(seq(1,dim(Loads[[1]])[1],1))
  L <- tibble(var = seq(1,25),
            loading = Loads[[j]][,i]) %>%
    group_by(var) %>%
    mutate(color_code = as.factor(ceiling(as.numeric(var)/5))) %>%
    ungroup() %>%

    ggplot(aes(x=var,
               y=loading, 
               fill=color_code))+ 
    geom_bar(
      stat = "identity",
      position = "dodge",
      width=0.8
    )+ 
    scale_fill_brewer(type = "qual", palette = 'Set3',direction = -1) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01)
    ) +
    geom_hline(yintercept = 0, linetype='dotted', color="grey") +
    
    theme_classic()+
    theme(plot.margin = margin(0,0.2,0,0, unit = "cm")) +

    labs(x="",
         y="")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y = element_text(color="black"),
          axis.ticks = element_line(color="black"),
          legend.background = element_rect(color = NA, fill = "transparent"),
          legend.box.background = element_rect(color = NA, fill = "transparent"),
          legend.key = element_rect(color = NA, fill = "transparent"),
          legend.title = element_blank(),
          legend.text = element_text(color='black', size=10),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA),
          strip.background = element_blank(),
          strip.text.x = element_blank()) + 
    scale_x_discrete(breaks = seq(1,25,2), limits = pos)
}

score_plot_avg <- function(Scores, Variances,label,j,i){
  S <- Scores[[j]] %>%
    filter(PC == i) %>%
    group_by(tp_L) %>%
    mutate(
      quantile_025 = quantile(score, .025),
      quantile_500 = quantile(score, .500),
      quantile_975 = quantile(score, .975)
    ) %>%
    
    ggplot(aes(x = tp_L,
               y = quantile_500)) +
    geom_hline(yintercept = 0, linetype='dotted', color="grey") +
    geom_ribbon(aes(ymin = quantile_025, 
                    ymax = quantile_975),
                alpha=0.3
    ) +
    geom_line(aes(x = tp_L,
                  y = quantile_500)) +

    scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
    scale_x_continuous(breaks = c(0,2,4,6,8,10))+
    theme_classic() +
    theme(
      plot.margin = margin(0, 0, 0, 0, unit = "cm"),
      axis.text.x = element_blank(),
      axis.title.y = element_text(vjust = -1),
      axis.text.y = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),
      legend.background = element_rect(color = NA, fill = "transparent"),
      legend.box.background = element_rect(color = NA, fill = "transparent"),
      legend.key = element_rect(color = NA, fill = "transparent"),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    labs(x = "", y = label[i])
}

score_plot <- function(Scores, Variances,label,j,i){
  S <- Scores[[j]] %>%
    filter(PC == i) %>%
    mutate(treatment = recode(as.factor(treatment), `1` = "B", `-1` = "A")) %>%
    group_by(tp_L, treatment) %>%
    mutate(
      quantile_025 = quantile(score, .025),
      quantile_500 = quantile(score, .500),
      quantile_975 = quantile(score, .975)
    ) %>%
    
    ggplot(aes(x = tp_L,
               y = quantile_500)) +
    geom_hline(yintercept = 0, linetype='dotted', color="grey") +
    geom_ribbon(aes(ymin = quantile_025, 
                    ymax = quantile_975,
                    fill = treatment),
                alpha=0.5
    ) +
    geom_line(aes(x = tp_L,
                  y = quantile_500,
                  linetype = treatment)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
    scale_x_continuous(breaks = c(0,2,4,6,8,10))+
    theme_classic() +
    theme(
      plot.margin = margin(0, 0, 0, 0, unit = "cm"),
      axis.text.x = element_blank(),
      axis.title.y = element_text(vjust = -1),
      axis.text.y = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),
      legend.background = element_rect(color = NA, fill = "transparent"),
      legend.box.background = element_rect(color = NA, fill = "transparent"),
      legend.key = element_rect(color = NA, fill = "transparent"),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    labs(x = "", y = label[i])
}

loads_plot_sim <- function(Loads, Variances, j, i){
  
  pos <- as.factor(seq(1,length(unique(as.numeric(Loads[[1]]$var))),1))
  
  L<-Loads[[j]] %>%
    filter(PC==i) %>%
    group_by(var) %>%
    mutate(quantile_025 = quantile(loading,.025),
           quantile_500 = quantile(loading,.500),
           quantile_975 = quantile(loading,.975)) %>%
    mutate(color_code = as.factor(ceiling(as.numeric(var)/5))) %>%
    ungroup() %>%
    select(-loading) %>%
    ggplot(aes(x=var, 
               y=quantile_500, 
               fill=color_code))+ 
    geom_bar(
      stat = "identity",
      position = "dodge",
      width = 0.8
    )+ 

    scale_fill_brewer(type = "qual", palette = 'Set3',direction = -1) +
    geom_hline(yintercept = 0, linetype='dotted', color="grey") +
    geom_errorbar(aes(ymin=quantile_025,
                      ymax=quantile_975),
                  size=0.3,
                  width=0.8,
                  stat = 'identity',
                  position = 'dodge')+
    theme_classic()+
    theme(plot.margin = margin(0,0.2,0,0, unit = "cm")) +
    labs(x="",
         y="")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y = element_text(color="black"),
          axis.ticks = element_line(color="black"),
          legend.background = element_rect(color = NA, fill = "transparent"),
          legend.box.background = element_rect(color = NA, fill = "transparent"),
          legend.key = element_rect(color = NA, fill = "transparent"),
          legend.title = element_blank(),
          legend.text = element_text(color='black', size=10),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA),
          strip.background = element_blank(),
          strip.text.x = element_blank()) + 
    scale_x_discrete(breaks = seq(1,25,2), limits = pos)
}

loads_plot <- function(Loads, Variances, j, i){
  
  # bar plot positions lipo_dat_VLDL_LDL.Rds
  pos <- c(
    "XXL_VLDL_CE",
    "XXL_VLDL_FC",
    "XXL_VLDL_PL",
    "XXL_VLDL_TG",
    "XL_VLDL_CE",
    "XL_VLDL_FC",
    "XL_VLDL_PL",
    "XL_VLDL_TG",
    "L_VLDL_CE",
    "L_VLDL_FC",
    "L_VLDL_PL",
    "L_VLDL_TG",
    "M_VLDL_CE",
    "M_VLDL_FC",
    "M_VLDL_PL",
    "M_VLDL_TG",
    "S_VLDL_CE",
    "S_VLDL_FC",
    "S_VLDL_PL",
    "S_VLDL_TG",
    "XS_VLDL_CE",
    "XS_VLDL_FC",
    "XS_VLDL_PL",
    "XS_VLDL_TG",
    "L_LDL_CE",
    "L_LDL_FC",
    "L_LDL_PL",
    "L_LDL_TG",
    "M_LDL_CE",
    "M_LDL_FC",
    "M_LDL_PL",
    "M_LDL_TG",
    "S_LDL_CE",
    "SLDL_FC",
    "S_LDL_PL",
    "S_LDL_TG"
  )
    

  L<-Loads[[j]] %>%
    filter(PC==i)%>%
    group_by(var)%>%
    mutate(quantile_025 = quantile(loading,.025),
           quantile_500 = quantile(loading,.500),
           quantile_975 = quantile(loading,.975)) %>%
    select(-loading)%>%
    distinct()%>%
      ggplot(aes(x=var,
                 y=quantile_500,
                 fill=paste0(size, lipoclass)))+ 
    geom_bar(stat="identity",
             position="dodge"
    )+ 
    scale_fill_brewer(type = "qual", palette = 'Set3',direction = -1) +
    geom_hline(yintercept = 0, linetype='dashed', color="grey") +
    geom_errorbar(aes(ymin=quantile_025,
                      ymax=quantile_975),
                  stat = 'identity',
                  position = 'dodge')+
    theme_classic()+
    theme(plot.margin = margin(0,0,0,0, unit = "cm")) +
    scale_y_continuous(
      labels = scales::number_format(accuracy = 0.1)) +
    labs(x="",
         y="")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y = element_text(color="black"),
          axis.ticks = element_line(color="black"),
          legend.background = element_rect(color = NA, fill = "transparent"),
          legend.box.background = element_rect(color = NA, fill = "transparent"),
          legend.key = element_rect(color = NA, fill = "transparent"),
          legend.title = element_blank(),
          legend.text = element_text(color='black', size=10),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA),
          strip.background = element_blank(),
          strip.text.x = element_blank()) + 
    scale_x_discrete(limits = pos)
}

fig_random_effects <- function(sim_models) {
  # create random effect effect matrices, then do bar plot, biplot
  
  n_obs <- dim(sim_models$M[[1]]@frame)[1]
  n_ind <- length(unique(sim_models$M[[1]]@frame$subject))
  
  time_timetreat <- list()
  timetreat <- list()
  time <- list()
  for (i in 1:25) {
    time_timetreat[[paste0('metab_',i)]] <-
      as.matrix(select(as.data.frame(model.matrix(sim_models$M[[i]])), contains("tp")))[seq(1, n_obs, n_ind),] %*% t(as.matrix(select(
        sim_models$random[[i]]$subject, contains('tp')
      )))
    timetreat[[paste0('metab_',i)]] <-
      as.matrix(select(as.data.frame(model.matrix(sim_models$M[[i]])), contains('treatment')))[seq(1, n_obs, n_ind),] %*% t(as.matrix(select(
        sim_models$random[[i]]$subject, contains('treatment')
      )))
    time[[paste0('metab_',i)]] <-
      as.matrix(select(as.data.frame(model.matrix(sim_models$M[[i]])), ends_with(c(
        "L", "Q", "C", "d", "L2", "Q2", "C2"
      ))))[seq(1, n_obs, n_ind),] %*% t(as.matrix(select(
        sim_models$random[[i]]$subject, ends_with(c("L", "Q", "C", "d", "L2", "Q2", "C2"))
      )))
  }
  
  timedf <- t(do.call(rbind, time))
  timetreatdf <- t(do.call(rbind, timetreat))
  time_timetreatdf <- t(do.call(rbind, time_timetreat))
  
  pcTime <- perform_pca(scale(timedf, center = T, scale = F))
  pcTimeTreat <- perform_pca(scale(timetreatdf, center = T, scale = F))
  pcTime_TimeTreat <- perform_pca(scale(time_timetreatdf, center = T, scale = F))
  
  cnames <- 
    c(
      "individual",
      paste0(c("Time/PC1 (","Time/PC2 ("), expl_var_from_svd(pcTime$singular_values)[1:2],"%)"),
      paste0(c("Time-Treatment/PC1 (","Time-Treatment/PC2 ("),expl_var_from_svd(pcTimeTreat$singular_values)[1:2],"%)"),
      paste0(c("Time + Time-Treatment/PC1 (","Time + Time-Treatment/PC2 ("),expl_var_from_svd(pcTime_TimeTreat$singular_values)[1:2],"%)")
    )
  
  scDF <- as_tibble(cbind(
    seq(1, n_ind, 1),
    pcTime$scores[, 1:2],
    pcTimeTreat$scores[, 1:2],
    pcTime_TimeTreat$scores[, 1:2]
  ))
  
  colnames(scDF) <- cnames
  
  scoresDF <-
    pivot_longer(
      scDF,
      cols = !individual,
      names_to = "PC",
      values_to = "Score"
    ) %>%
    separate(PC,
             into=c('Effect','PC'),
             sep = "/") %>%
    separate(PC,
             into=c("PC","var"),
             sep = " ") %>%
    mutate(Effect = factor(Effect, levels = c("Time", "Time-Treatment", "Time + Time-Treatment")),
           PC = factor(PC, levels = c("PC1", "PC2"))
           )
  
  p <- ggplot(scoresDF, aes(x=as.factor(individual), y=Score)) +
    geom_bar(stat="identity", width=0.8) +
    scale_x_discrete(guide = guide_axis(n.dodge=2)) +
    labs(x = "Individual")
  p <- p + facet_grid(Effect~PC) +
    theme(axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color="black"),
          axis.ticks = element_line(color="black"),
          panel.border = element_rect(color="black", fill = NA),
          strip.background = element_rect(color = "black")
    ) + 
    geom_text(
      mapping = aes(x = Inf, y = -Inf, label = var),
      hjust   = 1.5,
      vjust   = -1,
      size=3
    )
  show(p)
  
  ggsave(
    filename = paste0(Sys.Date(), "_sim_random_effects.png"),
    device = "png",
    plot = p,
    width = 6,
    height = 5
  )
}

fig_random_effects2 <- function(sim_models, design_matrix) {
  # create random effect effect matrices, then create bar plot
  
  n_ind <- length(unique(sim_models$M[[1]]@frame$subject))
  
  time_timetreat <- list()
  timetreat <- list()
  time <- list()
  for (i in 1:25) {
    des_t <-
      create_re_effect_matrices(design_matrix, sim_models, all_of(c("tp_L", "tp_Q")), i)
    time[[paste0('metab_', i)]] <- des_t$X %*% des_t$gamma
    
    des_ttr <-
      create_re_effect_matrices(design_matrix, sim_models, all_of(c("tp_L:treatment", "tp_Q:treatment")), i)
    timetreat[[paste0('metab_', i)]] <- des_ttr$X %*% des_ttr$gamma
    
    des_t_ttr <-
      create_re_effect_matrices(design_matrix, sim_models, all_of(c("tp_L", "tp_Q", "tp_L:treatment", "tp_Q:treatment")), i)
    time_timetreat[[paste0('metab_', i)]] <- des_t_ttr$X %*% des_t_ttr$gamma
    
  }
  
  timedf_t <- do.call(cbind, lapply(time, function(x) Reshape(x,880,1)))
  timetreatdf_t <- do.call(cbind, lapply(timetreat, function(x) Reshape(x,880,1)))
  time_timetreatdf_t <- do.call(cbind, lapply(time_timetreat, function(x) Reshape(x,880,1)))
  
  pcTime_t <- perform_pca(scale(timedf_t, center = T, scale = F))
  pcTimeTreat_t <- perform_pca(scale(timetreatdf_t, center = T, scale = F))
  pcTime_TimeTreat_t <- perform_pca(scale(time_timetreatdf_t, center = T, scale = F))
  
  cnames_t <- 
    c("individual",
      "time",
      "treatment",
      paste0(c("Time/PC1 (","Time/PC2 ("), expl_var_from_svd(pcTime_t$singular_values)[1:2],"%)"),
      paste0(c("Time-Treatment/PC1 (","Time-Treatment/PC2 ("),expl_var_from_svd(pcTimeTreat_t$singular_values)[1:2],"%)"),
      paste0(c("Time + Time-Treatment/PC1 (","Time + Time-Treatment/PC2 ("),expl_var_from_svd(pcTime_TimeTreat_t$singular_values)[1:2],"%)")
    )
  
  scDF_t <- as_tibble(cbind(
    as.matrix(design_matrix[,c("subject","tp_L","treatment")]),
    pcTime_t$scores[, 1:2],
    pcTimeTreat_t$scores[, 1:2],
    pcTime_TimeTreat_t$scores[, 1:2]
  ))
  
  colnames(scDF_t) <- cnames_t
  
  scoresDF_t <-
    pivot_longer(
      scDF_t,
      cols = !c(individual,time,treatment),
      names_to = "PC",
      values_to = "Score"
    ) %>%
    separate(PC,
             into=c('Effect','PC'),
             sep = "/") %>%
    separate(PC,
             into=c("PC","var"),
             sep = " ") %>%
    mutate(Effect = factor(Effect, levels = c("Time", "Time-Treatment", "Time + Time-Treatment")),
           PC = factor(PC, levels = c("PC1", "PC2"))
    )
  
  
  scoresDF_t_wide <- scoresDF_t %>%
    pivot_wider(names_from = PC, values_from = c('Score','var')) %>%
    mutate(individual = as.factor(individual),
           time = as.numeric(time),
           Score_PC1 = as.numeric(Score_PC1),
           Score_PC2 = as.numeric(Score_PC2),
           treatment = as.factor(treatment))
  
  eff <- c("Time + Time-Treatment", "Time", "Time-Treatment")
  pp <- list()
  pp_ <- list()
  for (i in 1:3) {
    df <- subset(scoresDF_t_wide, Effect == eff[i])
    
    pp[[i]] <- ggplot(df, aes(x = time, y = Score_PC1, group=interaction(individual,treatment), linetype=treatment)) +
      geom_line() +
      labs(
        x = "time, min",
        y = paste0("PC1", " ", df$var_PC1[[1]]),
        title = df$Effect[[1]]
      ) +
      theme_bw() +
      theme(
        plot.margin = margin(0, 0, 0, 0, unit = "cm"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        legend.background = element_rect(color = NA, fill = "transparent"),
        legend.box.background = element_rect(color = NA, fill = "transparent"),
        legend.key = element_rect(color = NA, fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        aspect.ratio = 1
      )
    
    pp_[[i]] <- ggplot(df, aes(x = time, y = Score_PC2, group=interaction(individual,treatment), linetype=treatment)) +
      geom_line() +
      labs(
        x = "time, min",
        y = paste0("PC2", " ", df$var_PC2[[1]])
      ) +
      theme_bw() +
      theme(
        plot.margin = margin(0, 0, 0, 0, unit = "cm"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        legend.background = element_rect(color = NA, fill = "transparent"),
        legend.box.background = element_rect(color = NA, fill = "transparent"),
        legend.key = element_rect(color = NA, fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        aspect.ratio = 1
      )
  }
  
  pl_re <- cowplot::plot_grid(
    pp[[1]],
    NULL,
    pp[[2]],
    NULL,
    pp[[3]]+theme(legend.position = c(0.2,0.1),
                  legend.text = element_text(color='black', 
                                             size=6,
                                             margin=margin(t=1,r=1,b=1,l=1, unit = "pt")),
                  legend.margin=margin(c(0,0,0,0)),
                  legend.spacing.x = unit(0, "mm"),
                  legend.spacing.y = unit(0, "mm"),
                  legend.key.size = unit(0.25, "cm"),
                  legend.title=element_text(size=8)),
    NULL,NULL,NULL,NULL,NULL,
    pp_[[1]],
    NULL,
    pp_[[2]],
    NULL,
    pp_[[3]],
    ncol = 5,
    nrow = 3,
    rel_widths = c(1,0.1,1,0.1,1),
    rel_heights = c(1,-0.5,1),
    align = 'hv',
    axis = 'lb'
  )
  
  ggsave(filename=paste0(Sys.Date(),"_sim_random_effects_ind_t",
                         ".png"),
         #fallback_resolution=600,
         plot=pl_re,
         width = 8,
         height= 8)
  
  timedf <- t(do.call(rbind, time))
  timetreatdf <- t(do.call(rbind, timetreat))
  time_timetreatdf <- t(do.call(rbind, time_timetreat))
  
  pcTime <- perform_pca(scale(timedf, center = T, scale = F))
  pcTimeTreat <- perform_pca(scale(timetreatdf, center = T, scale = F))
  pcTime_TimeTreat <- perform_pca(scale(time_timetreatdf, center = T, scale = F))
  
  cnames <- 
    c(
      "individual",
      paste0(c("Time/PC1 (","Time/PC2 ("), expl_var_from_svd(pcTime$singular_values)[1:2],"%)"),
      paste0(c("Time-Treatment/PC1 (","Time-Treatment/PC2 ("),expl_var_from_svd(pcTimeTreat$singular_values)[1:2],"%)"),
      paste0(c("Time + Time-Treatment/PC1 (","Time + Time-Treatment/PC2 ("),expl_var_from_svd(pcTime_TimeTreat$singular_values)[1:2],"%)")
    )
  
  scDF <- as_tibble(cbind(
    seq(1, n_ind, 1),
    pcTime$scores[, 1:2],
    pcTimeTreat$scores[, 1:2],
    pcTime_TimeTreat$scores[, 1:2]
  ))
  
  colnames(scDF) <- cnames
  
  scoresDF <-
    pivot_longer(
      scDF,
      cols = !individual,
      names_to = "PC",
      values_to = "Score"
    ) %>%
    separate(PC,
             into=c('Effect','PC'),
             sep = "/") %>%
    separate(PC,
             into=c("PC","var"),
             sep = " ") %>%
    mutate(Effect = factor(Effect, levels = c("Time", "Time-Treatment", "Time + Time-Treatment")),
           PC = factor(PC, levels = c("PC1", "PC2"))
    )
  
  p1 <- ggplot(scoresDF, aes(x=as.factor(individual), y=Score)) +
    geom_bar(stat="identity", width=0.8) +
    scale_x_discrete(guide = guide_axis(n.dodge=2)) +
    labs(x = "Individual")
  p1 <- p1 + facet_grid(Effect~PC) +
    theme(axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color="black"),
          axis.ticks = element_line(color="black"),
          panel.border = element_rect(color="black", fill = NA),
          strip.background = element_rect(color = "black")
    ) + 
    geom_text(
      mapping = aes(x = Inf, y = -Inf, label = var),
      hjust   = 1.5,
      vjust   = -1,
      size=3
    )
  show(p1)
  
  ggsave(
    filename = paste0(Sys.Date(), "sim_RE_scores.png"),
    device = "png",
    plot = p1,
    width = 6,
    height = 5
  )
  
  scoresDF_wide <- scoresDF %>%
    pivot_wider(names_from = PC, values_from = c('Score','var'))
  
  eff <- c("Time", "Time-Treatment")
  pp <- list()
  for (i in 1:2) {
    df <- subset(scoresDF_wide, Effect == eff[i])
    
    pp[[i]] <- ggplot(df, aes(x = Score_PC1, y = Score_PC2)) +
      geom_hline(yintercept = 0, linetype='dotted', color="grey") +
      geom_vline(xintercept = 0, linetype='dotted', color="grey") +
      geom_point() +
      ggrepel::geom_text_repel(
        aes(label = individual),
        size=3,
        box.padding   = 0.35,
        point.padding = 0.5,
        min.segment.length = 0.1
      ) +
      labs(
        x = paste0("PC1", " ", df$var_PC1[[1]]),
        y = paste0("PC2", " ", df$var_PC2[[1]])#,
      ) +
      theme_bw() +
      theme(
        plot.margin = margin(0, 0, 0, 0, unit = "cm"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.background = element_rect(color = NA, fill = "transparent"),
        legend.box.background = element_rect(color = NA, fill = "transparent"),
        legend.key = element_rect(color = NA, fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        aspect.ratio = 1
      )
  }
  
  pl_re <- cowplot::plot_grid(
    pp[[1]],
    NULL,
    pp[[2]],
    ncol = 3,
    labels = c("a","","b"),
    rel_widths = c(1,0.1,1),
    align = 'hv',
    axis = 'lb'
  )
  
  ggsave(filename=paste0(Sys.Date(),"_sim_random_effects_ind",
                         ".png"),
         #fallback_resolution=600,
         plot=pl_re,
         width = 5.33,
         height= 3)
  
  
  pdf <- subset(scoresDF_wide, Effect=="Time-Treatment") %>%
    mutate(individual = as.factor(individual))
  p <- ggplot(pdf, aes(x = Score_PC1, y = Score_PC2)) +
    geom_hline(yintercept = 0, linetype='dotted', color="grey") +
    geom_vline(xintercept = 0, linetype='dotted', color="grey") +
    geom_point() +
    ggrepel::geom_text_repel(
      aes(label = individual),
      size=3,
      box.padding   = 0.35,
      point.padding = 0.5,
      min.segment.length = 0.1
    ) +
    geom_point(data = subset(pdf, individual %in% c(7,8,18)),
               aes(x = Score_PC1,
                   y = Score_PC2,
                   color = individual)) +
    labs(
      x = paste0("PC1", " ", pdf$var_PC1[[1]]),
      y = paste0("PC2", " ", pdf$var_PC2[[1]]),
      title = pdf$Effect[[1]]
    ) +
    theme_bw() +
    theme(
      plot.margin = margin(0, 0, 0, 0, unit = "cm"),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),
      legend.background = element_rect(color = NA, fill = "transparent"),
      legend.box.background = element_rect(color = NA, fill = "transparent"),
      legend.key = element_rect(color = NA, fill = "transparent"),
      legend.position = "none",
      panel.background = element_rect(fill = "transparent"),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      aspect.ratio = 1
    )
  
  ggsave(filename=paste0(Sys.Date(),"_sim_re_highlighted",
                         ".png"),
         #fallback_resolution=600,
         plot=p,
         width = 4,
         height= 4)
  
}

fig_sim_data <- function(data){

  data$treatment <-
    recode(as.factor(data$treatment), `1` = "B", `-1` = "A")

  p1 = ggplot(data=data, aes(x = tp_L, 
                             y = concentration, 
                             linetype=treatment,
                             color=treatment,
                             group = interaction(subject, treatment))) +
    geom_line(alpha=0.6) +
    
    theme_bw() +
    facet_wrap(~metabolite, nrow = 3, ncol=9, 
               strip.position ="top") +
    labs(x = "time", y= "concentration") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    theme(
      plot.margin = unit(c(0.06,0.1,0.06,0.06), "cm"),
      axis.ticks = element_line(color="black"),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      legend.background = element_rect(color = NA, fill = "transparent"),
      legend.box.background = element_rect(color = NA, fill = "transparent"),
      legend.key = element_rect(color = NA, fill = "transparent"),
      legend.title = element_text(size=8,margin = margin(b = 1, unit = "pt")),
      legend.text = element_text(size=8, margin = margin(b = 1, unit = "pt")),
      legend.position = c(0.9,0.14),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
      strip.text = element_text(size=8)) +
    guides(linetype = guide_legend(nrow=2))
  
  
  ggsave(paste0(Sys.Date(),"sim_data.png"),
         p1, width = 7, height = 3.2)
  
  p2 = ggplot(data=data, aes(x = tp_L, 
                             y = concentration, 
                             linetype=treatment,
                             group = interaction(subject, treatment))) +
    geom_line(color="gray") +
    scale_x_continuous(breaks = c(0,2,4,6,8,10))+
    geom_line(data=subset(data, subject %in% c(36,37)),aes(color=subject))+
    theme_bw() +
    facet_wrap(~metabolite, nrow = 5, ncol=5) +
    labs(x = "time", y= "standardized concentration") +
    
    theme(
      plot.margin = unit(c(0.06,0.06,0.06,0.06), "cm"),
      axis.ticks = element_line(color="black"),
      axis.text = element_text(color="black"),
      axis.line = element_line(color="black"),
      legend.background = element_rect(color = NA, fill = "transparent"),
      legend.box.background = element_rect(color = NA, fill = "transparent"),
      legend.key = element_rect(color = NA, fill = "transparent"),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
      strip.background = element_blank(),
      strip.text = element_text(size=8),
      legend.direction = "vertical", 
      legend.box = "horizontal")
  show(p2)
  
  ggsave(paste0(Sys.Date(),"_sim_data_highlighted.png"),
         p2, width = 7, height = 7)
  
}

data_gen_process <- function(data){
  
  df <- subset(data, metabolite==1 & subject %in% c(1,2)) %>%
    select(metabolite, tp_L, treatment, subject, fxd_model, mxd_model, mxd_model_noise) %>%
    pivot_longer(cols=c(fxd_model,mxd_model,mxd_model_noise), names_to="var") %>%
    mutate(treatment = recode(as.factor(treatment), `1` = "B", `-1` = "A"))
  
  theme_options <- theme_classic() +
    theme(
      plot.margin = margin(0.2, 0.2, 0.2, 0.2, unit = "cm"),
      axis.text.x = element_text(color="black"),
      axis.line.x = element_line(color="black"),
      axis.text.y = element_text(color = "black"),
      axis.line.y = element_line(color="black"),
      axis.ticks = element_line(color = "black"),
      legend.title = element_text(size=8,margin = margin(b = 1, unit = "pt")),
      legend.text = element_text(size=8,margin = margin(b = 1, unit = "pt")),
      legend.spacing = unit(0.01, 'cm'),
      # legend.key.size = unit(0.3, "cm"),
      legend.key.width = unit(0.5, "cm"),
      legend.key.height = unit(0.3, "cm"),
      legend.background = element_rect(color = NA, fill = "transparent"),
      legend.box.background = element_rect(color = NA, fill = "transparent"),
      legend.key = element_rect(color = NA, fill = "transparent"),
      legend.position = "none",
      panel.background = element_rect(fill = "transparent"),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
      strip.background = element_blank(),
      strip.text.x = element_blank()
    )
  
  
  p1 = ggplot(subset(df,var=="fxd_model")) +
    geom_line(aes(
      x = tp_L,
      y = value,
      linetype = treatment
    )) + 
    labs(x="time", y="concentration") +
    scale_linetype_manual(values = c('solid','dashed'), labels = c("A", "B"))+
    scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
    ylim(4.4,5.4) +
    theme_options
  
  p2 <- ggplot(subset(df,var=="mxd_model")) +
    geom_line(aes(
      x = tp_L,
      y = value,
      linetype = treatment,
      color = subject,
      group = interaction(subject, treatment)
    )) + 
    scale_linetype_manual(values = c('solid','dashed'), labels = c("A", "B")) +
    scale_color_manual(values = c("#fc8d62", "#8da0cb")) +
    labs(x="time", y="") +
    scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
    ylim(4.4,5.4) +
    theme_options + 
    theme(legend.position = c(0.77,0.71),
          legend.box = "horizontal")
  
  p3 <- ggplot(subset(df,var=="mxd_model_noise")) +
    geom_line(aes(
      x = tp_L,
      y = value,
      linetype = treatment,
      color = subject,
      group = interaction(subject, treatment)
    )) + 
    scale_linetype_manual(values = c('solid','dashed'), labels = c("A", "B"))+
    scale_color_manual(values = c("#fc8d62", "#8da0cb")) +
    scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
    ylim(4.4,5.4) +
    labs(x="time", y="") +
    theme_options 

  p = cowplot::plot_grid(
    p1,p2,p3,
    ncol=3,
    align='hv',
    axis='lb',
    labels=c("a","b","c")
    )
  
  ggsave(filename=paste0(Sys.Date(),"_sim_example_metab",
                         ".png"),
         plot=p,
         width = 8,
         height = 2)
}


fig_enc_sim_full <- function(titles, Scores, Loads, Variances, design_matrix_enc,
                             titles_est, Scores_est, Loads_est, Variances_est){
  list_of_plots_vert <- list()
  list_of_plots_horz <- list()
  panel <-list()
  panel[[1]] <- c("a","c")
  panel[[2]] <- c("b","d")

  for (j in 1:length(Scores)){
    n_PC <- 3
    labels_ <- list()
    labels_est <-list()
    for (i in 1:n_PC){
      labels_[[i]] <- paste0("PC",i,"\n(",
                             Variances[[j]][i] %>% 
                               round(1) %>%
                               format(nsmall = 1),
                             "%)")
      
      labels_est[[i]] <- paste0(
        "PC",
        i,
        "\n(",
        filter(Variances_est[[j]],
               PC == i)$perc_expl %>%
          quantile(.025) %>%
          round(1) %>%
          format(nsmall = 1),
        "-",
        filter(Variances_est[[j]],
               PC == i)$perc_expl %>%
          quantile(.975) %>%
          round(1) %>%
          format(nsmall = 1),
        "%)") 
    }
    
    labelss <- unlist(labels_, recursive=FALSE)
    labelss_est <- unlist(labels_est, recursive=FALSE)
    
    Stitle <- ggdraw() +
      draw_label(
        "Scores",
        x = 0.45,
        hjust = 0.5
      ) +
      theme(plot.margin = margin(0, 0, 0, 0))
    
    Ltitle <- ggdraw() +
      draw_label(
        "Loadings",
        x = 0.45,
        hjust = 0.5
      ) +
      theme(plot.margin = margin(0, 0, 0, 0))
    
    title_row <- cowplot::plot_grid(
      Stitle,
      Ltitle,
      ncol = 2,
      nrow = 1,
      align = 'hv',
      axis = 'lb'
    )
    
    # ENCODED
    S1a <- score_plot_enc_avg(design_matrix_enc, Scores, Variances, labelss, j, 1) +
      theme(
        legend.position = 'none',
        axis.title.y = element_text(margin=margin(r=7)),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank())
    
    S1 <- score_plot_enc(design_matrix_enc, Scores, Variances, labelss, j, 1) +
      theme(
        legend.position = 'none',
        axis.title.y = element_text(margin=margin(r=7)),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank())

    L1 <- loads_plot_enc(Loads, Variances, j, 1) +

      theme(
        axis.line.x = element_blank(),
        legend.position = 'none') #+
    
    S2a <- score_plot_enc_avg(design_matrix_enc, Scores, Variances, labelss, j, 2) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = 'none',
            axis.title.y = element_text(margin=margin(r=7)),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank()) #+
    
    S2 <- score_plot_enc(design_matrix_enc, Scores, Variances, labelss, j, 2) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = 'none',
            axis.title.y = element_text(margin=margin(r=7)),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank()) #+

    L2 <- loads_plot_enc(Loads, Variances, j, 2) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.line.x = element_blank(),
            legend.position = 'none') #+

    S3a <- score_plot_enc_avg(design_matrix_enc, Scores, Variances, labelss, j, 3) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = 'none',
            axis.title.y = element_text(margin=margin(r=7)),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank())# +
    
    S3 <- score_plot_enc(design_matrix_enc, Scores, Variances, labelss, j, 3) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = 'none',
            axis.title.y = element_text(margin=margin(r=7)),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank())# +

    L3 <- loads_plot_enc(Loads, Variances, j, 3) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.line.x = element_blank(),
            legend.position = 'none')#+
    
    if(j==1){
      S1 <- S1a
      S2 <- S2a
      S3 <- S3a
      
      S1 <- S1 + theme(
        legend.title = element_text(size=8,margin = margin(b = 1, unit = "pt")),
        legend.text = element_text(size=8, margin = margin(b = 1, unit = "pt")),
        legend.spacing = unit(0.01, 'cm'),
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.position = c(0.3, 0.9)) +
        guides(linetype = guide_legend(override.aes = list(size = 0.4)))
    }
    
    if (j==2){
      S2_ <- S2 + 
        labs(x='time') +
        theme(axis.text.x = element_text(color="black"),
              axis.ticks.x = element_line(color="black"),
              axis.line.x = element_line(color="black"),
              axis.title.y = element_text(margin=margin(r=7)))
      
      L2_ <- L2 +
        labs(x='metabolite') +
        theme(axis.ticks.x = element_line(color="black"), 
              axis.line.x = element_line(color="black"),
              axis.text.x = element_text(color='black'),
              legend.position = 'none')
    } else {
      S2_ <- S2 + 
        labs(x='') +
        theme(axis.ticks.x = element_blank(), 
              axis.line.x = element_blank(),
              axis.text.x = element_text(color='black', size=5),
              legend.position = 'none') +
        scale_x_continuous(labels = function(breaks) {rep_along(breaks, "")})
      
      L2_ <- L2 +
        labs(x='') +
        theme(axis.ticks.x = element_blank(), 
              axis.line.x = element_blank(),
              axis.text.x = element_text(color='black'),
              legend.position = 'none') +
        scale_x_discrete(labels = rep("", 25), breaks = 1:25)
    }
    
    # ESTIMATED
    
    S1aest <- score_plot_avg(Scores_est, Variances_est, labelss_est, j, 1) +
      theme(
        legend.position = "none",
        axis.title.y = element_text(margin=margin(r=7)),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()) #+
    
    S1est <- score_plot(Scores_est, Variances_est, labelss_est, j, 1) +
      theme(
        legend.position = "none",
        axis.title.y = element_text(margin=margin(r=7)),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()) #+

    L1est <- loads_plot_sim(Loads_est, Variances_est, j, 1) +
      theme(
        axis.line.x = element_blank(),
        legend.position = 'none')+
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01),
                         breaks=c(-0.1,0.0,0.1))
    
    S2aest <- score_plot_avg(Scores_est, Variances_est, labelss_est, j, 2) +
      theme(legend.position = 'none',
            axis.title.y = element_text(margin=margin(r=7)),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank())# +
    
    S2est <- score_plot(Scores_est, Variances_est, labelss_est, j, 2) +
      theme(legend.position = 'none',
            axis.title.y = element_text(margin=margin(r=7)),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank())# +
    L2est <- loads_plot_sim(Loads_est, Variances_est, j, 2) +
      guides(fill = guide_legend(nrow = 3, byrow=TRUE,
                                 override.aes = list(size = 0.1, 'mm'))) +
      theme(axis.line.x = element_blank(),
            legend.position = 'none')+
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
                         #breaks=c(-0.05,0.0,0.05))
    S3aest <- score_plot_avg(Scores_est, Variances_est, labelss_est, j, 3) +
      theme(legend.position = 'none',
            axis.title.y = element_text(margin=margin(r=7)),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank())
    
    S3est <- score_plot(Scores_est, Variances_est, labelss_est, j, 3) +
      theme(legend.position = 'none',
            axis.title.y = element_text(margin=margin(r=7)),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank())
    L3est <- loads_plot_sim(Loads_est, Variances_est, j, 3) +
      guides(fill = guide_legend(nrow = 3, byrow=TRUE,
                                 override.aes = list(size = 0.1, 'mm'))) +
      theme(axis.line.x = element_blank(),
            legend.position = 'none') +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01),
                         breaks = c(-5e-16,0.0,5e-16))#+
    if(j==1){
      S1est <- S1aest
      S2est <- S2aest
      S3est <- S3aest
      
      S1est <- S1est + theme(
        legend.title = element_text(size=8,margin = margin(b = 1, unit = "pt")),
        legend.text = element_text(size=8, margin = margin(b = 1, unit = "pt")),
        legend.spacing = unit(0.01, 'cm'),
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.position = c(0.3, 0.88)) +
        guides(linetype = guide_legend(override.aes = list(size = 0.4)))
    }
    
    if (j==2){
      S1est <- S1est + theme(
        legend.title = element_text(size=8,margin = margin(b = 1, unit = "pt")),
        legend.text = element_text(size=8, margin = margin(b = 1, unit = "pt")),
        legend.spacing = unit(0.01, 'cm'),
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.position = c(0.3, 0.9)) +
        guides(linetype = guide_legend(override.aes = list(size = 0.4)))
      
      S3est_ <- S3est + 
        labs(x='time') +
        theme(axis.text.x = element_text(color="black"),
              axis.ticks.x = element_line(color="black"),
              axis.line.x = element_line(color="black"),
              axis.title.y = element_text(margin=margin(r=7)))
      
      L3est_ <- L3est +
        labs(x='metabolite') +
        theme(axis.ticks.x = element_line(color="black"), 
              axis.line.x = element_line(color="black"),
              axis.text.x = element_text(color='black'),
              legend.position = 'none')
    } else {
      S3est_ <- S3est + 
        labs(x='') +
        theme(axis.ticks.x = element_blank(), 
              axis.line.x = element_blank(),
              axis.text.x = element_text(color='black'),
              legend.position = 'none') +
        scale_x_continuous(labels = function(breaks) {rep_along(breaks, "")})

      L3est_ <- L3est +
        labs(x='') +
        theme(axis.ticks.x = element_blank(), 
              axis.line.x = element_blank(),
              axis.text.x = element_text(color='black'),
              legend.position = 'none') +
        scale_x_discrete(labels = rep("", 25), breaks = 1:25)
      
    }
    
    
    g <-
      cowplot::plot_grid(
        S1,
        L1,
        NULL,
        S1est,
        L1est,
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        S2_,
        L2_,
        NULL,
        S2est,
        L2est,
        ncol = 5,
        nrow = 3 ,
        rel_heights = c(1,
                        -0.1,
                        1),
        rel_widths = c(0.5,1,0.1,0.5,1),
        align = 'hv',
        axis = 'lb'
      )
    g + theme(plot.margin=margin(0.7,0.7,0.7,0.7, unit="cm"))
    
    h <-
      cowplot::plot_grid(
        S1,
        NULL,
        L1,
        NULL,
        NULL,
        NULL,
        S2 +
          labs(x='time, min') +
          theme(axis.text.x = element_text(color="black"),
                axis.ticks.x = element_line(color="black"),
                axis.line.x = element_line(color="black"),
                axis.title.y = element_text(margin=margin(r=7))),
        NULL,
        L2 +
          labs(x='metabolite') +
          theme(axis.ticks.x = element_line(color="black"), 
                axis.line.x = element_line(color="black"),
                axis.text.x = element_text(color='black'),
                legend.position = 'none'),
        NULL,
        NULL,
        NULL,
        S1est,
        NULL,
        L1est,
        NULL,
        NULL,
        NULL,
        S2est,
        NULL,
        L2est,
        NULL,
        NULL,
        NULL,
        S3est+ 
          labs(x='time') +
          theme(axis.text.x = element_text(color="black"),
                axis.ticks.x = element_line(color="black"),
                axis.line.x = element_line(color="black"),
                axis.title.y = element_text(margin=margin(r=7))),
        NULL,
        L3est+
          labs(x='metabolite') +
          theme(axis.ticks.x = element_line(color="black"), 
                axis.line.x = element_line(color="black"),
                axis.text.x = element_text(color='black'),
                legend.position = 'none'),
        ncol = 3,
        nrow = 9,
        rel_heights = c(1,
                        -0.1,
                        1,
                        0.2,
                        1,
                        -0.1,
                        1,
                        -0.1,
                        1),
        rel_widths = c(0.5,-0.04,1),
        align = 'hv',
        axis = 'lb',
        labels = c(panel[[j]][1],"","",
                   "","","",
                   "","","",
                   "","","",
                   panel[[j]][2],"","",
                   "","","",
                   "","","",
                   "","","",
                   "","",""),
        label_y=1.1
      )
    h + theme(plot.margin=margin(0.7,0.7,0.7,0.7, unit="cm"))
    
    list_of_plots_vert[[j]] <- g
    list_of_plots_horz[[j]] <- h
  }
  p2 <-
    cowplot::plot_grid(
      title_row,
      NULL,
      list_of_plots_vert[[1]],
      NULL,
      list_of_plots_vert[[2]],
      ncol = 1,
      nrow = 5 ,
      rel_heights = c(0.3,-0.1,1,0.1,1),
      align = 'h',
      axis = 'b')
  p2 + theme(plot.margin=margin(0.7,0.7,0.7,0.7, unit="cm"))
  

  p1 <-
    cowplot::plot_grid(
      title_row,
      NULL,
      title_row,
      NULL,
      NULL,
      NULL,
      list_of_plots_horz[[1]],
      NULL,
      list_of_plots_horz[[2]],
      ncol = 3,
      nrow = 3 ,
      rel_widths = c(1,0.1,1),
      rel_heights=c(0.3,-0.1,1),
      align = 'h',
      axis = 'b',
      labels = c("","","",
                 "","","",
                 "","","")
    )
  p1 + theme(plot.margin=margin(0.7,0.7,0.7,0.7, unit="cm"))
  
  ggsave(filename=paste0(Sys.Date(),"_sim_fixed",
                         ".png"),
         #fallback_resolution=600,
         plot=p1,
         width = 9,
         height= 9)
  
  
}

fig_enc_est <- function(){
  
  # encoded vs estimated loadings
  
  Loads_estA <- bootstrap_results$Pa %>%
    as.data.frame() %>%
    mutate(metab = as.numeric(var),
           var=NULL) %>%
    filter(PC<3)%>%
    group_by(PC, metab) %>%
    summarise(loadi = quantile(loading, c(0.025, 0.5, 0.975)), prob = c(0.025,0.5,0.975)
    ) %>%
    ungroup() %>%
    mutate(eff = "Time") %>%
    pivot_wider(id_cols = c(PC, metab, eff), values_from = loadi, names_from = prob)
  
  
  Loads_estAB <- bootstrap_results$Pab %>%
    as.data.frame() %>%
    mutate(metab = as.numeric(var)) %>%
    filter(PC<3)%>%
    group_by(PC, metab) %>%
    summarise(
      loadi = quantile(loading, c(0.025, 0.5, 0.975)), prob = c(0.025,0.5,0.975)
    ) %>%
    ungroup() %>%
    mutate(eff = "Time-Treatment") %>%
    pivot_wider(id_cols = c(PC, metab, eff), values_from = loadi, names_from = prob)
  
  
  Loads_est <- rbind(Loads_estA, Loads_estAB)
  
  
  LL <- do.call("cbind", loads) %>%
    as_tibble() %>%
    mutate(V3=NULL,
           V6=NULL) %>%
    rename(`Time_1`=V1,`Time_2`=V2,`Time-Treatment_1`=V4,`Time-Treatment_2`=V5) %>%
    rowid_to_column("metab") %>%
    pivot_longer(cols=!metab, names_to="eff") %>%
    separate(eff,
             into=c('eff','PC'),
             sep ="_")
  
  loa <- full_join(LL, Loads_est, by=c("metab", "PC", "eff")) %>%
    mutate(PC = paste0("PC", PC))
  
  Lp <- ggplot(loa) +
    geom_errorbar(aes(x=value,ymin=`0.025`,ymax=`0.975`)) +
    geom_point(aes(x=value,y=`0.5`)) +
    geom_abline() +
    theme(aspect.ratio = 1) +
    facet_wrap(eff~PC, scales = "free", nrow = 2)+
    theme(axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color="black"),
          axis.ticks = element_line(color="black"),
          panel.border = element_rect(color="black", fill = NA),
          strip.background = element_rect(color = "black")) +
    xlab("Encoded Loading") +
    ylab("Estimated Loading")
  
  ggsave(filename=paste0(Sys.Date(),
                         "_enc_est_loadings.png"),
         plot=Lp,
         width = 6, 
         height= 6)
  
  # encoded vs estimated scores 
  
  scores_estA <- bootstrap_results$Ta %>%
    as.data.frame() %>%
    filter(PC!=3)%>%
    group_by(PC, treatment, tp_L) %>%
    summarise(
      score = quantile(score, c(0.025, 0.5, 0.975)), prob = c(0.025,0.5,0.975)
      
    ) %>%
    ungroup() %>%
    mutate(eff="Time") %>%
    pivot_wider(id_cols = c(PC, treatment, tp_L, eff), values_from = score, names_from = prob)
  
  
  scores_estAB <- bootstrap_results$Tab %>%
    as.data.frame() %>%
    filter(PC!=3)%>%
    group_by(PC, treatment, tp_L) %>%
    summarise(
      score = quantile(score, c(0.025, 0.5, 0.975)), prob = c(0.025,0.5,0.975)
      
    ) %>%
    ungroup() %>%
    mutate(eff="Time-Treatment") %>%
    pivot_wider(id_cols = c(PC, treatment, tp_L, eff), values_from = score, names_from = prob)
  
  scores_est <- rbind(scores_estA, scores_estAB)
  
  SS <- do.call("cbind", scores) %>%
    as_tibble() %>%
    mutate(V3=NULL,
           V6=NULL) %>%
    rename(`Time_1`=V1,`Time_2`=V2,`Time-Treatment_1`=V4,`Time-Treatment_2`=V5) %>%
    head(22) %>%
    mutate(tp_L = c(seq(0,10),seq(0,10)),
           treatment = as.factor(c(rep(-1,11),rep(1,11)))) %>%
    pivot_longer(cols=!c(tp_L, treatment), names_to="eff") %>%
    separate(eff,
             into=c('eff','PC'),
             sep ="_") %>%
    rename(score=value)
  
  
  sco <- full_join(SS, scores_est, by=c("tp_L", "treatment", "eff", "PC")) %>%
    mutate(PC = paste0("PC", PC))
  
  Sp <- ggplot(sco) +
    geom_errorbar(aes(x=score,ymin=`0.025`,ymax=`0.975`)) +
    geom_point(aes(x=score,`0.5`)) +
    geom_abline() +
    theme(aspect.ratio = 1) +
    facet_wrap(eff~PC, scales = "free", nrow=2) +
    theme(axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color="black"),
          axis.ticks = element_line(color="black"),
          panel.border = element_rect(color="black", fill = NA),
          strip.background = element_rect(color = "black")) +
    xlab("Encoded Score") +
    ylab("Estimated Score")
  
  ggsave(filename=paste0(Sys.Date(),
                         "_enc_est_scores.png"),
         plot=Sp,
         width = 6, 
         height= 6)
  
  
}
