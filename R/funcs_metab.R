
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



data_design <- function(n_treatments,
                        n_timepoints,
                        timepoints,
                        treatments,
                        n_PC,
                        n_variables,
                        n_subjects) {
  # list of design properties of data
  
  output <- list()
  output$n_timepoints <- n_timepoints
  output$timepoints <- timepoints
  output$n_treatments <- n_treatments
  output$treatments <- treatments
  output$n_PC <- n_PC
  output$n_variables <- n_variables
  output$variables <- n_variables
  output$n_subjects <- n_subjects
  output$n_observations <- prod(n_timepoints,
                                n_treatments,
                                n_subjects)
  return(output)
}

# Create Design matrix 
gen_design <- function(design_info_list){
  ## functional forms for time model
  
  t <- seq(1, 11, 1)
  spline <- ns(t,2)


  output_df <- tibble(
    # column of subjects
    subject = 
      seq(design_info_list$n_subjects) %>% 
      rep(design_info_list$n_treatments) %>%
      rep(design_info_list$n_timepoints) %>% 
      factor(),
    
    # time keeping
    t = t %>%
      rep(each=design_info_list$n_treatments*
            design_info_list$n_subjects),
    # spline
    tp_L = spline[,1] %>%
      rep(each=design_info_list$n_treatments*
            design_info_list$n_subjects),
    tp_Q = spline[,2] %>%
      rep(each=design_info_list$n_treatments*
            design_info_list$n_subjects),
    # column of treatments
    treatment = c('A','B') %>%
      rep(each=design_info_list$n_subjects) %>%
      rep(design_info_list$n_timepoints) %>%
      factor())
  # set contrasts
  contrasts(output_df$treatment) <- contr.sum(design_info_list$n_treatments)
  return(output_df)}


fit_lmm <- function(long_response_df) {
  # fit LMM
  
  cl <- parallel::makeCluster(12)
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
            subject),
        data = x,
        na.action = na.exclude,
        REML = TRUE
      )))
  
  parallel::stopCluster(cl)
  future:::ClusterRegistry("stop")
  plan(sequential)
  return(res)
}

extract_coefs <- function(model_df,
                           include_residuals = FALSE){
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
  USV <- svd(input_matrix)
  output <- list()
  output$scores <- (USV$u*sqrt(dim(input_matrix)[1]-1))[,1:10]
  output$loadings <- (USV$v %*% diag(USV$d) / sqrt(dim(input_matrix)[1]-1))[,1:10]
  output$singular_values = USV$d[1:10]
  return(output)
}

pca_effect_matrices <- function(input_matrix) {
  return(
    lapply(input_matrix,
           FUN = "col_center") %>%
      lapply(FUN = "perform_pca")
  )
}

expl_var_from_svd <- function(singular_values,
                              decimals = 1){
  return(round(singular_values^2/sum(singular_values^2)*100,decimals))
}

col_center <- function(data_matrix){
  return(sweep(x      = data_matrix,
               MARGIN = 2,
               STATS  = colMeans(data_matrix),
               FUN    = "-"))}

# scale data by SD at baseline
SDscale <- function(data){
  data <- data %>%
    group_by(metabolite) %>%
    mutate(concentration=concentration_raw/sd(concentration_raw[Time==1],na.rm = T))
  return(data)
}

orthprocr_pca <- function(target,query){
  sol <- query
  sol$loadings <- cds::orthprocr(Z = target$loadings,
                                 X = query$loadings)$XQ
  sol$scores <- sol$scores %*% (cds::orthprocr(Z = target$loadings,
                                               X = query$loadings)$Q %>% solve())
  return(sol)
}

fit_bootstrap_models <- function(fitted_models,
                                 design_info_list){

  original_data <- do.call("rbind",fitted_models$data) %>%
    mutate(metabolite = rep(fitted_models$metabolite, each = nrow(fitted_models$data[[1]]))) %>%
    arrange(metabolite, subject, treatment, t)

  bootstrap_models <- list()
  bootstrap_metabolites  <- list()
  bootstrap_PCA <- list()
  bootstrap_res <- list()

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
      group_by(subject,treatment,t) %>%
      mutate(dupeID = row_number()) %>%
      ungroup() %>%
      group_by(metabolite) %>%
      arrange(dupeID,subject,treatment,t) %>%
      mutate(concentration=concentration_raw/sd(concentration_raw[t==1],na.rm = T)) %>% # scaling
      group_nest() %>%
      mutate(sel = fitted_models$sel) %>%
      fit_lmm(long_response_df = .) %>%
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
  return(output)
  }


process_bootstrap_models <- function(bootstrap_models_list,
                                     design_info_list,
                                     effects){
  shared_fitted_metabolites <- bootstrap_models_list$metabolites
  bootstrap_models <- bootstrap_models_list$models
  fitted_models <- bootstrap_models_list$full_model 
  bootstraps <- bootstrap_models_list$bootstraps
  bootstrap_res <- list()

  effect_matrices <- compose_effmat(fitted_models,design_info_list$design)
 
  reference_pca <- effect_matrices %>%
    pca_effect_matrices()

  for(i in 1:length(bootstrap_models)){

    bootstrap_fold_design <- design_info_list$design

    fold_coef_mat <-compose_effmat(bootstrap_models[[i]],design_info_list$design)
    
    fold_pca_res <- fold_coef_mat %>%
      pca_effect_matrices()
    
    
    fold_pca_res$Ma <- orthprocr_pca(reference_pca$Ma,
                                     fold_pca_res$Ma)
    fold_pca_res$Mab <- orthprocr_pca(reference_pca$Mab,
                                      fold_pca_res$Mab)
    fold_pca_res$Mfull <- orthprocr_pca(reference_pca$Mfull,
                                        fold_pca_res$Mfull)


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
      cbind(bootstrap_fold_design) %>%
      pivot_longer(cols=-colnames(bootstrap_fold_design),
                   names_to = "PC",
                   values_to = "score",
                   names_prefix = "V") %>%
      mutate(PC = as.factor(PC),
             treatment = as.factor(treatment)) %>%
      filter(PC %in% seq(1:3))

    fold_Tab <- fold_pca_res$Mab$scores %>%
      as_tibble(.name_repair)%>%
      cbind(bootstrap_fold_design) %>%
      pivot_longer(cols=-colnames(bootstrap_fold_design),
                   names_to = "PC",
                   values_to = "score",
                   names_prefix = "V")%>%
      mutate(PC = as.factor(PC),
             treatment = as.factor(treatment)) %>%
      filter(PC %in% seq(1:3))
   
    fold_Tfull <- fold_pca_res$Mfull$scores %>%
      as_tibble(.name_repair)%>%
      cbind(bootstrap_fold_design) %>%
      pivot_longer(cols=-colnames(bootstrap_fold_design),
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
    output$Mfull_ <- do.call("cbind",output_$Mfull)
    if (include_residuals) {
      output$Mres <-
        do.call("cbind", ref_mod$resids)
    }
    return(output)
  }

score_plot <- function(Scores, Variances,label,j,i){
  S <- Scores[[j]] %>%
    filter(PC == i) %>%
    group_by(t) %>%
    mutate(
      quantile_025 = quantile(score, .025),
      quantile_500 = quantile(score, .500),
      quantile_975 = quantile(score, .975)
    ) %>%
    ungroup() %>%
    ggplot(aes(x = t * 30,
               y = quantile_500)) +
    
    geom_hline(yintercept = 0, linetype='dotted', color="grey") +
    geom_ribbon(aes(ymin = quantile_025, 
                    ymax = quantile_975),
                alpha=0.3
    ) +
    
    geom_line(aes(x = t * 30,
                  y = quantile_500)) +
    
    scale_linetype_manual(values = c('solid','dashed'), labels = c("A", "B")) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    theme_classic() +
    theme(
      plot.margin = margin(0, 0, 0, 0, unit = "cm"),
      axis.text.x = element_blank(),
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

score_plot_treat <- function(Scores, Variances,label,j,i){
  S <- Scores[[j]] %>%
    filter(PC == i) %>%
    mutate(treatment = as.factor(treatment)) %>%
    group_by(t, treatment) %>%
    mutate(
      quantile_025 = quantile(score, .025),
      quantile_500 = quantile(score, .500),
      quantile_975 = quantile(score, .975)
    ) %>%
    ungroup() %>%
    ggplot(aes(x = t * 30,
               y = quantile_500)) +
    
    geom_hline(yintercept = 0, linetype='dotted', color="grey") +
    geom_ribbon(aes(ymin = quantile_025, 
                    ymax = quantile_975,
                    fill = treatment),
                alpha=0.5
                ) +
    
    geom_line(aes(x = t * 30,
                  y = quantile_500,
                  linetype = treatment)) +

    scale_linetype_manual(values = c('solid','dashed'), labels = c("A", "B")) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    theme_classic() +
    theme(
      plot.margin = margin(0, 0, 0, 0, unit = "cm"),
      axis.text.x = element_blank(),
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


loads_plot <- function(Loads, Variances, j, i){
  
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
    "S_LDL_FC",
    "S_LDL_PL",
    "S_LDL_TG"
  )

    
  cols_ <- c(
    "#0c2c84",
    "#225ea8",
    "#1d91c0",
    "#41b6c4",
    "#7fcdbb",
    "#c7e9b4",
    "#ae017e",
    "#f768a1",
    "#fbb4b9"
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
                 fill=interaction(size, lipoclass)))+ 

    geom_bar(stat="identity",
             position="dodge"
    )+ 
    scale_fill_manual(values = cols_)+
    guides(fill = guide_legend(override.aes = list(size = 0.01),
                                nrow=6,ncol = 2))+
    geom_errorbar(aes(ymin=quantile_025,
                      ymax=quantile_975),
                  size=0.4,
                  stat = 'identity',
                  position = 'dodge')+

    geom_hline(yintercept = 0, linetype='dotted', color="grey") +
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
          legend.title = element_blank(),
          legend.text = element_text(color='black', size=3),
          legend.spacing = unit(0.01, 'cm'),
          legend.key.size = unit(0.2, "cm"),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA),
          strip.background = element_blank(),
          strip.text.x = element_blank()) + 
    scale_x_discrete(limits = pos, labels = rep(c("CE","FC","PL","TG"),9)) 
}

fig_est_scores_loadings_JK_re2 <- function(titles, Scores, Loads, Variances){
  list_of_plots_horz <- list()
  for (j in 1:length(Scores)){
    n_PC <- 3
    labels_ <- list()
    explvars <- list()
    for (i in 1:n_PC){
      labels_[[i]] <- paste0(
        "PC",
        i,
        "\n(",
        filter(Variances[[j]],
               PC == i)$perc_expl %>%
          quantile(.025) %>%
          round(1) %>%
          format(nsmall = 1),
        "-",
        filter(Variances[[j]],
               PC == i)$perc_expl %>%
          quantile(.975) %>%
          round(1) %>%
          format(nsmall = 1),
        "%)")
      
      explvars[[i]] <- paste0(
        filter(Variances[[j]],
               PC == i)$perc_expl %>%
          quantile(.025) %>%
          round(1) %>%
          format(nsmall = 1),
        "-",
        filter(Variances[[j]],
               PC == i)$perc_expl %>%
          quantile(.975) %>%
          round(1) %>%
          format(nsmall = 1),
        "%)"
      ) 
      
      titles_ <- if(i==1) {
        c("Scores", "Loadings")
      } else{
        c("", "")
      }
    }
    
    labelss <- unlist(labels_, recursive=FALSE)
    
    
    S1a <- score_plot(Scores, Variances, labelss, j, 1) +
      guides(linetype = guide_legend(nrow = 1)) +
      theme(
        legend.position = "none",
        axis.title.y = element_text(margin=margin(r=5)),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank())
    
    S1 <- score_plot_treat(Scores, Variances, labelss, j, 1) +
      guides(linetype = guide_legend(nrow = 1)) +
      theme(
        legend.position = "none",
        axis.title.y = element_text(margin=margin(r=5)),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank())
    
    L1 <- loads_plot(Loads, Variances, j, 1) +
      theme(
        axis.line.x = element_blank(),
        legend.position = "none")
    
    S2a <- score_plot(Scores, Variances, labelss, j, 2) +
      theme(legend.position = 'none',
            axis.title.y = element_text(margin=margin(r=5)),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank())
    
    S2 <- score_plot_treat(Scores, Variances, labelss, j, 2) +
      theme(legend.position = 'none',
            axis.title.y = element_text(margin=margin(r=5)),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank())
    
    if(j==1){
      L1 <- L1 + theme(legend.title = element_blank(),
                       legend.text = element_text(color='black', 
                                                  size=6,
                                                  margin=margin(t=1,r=1,b=1,l=1, unit = "pt")),
                       legend.margin=margin(c(0,0,0,0)),
                       legend.spacing.x = unit(0, "mm"),
                       legend.spacing.y = unit(0, "mm"),
                       legend.key.size = unit(0.25, "cm"),
                       legend.position = c(0.8, 0.8)) +
        guides(fill = guide_legend(override.aes = list(size = 0.01),
                                                nrow=6,ncol = 2))
      
      S1 <- S1 + theme(
        legend.title = element_text(size=8,margin = margin(b = 1, unit = "pt")),
        legend.text = element_text(size=8, margin = margin(b = 1, unit = "pt")),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.position = c(0.3, 0.95)) +
        guides(linetype = guide_legend(byrow=TRUE,nrow = 2,override.aes = list(size = 0.1)))
    }
    L2 <- loads_plot(Loads, Variances, j, 2) +
      guides(fill = guide_legend(nrow = 3, byrow=TRUE,
                                 override.aes = list(size = 0.1, 'mm'))) +
      theme(axis.line.x = element_blank(),
            legend.position = 'none')
    
    S3a <- score_plot(Scores, Variances, labelss, j, 3) +
      theme(legend.position = 'none',
            axis.title.y = element_text(margin=margin(r=5)),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank())
    
    S3 <- score_plot_treat(Scores, Variances, labelss, j, 3) +
      theme(legend.position = 'none',
            axis.title.y = element_text(margin=margin(r=5)),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank())
    
    L3 <- loads_plot(Loads, Variances, j, 3) +
      guides(fill = guide_legend(nrow = 3, byrow=TRUE,
                                 override.aes = list(size = 0.1, 'mm'))) +
      theme(axis.line.x = element_blank(),
            legend.position = 'none')
    
    if (j==2) {
      S1 <- S1a
      S2 <- S2a
      S3 <- S3a
      
      L1 <- L1 +
        theme(legend.title = element_blank(),
              legend.text = element_text(color='black', 
                                         size=6,
                                         margin=margin(t=1,r=1,b=1,l=1, unit = "pt")),
              legend.margin=margin(c(0,0,0,0)),
              legend.spacing.x = unit(0, "mm"),
              legend.spacing.y = unit(0, "mm"),
              legend.key.size = unit(0.25, "cm"),
              legend.position = c(0.8, 0.8)) +
        guides(fill = guide_legend(override.aes = list(size = 0.01),
                                   nrow=6,ncol = 2))
    }
    if (j==3){
      S1 <- S1 + 
        theme(legend.title = element_text(size=8,margin = margin(b = 1, unit = "pt")),
              legend.text = element_text(size=8, margin = margin(b = 1, unit = "pt")),
              legend.spacing.x = unit(0, "mm"),
              legend.spacing.y = unit(0, "mm"),
              legend.key.size = unit(0.2, "cm"),
              legend.key.width = unit(0.2, "cm"),
              legend.position = c(0.33, 0.95))+
        guides(linetype = guide_legend(byrow=TRUE,nrow = 2,override.aes = list(size = 0.1)))
      
      S3_ <- S3 + 
        labs(x='time, min') +
        theme(axis.text.x = element_text(color="black"),
              axis.ticks.x = element_line(color="black"),
              axis.line.x = element_line(color="black"),
              axis.title.y = element_text(margin=margin(r=5)))
      
      L3_ <- L3 +
        labs(x='metabolite') +
        theme(axis.ticks.x = element_line(color="black"), 
              axis.line.x = element_line(color="black"),
              axis.text.x = element_text(color='black', size=6, angle = 90),
              legend.title = element_blank())
    } else {
      S3_ <- S3 + 
        labs(x='') +
        theme(axis.ticks.x = element_blank(), 
              axis.line.x = element_blank(),
              axis.text.x = element_text(color='black', size=5),
              legend.position = 'none') +
        scale_x_continuous(labels = function(breaks) {rep_along(breaks, "")})

      L3_ <- L3 +
        labs(x='') +
        theme(axis.ticks.x = element_blank(), 
              axis.line.x = element_blank(),
              axis.text.x = element_text(color='black', size=6, angle = 90),
              legend.position = 'none') +
        scale_x_discrete(labels = rep("", 25), breaks = 1:25)

    }
    
    Stitle <- ggdraw() +
      draw_label(
        "Scores",
        x = 0.5,
        hjust = 0.5
      ) +
      theme(
        plot.margin = margin(0, 0, 0, 0)
      )
    
    Ltitle <- ggdraw() +
      draw_label(
        "Loadings",
        x = 0.5,
        hjust = 0.5
      ) +
      theme(
        plot.margin = margin(0, 0, 0, 0)
      )
    
    title_row <- cowplot::plot_grid(
      Stitle,
      Ltitle,
      ncol = 2,
      nrow = 1 ,
      align = 'hv',
      axis = 'lb'
    )
    
    h <-
      cowplot::plot_grid(
        S1,
        NULL,
        L1,
        NULL,
        NULL,
        NULL,
        S2,
        NULL,
        L2,
        NULL,
        NULL,
        NULL,
        S3+ 
          labs(x='time, min') +
          theme(axis.text.x = element_text(color="black"),
                axis.ticks.x = element_line(color="black"),
                axis.line.x = element_line(color="black"),
                axis.title.y = element_text(margin=margin(r=5))),
        NULL,
        L3+
          labs(x='metabolite') +
          theme(axis.ticks.x = element_line(color="black"), 
                axis.line.x = element_line(color="black"),
                axis.text.x = element_text(color='black', size=6, angle = 90)),
        ncol = 3,
        nrow = 5,
        rel_heights = c(1,
                        -0.1,
                        1,
                        -0.1,
                        1),
        rel_widths = c(0.5,-0.1,1),
        align = 'hv',
        axis = 'lb'
      )
    
    h + theme(plot.margin=margin(0.7,0.7,0.7,0.7, unit="cm"))

    list_of_plots_horz[[j]] <- h
  }
  
  
  p_ <- cowplot::plot_grid(
    title_row,
    NULL,
    list_of_plots_horz[[1]],
    ncol = 1,
    nrow = 3,
    rel_heights = c(0.3,-0.05,1),
    align = 'hv',
    axis = 'lb'
  )
  p_ + theme(plot.margin=margin(0.7,0.7,0.7,0.7, unit="cm"))
  
  ggsave(filename=paste0(Sys.Date(),"melc_full_fixed",
                         ".png"),
         #fallback_resolution=600,
         plot=p_,
         width = 4.5,
         height= 4.5)
  
  p2 <-
    cowplot::plot_grid(
      title_row,
      NULL,
      title_row,
      NULL,
      NULL,
      NULL,
      list_of_plots_horz[[2]],
      NULL,
      list_of_plots_horz[[3]] + theme(legend.position = 'none'),
      ncol = 3,
      nrow = 3 ,
      rel_widths = c(1,0.1,1),
      rel_heights=c(0.3,-0.05,1),
      align = 'hv',
      axis = 'lb',
      labels = c("a"," ","b" )
    )
  p2 + theme(plot.margin=margin(0.7,0.7,0.7,0.7, unit="cm"))
  
  ggsave(filename=paste0(Sys.Date(),"_melc_time_treat_fixed",
                         ".png"),
         #fallback_resolution=600,
         plot=p2,
         width = 9,
         height= 4.52)
  
  return(list_of_plots_horz)
}

create_effect_matrices <- function(design_matrix, sim_models, coeffs_, i) {
  
  des <- list()
  
  X <- design_matrix %>%
    filter(subject==1) %>%
    model.matrix(
      ~ tp_L +
        tp_Q +
        tp_L:treatment +
        tp_Q:treatment
      ,
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

fig_random_effects <- function(sim_models, design_matrix) {
  # create random effect effect matrices, then do bar plot, biplot
  
  n_ind <- length(unique(sim_models$M[[1]]@frame$subject))
  
  time_timetreat <- list()
  timetreat <- list()
  time <- list()
  for (i in 1:36) {
    des_t <-
      create_effect_matrices(design_matrix, sim_models, all_of(c("tp_L", "tp_Q")), i)
    time[[paste0('metab_', i)]] <- des_t$X %*% des_t$gamma
    
    des_ttr <-
      create_effect_matrices(design_matrix, sim_models, all_of(c("tp_L:treatment1", "tp_Q:treatment1")), i)
    timetreat[[paste0('metab_', i)]] <- des_ttr$X %*% des_ttr$gamma
    
    des_t_ttr <-
      create_effect_matrices(design_matrix, sim_models, all_of(c("tp_L", "tp_Q", "tp_L:treatment1", "tp_Q:treatment1")), i)
    time_timetreat[[paste0('metab_', i)]] <-
      des_t_ttr$X %*% des_t_ttr$gamma
  }
  
  timedf_t <- do.call(cbind, lapply(time, function(x) Reshape(x,440,1)))
  timetreatdf_t <- do.call(cbind, lapply(timetreat, function(x) Reshape(x,440,1)))
  time_timetreatdf_t <- do.call(cbind, lapply(time_timetreat, function(x) Reshape(x,440,1)))
  
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
    as.matrix(design_matrix[,c("subject","t","treatment")]),
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
  
  eff <- c("Time", "Time-Treatment")
  
  pp <- list()
  pp_ <- list()
  for (i in 1:2) {
    df <- subset(scoresDF_t_wide, Effect == eff[i])
    
    pp[[i]] <- ggplot(df, aes(x = time*30, y = Score_PC1, group=interaction(individual,treatment), linetype=treatment)) +
      geom_line() +
      labs(
        x = "",
        y = paste0("PC1", " ", df$var_PC1[[1]])#,
      ) +
      theme_bw() +
      theme(
        plot.margin = margin(0, 0, 0, 0, unit = "cm"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
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
    
    pp_[[i]] <- ggplot(df, aes(x = time*30, y = Score_PC2, group=interaction(individual,treatment), linetype=treatment)) +
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
    pp[[2]]+theme(legend.position = c(0.2,0.1),
                   legend.text = element_text(color='black', 
                                              size=6,
                                              margin=margin(t=1,r=1,b=1,l=1, unit = "pt")),
                   legend.margin=margin(c(0,0,0,0)),
                   legend.spacing.x = unit(0, "mm"),
                   legend.spacing.y = unit(0, "mm"),
                   legend.key.size = unit(0.25, "cm"),
                   legend.title=element_text(size=8)),

    NULL,NULL,NULL,
    pp_[[1]],
    NULL,
    pp_[[2]],
    ncol = 3,
    nrow = 3,
    rel_widths = c(1,0.1,1),
    rel_heights = c(1,-0.05,1),
    labels = c("a","","b"),
    align = 'hv',
    axis = 'lb'
  )
  
  ggsave(filename=paste0(Sys.Date(),"_melc_random_effects_ind_t",
                         ".png"),
         #fallback_resolution=600,
         plot=pl_re,
         width = 5.33,
         height= 5.1)

  
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
    filename = paste0(Sys.Date(), "melc_RE_scores.png"),
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
    rel_widths = c(1,0.1,1),
    labels = c("a","","b"),
    align = 'hv',
    axis = 'lb'
  )
  
  ggsave(filename=paste0(Sys.Date(),"melc_random_effects_ind",
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
  
  ggsave(filename=paste0(Sys.Date(),"_melc_re_highlighted",
                         ".png"),
         #fallback_resolution=600,
         plot=p,
         width = 4,
         height= 4)
  
  
}

fig_data <- function(data){
  
  p1 = ggplot(data=data, aes(x = Time*30, 
                                     y = concentration_raw, 
                                     linetype=treatment,
                                     color=treatment,
                                     group = interaction(subject, treatment))) +
    geom_line() +
    theme_bw() +
    facet_wrap(~metabolite, nrow = 9, ncol=4, 
               strip.position ="top", scales = "free_y") +
    labs(x = "time, min", y= "concentration, mmol/L") +
   
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
      strip.text = element_text(size=8)) +
    guides(linetype = guide_legend(nrow=1))

  
  ggsave(paste0(Sys.Date(),"_melc_data.png"),
         p1, width = 6.5, height = 8)
  
  p2 = ggplot(data=data, aes(x = Time*30, 
                            y = concentration, 
                            linetype=treatment,
                            group = interaction(subject, treatment))) +
    geom_line(color="gray") +
    
    geom_line(data=subset(lipo_dat, subject %in% c(7,14,17)),aes(color=subject))+
    theme_bw() +
    facet_wrap(~metabolite, nrow = 5, ncol=8, 
               strip.position ="top", as.table = F) +
    labs(x = "time, min", y= " standardized concentration") +
    theme(
      plot.margin = unit(c(0.06,0.06,0.06,0.06), "cm"),
      axis.ticks = element_line(color="black"),
      axis.text = element_text(color="black"),
      axis.line = element_line(color="black"),
      legend.background = element_rect(color = NA, fill = "transparent"),
      legend.box.background = element_rect(color = NA, fill = "transparent"),
      legend.key = element_rect(color = NA, fill = "transparent"),
      legend.position = c(0.74,0.91),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
      strip.background = element_blank(),
      strip.text = element_text(size=8),
      legend.direction = "vertical", 
      legend.box = "horizontal")
  show(p2)
  
  ggsave(paste0(Sys.Date(),"_melc_data_highlighted.png"),
         p2, width = 7, height = 7)
  
  data_tg <- lipo_dat %>%
    filter(grepl("TG",metabolite))
  
  p3 = ggplot(data=data_tg, aes(x = Time*30, 
                             y = concentration_raw, 
                             linetype=treatment,
                             color=treatment,
                             group = interaction(subject, treatment))) +
    geom_line() +
    
    theme_bw() +
    facet_wrap(~metabolite, scales = "free_y", nrow = 2, ncol=5, 
               strip.position ="top") +
    labs(x = "time, min", y= "concentration, mmol/L") +
    
    theme(
      plot.margin = unit(c(0.06,0.06,0.06,0.06), "cm"),
      axis.ticks = element_line(color="black"),
      axis.text = element_text(color="black"),
      axis.line = element_line(color="black"),
      legend.background = element_rect(color = NA, fill = "transparent"),
      legend.box.background = element_rect(color = NA, fill = "transparent"),
      legend.key = element_rect(color = NA, fill = "transparent"),
      legend.position = c(0.91,0.15),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
      strip.background = element_blank(),
      strip.text = element_text(size=8)) +
    guides(linetype = guide_legend(nrow=1))
  
  ggsave(paste0(Sys.Date(),"_melc_data_TG.png"),
         p3, width = 7, height = 2.5)
  
  
}
