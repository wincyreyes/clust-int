.libPaths("/home/wincy.reyes@IUF.LAN/R/x86_64-pc-linux-gnu-library/4.3")
# install.packages("purrr")
library("scrime")
library("MASS")
library("glmnet")
library("GMCM")
library("scales")
library("dplyr")
library("tidyr")
library("purrr")
library(mclust)
library("ks")

library("ggplot2")
prepareSNPData <- function(snpData, fill = 0) {
  x <- snpData["x"]$x
  y <- snpData["y"]$y
  fill_matrix <- matrix(NA, nrow = length(y), ncol = fill)
  preparedSNPData <- cbind(y = y, x, fill_matrix)
  mode(preparedSNPData) <- "integer"
  return(preparedSNPData)
}
set.seed(12345)
prepareENVData <- function(N_samp,p_env, mu_env = 10,corr_env,var_env=10){
  mu_vec <- rep(mu_env, p_env)                  
  Sigma_env <- matrix(corr_env*10, p_env, p_env)   # Off-diagonal correlation 
  diag(Sigma_env) <- var_env                     # Diagonal variance >=1
  ENVdata<-  mvrnorm(N_samp,mu_vec,Sigma_env)#rlnorm.rplus(N_samp, meanlog = mu_vec, varlog = Sigma_env)
  colnames(ENVdata) <- paste("ENV", 1:p_env, sep="")
  # ENVdata <- scale(ENVdata)
  return(ENVdata)
}

elnetweights <- function(foldid, Xmat, y, alpha=1, family="binomial"){
  as.numeric(coef(cv.glmnet(x=data.matrix(Xmat),y=y,foldid=foldid,alpha=alpha, standardize=FALSE,family="binomial"), s = "lambda.min"))[-1]
}
riskscore <- function(weights, Xmat){
  Xmat%*%weights
}
riskscores_sub_abs <- function(subset, weights.all, Xmat){
  if(sum(weights.all[subset])!=0){
    weights.all.norm <- weights.all[subset]/sum(abs(weights.all[subset]))
  }else{ if(all(weights.all[subset]==0)){weights.all.norm <- rep(1/length(subset),length(subset))}}
  risk <- riskscore(weights.all.norm, as.matrix(Xmat[,subset,drop=FALSE]))
  return(list(risk=risk,weights.norm=weights.all.norm))
}


to_quartiles <- function(x) {
  if (length(unique(x)) < 4) {
    factor(rep("Medium-High", length(x)), 
           levels = c("Low", "Medium-Low", "Medium-High", "High"))
  } else {
    # Equal-count groups (quantiles)
    f <- cut_number(x, n = 4, labels = c("Low", "Medium-Low", "Medium-High", "High"))
    return(f)
  }
}

to_tertiles <- function(x) {
  if (length(unique(x)) < 3) {
    factor(rep("Medium", length(x)), 
           levels = c("Low", "Medium", "High"))
  } else {
    f <- cut_number(x, n = 3, labels = c("Low", "Medium", "High"))
    return(f)
  }
}

to_halves <- function(x) {
  if (length(unique(x)) < 2) {
    factor(rep("High", length(x)), 
           levels = c("Low", "High"))
  } else {
    f <- cut_number(x, n = 2, labels = c("Low", "High"))
    return(f)
  }
}

split_mid90 <- function(x) {
  # Ensure consistent factor levels for output
  out_levels <- c("Low", "Medium", "High")
  
  n <- length(x)
  if (n == 0) {
    return(factor(character(0), levels = out_levels))
  }
  
  # If fewer than 3 unique non-NA values, mimic original behavior and return all "Medium"
  if (length(unique(x[!is.na(x)])) < 3) {
    return(factor(rep("Medium", n), levels = out_levels))
  }
  
  # Compute 5th and 95th percentiles on non-NA values
  q <- stats::quantile(x, probs = c(0.05, 0.95), na.rm = TRUE, names = FALSE, type = 7)
  lower <- q[1]
  upper <- q[2]
  
  # If quantiles coincide (rare), fall back to terciles based on 1/3 and 2/3 probs
  if (is.na(lower) || is.na(upper) || lower == upper) {
    q2 <- stats::quantile(x, probs = c(1/3, 2/3), na.rm = TRUE, names = FALSE, type = 7)
    lower <- q2[1]
    upper <- q2[2]
    # If still equal (extremely degenerate), return all "Medium"
    if (is.na(lower) || is.na(upper) || lower == upper) {
      return(factor(ifelse(is.na(x), NA, "Medium"), levels = out_levels))
    }
  }
  
  # Cut into three bins: (-Inf, lower], (lower, upper), [upper, Inf)
  # include.lowest=TRUE ensures values equal to the global minimum are included
  res <- cut(x,
             breaks = c(-Inf, lower, upper, Inf),
             labels = out_levels,
             right = TRUE,
             include.lowest = TRUE)
  
  # Ensure factor has the desired levels even if some are empty
  res <- factor(as.character(res), levels = out_levels)
  
  return(res)
}

split_mid80 <- function(x) {
  # Ensure consistent factor levels for output
  out_levels <- c("Low", "Medium", "High")
  
  n <- length(x)
  if (n == 0) {
    return(factor(character(0), levels = out_levels))
  }
  
  # If fewer than 3 unique non-NA values, mimic original behavior and return all "Medium"
  if (length(unique(x[!is.na(x)])) < 3) {
    return(factor(rep("Medium", n), levels = out_levels))
  }
  
  # Compute 10th and 90th percentiles on non-NA values
  q <- stats::quantile(x, probs = c(0.1, 0.9), na.rm = TRUE, names = FALSE, type = 7)
  lower <- q[1]
  upper <- q[2]
  
  # If quantiles coincide (rare), fall back to terciles based on 1/3 and 2/3 probs
  if (is.na(lower) || is.na(upper) || lower == upper) {
    q2 <- stats::quantile(x, probs = c(1/3, 2/3), na.rm = TRUE, names = FALSE, type = 7)
    lower <- q2[1]
    upper <- q2[2]
    # If still equal (extremely degenerate), return all "Medium"
    if (is.na(lower) || is.na(upper) || lower == upper) {
      return(factor(ifelse(is.na(x), NA, "Medium"), levels = out_levels))
    }
  }
  
  # Cut into three bins: (-Inf, lower], (lower, upper), [upper, Inf)
  # include.lowest=TRUE ensures values equal to the global minimum are included
  res <- cut(x,
             breaks = c(-Inf, lower, upper, Inf),
             labels = out_levels,
             right = TRUE,
             include.lowest = TRUE)
  
  # Ensure factor has the desired levels even if some are empty
  res <- factor(as.character(res), levels = out_levels)
  
  return(res)
}


weighted.clus <- function(data,varnames,coefs,clus=3){
  weighted_ERS <- sweep(data[,varnames], 2, coefs, FUN = "*")
  # print(coefs)
  weighted_clusters <- kmeans(weighted_ERS, centers = clus)$cluster
  return(weighted_clusters)
}


##ALL CLUSTERING
clus.inter <- function(y,grs,ers,ersvars=NA,coefs=NA,data,binary=F){
  
  grs_bin <- to_halves(data[,"grs"])
  ers_bin <- to_halves(data[,"ers"])
  ers_tert <- to_tertiles(data[,"ers"])
  ers_quart <- to_quartiles(data[,"ers"])
  ers_mid80 <- split_mid80(data[,"ers"])
  ers_mid90 <- split_mid90(data[,"ers"])
  ersmat <- data[,ersvars]
  ers_clus <- kmeans(ersmat,centers=3)$cluster
  ers_wclus <- weighted.clus(data,ersvars,coefs,clus=3)
  # ADD GMM CLUSTERING HERE
  
  gmm_fit <- Mclust(ersmat, G = 3, verbose=F)
  ers_gmm <- gmm_fit$classification
  
  if(binary==T){
    glm_bin <- glm(y ~ grs_bin + ers_bin+grs_bin*ers_bin, data, family="binomial")
    glm_tert <- glm(y ~ grs_bin + ers_tert+grs_bin*ers_tert, data, family="binomial")
    glm_quart <- glm(y ~ grs_bin + ers_quart+grs_bin*ers_quart, data, family="binomial")
    glm_mid80 <- glm(y ~ grs_bin + ers_mid80+grs_bin*ers_mid80, data, family="binomial")
    glm_mid90 <- glm(y ~ grs_bin + ers_mid90+grs_bin*ers_mid90, data, family="binomial")
    glm_clus <- glm(y ~ grs_bin + factor(ers_clus)+grs_bin* factor(ers_clus), data, family="binomial")
    glm_wclus <- glm(y ~ grs_bin +  factor(ers_wclus)+grs_bin* factor(ers_wclus), data, family="binomial")
    glm_gmm <- glm(y ~ grs_bin + factor(ers_gmm)+grs_bin*factor(ers_gmm), data, family="binomial")
  }else{ 
    glm_bin <- glm(y ~ grs_bin + ers_bin+grs_bin*ers_bin, data=data)
    glm_tert <- glm(y ~ grs_bin + ers_tert+grs_bin*ers_tert, data=data)
    glm_quart <- glm(y ~ grs_bin + ers_quart+grs_bin*ers_quart, data=data)
    glm_mid80 <- glm(y ~ grs_bin + ers_mid80+grs_bin*ers_mid80, data=data)
    glm_mid90 <- glm(y ~ grs_bin + ers_mid90+grs_bin*ers_mid90, data=data)
    glm_clus <- glm(y ~ grs_bin + factor(ers_clus)+grs_bin* factor(ers_clus), data=data)
    glm_wclus <- glm(y ~ grs_bin +  factor(ers_wclus)+grs_bin* factor(ers_wclus), data=data)
    glm_gmm <- glm(y ~ grs_bin + factor(ers_gmm)+grs_bin*factor(ers_gmm), data=data)
  }
  
  results <- cbind(summarize_glm_terms(glm_bin, alpha = 0.05, keep_details = F),
                   summarize_glm_terms(glm_tert, alpha = 0.05, keep_details = F),
                   summarize_glm_terms(glm_quart, alpha = 0.05, keep_details = F),
                   summarize_glm_terms(glm_mid80, alpha = 0.05, keep_details = F),
                   summarize_glm_terms(glm_mid90, alpha = 0.05, keep_details = F),
                   summarize_glm_terms(glm_clus, alpha = 0.05, keep_details = F),
                   summarize_glm_terms(glm_wclus, alpha = 0.05, keep_details = F),
                   summarize_glm_terms(glm_gmm, alpha = 0.05, keep_details = F))
  dim(results)
  colnames(results) <- c("median","tertile","quartile","mid80","mid90","clus","wclus","gmm")
  rownames(results) <- c("ers","gxe","grs")
  #out <- list(pvals_bin,pvals_tert,pvals_quart,pvals_mid80,pvals_mid90,pvals_clus,pvals_wclus)
  return(results)
}


clus.inter2 <- function(y,grs,ers,ersvars=NA,coefs=NA,data,binary=F){ ##raw data for ERS, no SNP - GRS only available
  
  grs_bin <- to_halves(grs)
  ers_bin <- to_halves(ers)
  ers_tert <- to_tertiles(ers)
  ers_quart <- to_quartiles(ers)
  ers_mid80 <- split_mid80(ers)
  ers_mid90 <- split_mid90(ers)
  ersmat <- data[,ersvars]
  ers_clus <- kmeans(ersmat,centers=3)$cluster
  ers_wclus <- weighted.clus(data,ersvars,coefs,clus=3)
  # ADD GMM CLUSTERING HERE
  
  gmm_fit <- Mclust(ersmat, G = 3, verbose=F)
  ers_gmm <- gmm_fit$classification
  
  if(binary==T){
    glm_bin <- glm(y ~ grs_bin + ers_bin+grs_bin*ers_bin, data, family="binomial")
    glm_tert <- glm(y ~ grs_bin + ers_tert+grs_bin*ers_tert, data, family="binomial")
    glm_quart <- glm(y ~ grs_bin + ers_quart+grs_bin*ers_quart, data, family="binomial")
    glm_mid80 <- glm(y ~ grs_bin + ers_mid80+grs_bin*ers_mid80, data, family="binomial")
    glm_mid90 <- glm(y ~ grs_bin + ers_mid90+grs_bin*ers_mid90, data, family="binomial")
    glm_clus <- glm(y ~ grs_bin + factor(ers_clus)+grs_bin* factor(ers_clus), data, family="binomial")
    glm_wclus <- glm(y ~ grs_bin +  factor(ers_wclus)+grs_bin* factor(ers_wclus), data, family="binomial")
    glm_gmm <- glm(y ~ grs_bin + factor(ers_gmm)+grs_bin*factor(ers_gmm), data, family="binomial")
  }else{ 
    glm_bin <- glm(y ~ grs_bin + ers_bin+grs_bin*ers_bin, data=data)
    glm_tert <- glm(y ~ grs_bin + ers_tert+grs_bin*ers_tert, data=data)
    glm_quart <- glm(y ~ grs_bin + ers_quart+grs_bin*ers_quart, data=data)
    glm_mid80 <- glm(y ~ grs_bin + ers_mid80+grs_bin*ers_mid80, data=data)
    glm_mid90 <- glm(y ~ grs_bin + ers_mid90+grs_bin*ers_mid90, data=data)
    glm_clus <- glm(y ~ grs_bin + factor(ers_clus)+grs_bin* factor(ers_clus), data=data)
    glm_wclus <- glm(y ~ grs_bin +  factor(ers_wclus)+grs_bin* factor(ers_wclus), data=data)
    glm_gmm <- glm(y ~ grs_bin + factor(ers_gmm)+grs_bin*factor(ers_gmm), data=data)
  }
  
  results <- cbind(summarize_glm_terms(glm_bin, alpha = 0.05, keep_details = F),
                   summarize_glm_terms(glm_tert, alpha = 0.05, keep_details = F),
                   summarize_glm_terms(glm_quart, alpha = 0.05, keep_details = F),
                   summarize_glm_terms(glm_mid80, alpha = 0.05, keep_details = F),
                   summarize_glm_terms(glm_mid90, alpha = 0.05, keep_details = F),
                   summarize_glm_terms(glm_clus, alpha = 0.05, keep_details = F),
                   summarize_glm_terms(glm_wclus, alpha = 0.05, keep_details = F),
                   summarize_glm_terms(glm_gmm, alpha = 0.05, keep_details = F))
  dim(results)
  colnames(results) <- c("median","tertile","quartile","mid80","mid90","clus","wclus","gmm")
  rownames(results) <- c("ers","gxe","grs")
  #out <- list(pvals_bin,pvals_tert,pvals_quart,pvals_mid80,pvals_mid90,pvals_clus,pvals_wclus)
  return(results)
}

summarize_glm_terms <- function(model, alpha = 0.05, keep_details = FALSE) {
  if (!inherits(model, "glm") && !inherits(model, "lm")) {
    stop("model must be a glm or lm object (or derived).")
  }
  
  # Coefficient table (may have column "Pr(>|z|)" or "Pr(>|t|)")
  sumtab <- summary(model)$coefficients
  if (is.null(sumtab)) stop("Could not extract coefficients from model summary.")
  pcol <- grep("^Pr\\(", colnames(sumtab), value = TRUE)
  if (length(pcol) == 0) {
    stop("Could not find p-value column in summary coefficients.")
  }
  # Use the first matching p-value column
  pcol <- pcol[1]
  coef_pvals <- sumtab[, pcol]
  coef_names <- rownames(sumtab)
  
  # model.matrix mapping from columns to terms
  mm <- model.matrix(model)
  mm_colnames <- colnames(mm)
  assign_idx <- attr(mm, "assign")  # 0 = intercept; others map to terms
  term_labels <- attr(terms(model), "term.labels") # in same order as assign (1..)
  
  # Build mapping from coef name -> term label
  # Note: mm columns correspond to coefficient rows in the same order as coef_names,
  # but summary()'s rows include "(Intercept)" and possibly aliased coefficients.
  # We'll match by name to be robust.
  # Create mapping from mm_colnames -> term
  mm_map <- data.frame(mm_col = mm_colnames, assign = assign_idx, stringsAsFactors = FALSE)
  mm_map$term <- ifelse(mm_map$assign == 0, "(Intercept)", term_labels[mm_map$assign])
  # Exclude intercept column when summarizing terms
  mm_map <- mm_map[mm_map$mm_col != "(Intercept)", , drop = FALSE]
  
  # For each mm column, find matching coef row. If not present (aliased/NA), include NA p-value
  mm_map$coef_name <- sapply(mm_map$mm_col, function(nm) {
    # Names sometimes match exactly to coef_names; otherwise try pattern match.
    if (nm %in% coef_names) return(nm)
    # Try escaping special characters for safe grepping
    matches <- grep(paste0("^", gsub("([\\^\\$\\.\\|\\(\\)\\[\\]\\*\\+\\?\\\\])", "\\\\\\1", nm)), coef_names, value = TRUE)
    if (length(matches) >= 1) return(matches[1])
    # Fallback to NA
    return(NA_character_)
  })
  
  # Create a data.frame keyed by term
  terms <- unique(mm_map$term)
  res_list <- lapply(terms, function(tnm) {
    rows <- mm_map[mm_map$term == tnm, , drop = FALSE]
    coef_names_term <- rows$coef_name
    # get p-values, allowing NAs
    pvals <- sapply(coef_names_term, function(cn) {
      if (is.na(cn)) return(NA_real_)
      if (cn %in% names(coef_pvals)) return(as.numeric(coef_pvals[cn]))
      return(NA_real_)
    }, USE.NAMES = TRUE)
    ncoefs <- length(pvals)
    minp <- if (all(is.na(pvals))) NA_real_ else min(pvals, na.rm = TRUE)
    sig_mask <- !is.na(pvals) & (pvals < alpha)
    any_sig <- any(sig_mask)
    sig_coefs <- if (any_sig) paste(names(pvals)[sig_mask], collapse = ", ") else ""
    # Determine term type
    term_type <- if (grepl(":", tnm, fixed = TRUE)) {
      "interaction"
    } else if (tnm %in% names(model$xlevels)) {
      "factor"
    } else {
      "numeric/other"
    }
    list(
      term = tnm,
      n_coefs = ncoefs,
      min_pvalue = if (is.infinite(minp)) NA_real_ else minp,
      significant = any_sig,
      significant_coefs = sig_coefs,
      coef_details = pvals,
      term_type = term_type
    )
  })
  
  # Turn into data.frame
  df <- do.call(rbind, lapply(res_list, function(x) {
    data.frame(
      #  term = x$term,
      #    n_coefs = x$n_coefs,
      min_pvalue = x$min_pvalue,
      #   significant = x$significant,
      #    significant_coefs = x$significant_coefs,
      #   term_type = x$term_type,
      stringsAsFactors = FALSE
    )
  }))
  # Attach details as list column if requested
  if (keep_details) {
    df$coef_details <- lapply(res_list, function(x) x$coef_details)
  }
  
  # Order rows: significant first, then by min_pvalue
  #df <- df[order(-as.integer(df$significant), df$min_pvalue, na.last = TRUE), ]
  rownames(df) <- NULL
  return(df)
}

summarize_pvalues <- function(pval_list, alpha = 0.05, replace_zero = 1e-16, comparator = NULL) {
  # Validate input
  if (!is.list(pval_list) || length(pval_list) == 0) {
    stop("pval_list must be a non-empty list of data.frames or matrices.")
  }
  
  # Convert each element into a data.frame and long format; keep replication index
  long_df <- purrr::map2_dfr(pval_list, seq_along(pval_list),
                             ~ {
                               df <- as.data.frame(.x)
                               # ensure variables column
                               if (is.null(rownames(df))) {
                                 df$variable <- paste0("V", seq_len(nrow(df)))
                               } else {
                                 df$variable <- rownames(df)
                               }
                               tidyr::pivot_longer(df, cols = -variable, names_to = "method", values_to = "pvalue") %>%
                                 mutate(replication = .y)
                             })
  
  # Clean pvalues: coerce to numeric and handle zeros or NAs
  long_df <- long_df %>%
    mutate(pvalue = as.numeric(pvalue),
           pvalue = ifelse(is.na(pvalue), NA_real_, pvalue),
           pvalue = ifelse(!is.na(pvalue) & pvalue == 0, replace_zero, pvalue),
           sig = ifelse(!is.na(pvalue) & pvalue < alpha, TRUE, FALSE))
  
  # Basic summary stats per variable-method across replications
  base_summary <- long_df %>%
    group_by(variable, method) %>%
    summarize(
      n_rep = sum(!is.na(pvalue)),
      prop_sig = ifelse(n_rep > 0, mean(pvalue < alpha, na.rm = TRUE), NA_real_),
      median_p = ifelse(n_rep > 0, median(pvalue, na.rm = TRUE), NA_real_),
      mean_p = ifelse(n_rep > 0, mean(pvalue, na.rm = TRUE), NA_real_),
      n_sig = sum(sig, na.rm = TRUE),
      .groups = "drop"
    )
  
  summary_df <- base_summary
  
  # If comparator provided and exists in the data, compute comparator-based metrics:
  # - sensitivity = P(method sig | comparator sig) = n_both_sig / n_comp_sig
  # - ppv = P(comparator sig | method sig) = n_both_sig / n_method_sig
  # - agreement = proportion of replications where method and comparator have same sig status
  # - fpr = P(method sig | comparator NOT sig) = n_fp / n_comp_notsig
  # - specificity = 1 - fpr (when defined)
  if (!is.null(comparator)) {
    methods_present <- unique(long_df$method)
    if (!(comparator %in% methods_present)) {
      warning(sprintf("Comparator '%s' not found among methods. Ignoring comparator.", comparator))
    } else {
      # comparator flags per variable x replication
      comp_flags <- long_df %>%
        filter(method == comparator) %>%
        select(variable, replication, comp_sig = sig)
      
      # join comparator flag to the full long_df
      joined <- long_df %>%
        left_join(comp_flags, by = c("variable", "replication")) %>%
        # comp_sig may be NA if comparator pvalue is missing; we avoid counting those replications in joint counts
        mutate(method_sig = sig,
               both_sig = ifelse(!is.na(method_sig) & !is.na(comp_sig), method_sig & comp_sig, FALSE),
               both_notsig = ifelse(!is.na(method_sig) & !is.na(comp_sig), (!method_sig & !comp_sig), FALSE),
               fp = ifelse(!is.na(method_sig) & !is.na(comp_sig) & method_sig & !comp_sig, 1L, 0L)
        )
      
      # compute counts for each variable-method
      comp_counts <- joined %>%
        group_by(variable, method) %>%
        summarize(
          n_rep_joint = sum(!is.na(pvalue) & !is.na(comp_sig)),
          n_comp_sig = sum(comp_sig, na.rm = TRUE),
          n_comp_notsig = sum((!comp_sig) & !is.na(comp_sig), na.rm = TRUE),
          n_method_sig = sum(method_sig, na.rm = TRUE),
          n_both_sig = sum(both_sig, na.rm = TRUE),
          n_both_notsig = sum(both_notsig, na.rm = TRUE),
          n_fp = sum(fp, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(
          sensitivity = ifelse(n_comp_sig > 0, n_both_sig / n_comp_sig, NA_real_),
          ppv = ifelse(n_method_sig > 0, n_both_sig / n_method_sig, NA_real_),
          agreement = ifelse(n_rep_joint > 0, (n_both_sig + n_both_notsig) / n_rep_joint, NA_real_),
          fpr = ifelse(n_comp_notsig > 0, n_fp / n_comp_notsig, NA_real_),
          specificity = ifelse(!is.na(fpr), 1 - fpr, NA_real_)
        )
      
      # join comparator counts back to summary_df
      summary_df <- summary_df %>%
        left_join(comp_counts %>% select(variable, method, n_rep_joint, n_comp_sig, n_comp_notsig, n_method_sig, n_both_sig, n_fp, sensitivity, ppv, agreement, fpr, specificity),
                  by = c("variable", "method"))
    }
  }
  
  # method_summary: per-variable summary (one row per variable-method)
  # This provides an easy-to-inspect table keyed by variable and method with the main metrics.
  method_summary <- summary_df %>%
    arrange(variable, method) %>%
    select(variable, method, n_rep, prop_sig, n_sig, median_p, mean_p,
           starts_with("n_"), sensitivity, ppv, agreement, fpr, specificity)
  print(method_summary)
  list(long = long_df, summary = summary_df, method_summary = method_summary, comparator = comparator)
}

interact.clust.test <- function(alpha,corr,beta_GxE,beta_GRS,beta_ERS,p_env){
  
  #Scenarios
  # noise_snp <- 0
  
  # p_env <- 5
  # beta_GRS <- 1.3
  # beta_ERS <- 1
  # beta_GxE <- 2.5
  # corr <- 0.5 #correlation within env vars
  # alpha <- 0
  # 
  effect_scale_env <- 0.05
  
  #CONSTANTS
  n_reps <- 50
  # maf_low <- 0.25
  # maf_high <- 0.45
  # effect_scale_snp <- 1
  # p_snp <- 50
  # p_causal <- p_snp-noise_snp
  N_samp <- 5000
  # interaction_mode <- "all"
  
  #prevalence <- 0.5#c(0.115, 0.071, 0.047, 0.033, 0.024)
  pvals_bin <- matrix(NA,nrow=n_reps,ncol=3)
  pvals_cont <- matrix(NA,nrow=n_reps,ncol=3)
  colnames(pvals_bin) <- c("grs","ers","grs:ers")
  colnames(pvals_cont) <- c("grs","ers","grs:ers")
  
  load("salia_emp.Rdata")
  outbin <- list()
  outcont <- list()
  for(l in 1:n_reps) {
    if (l %% 10 == 0) {cat("data: ", l,"\n")}
    # Here, simulateSNPglm is only used for simulating the genotypes
    # The phenotype is simulated afterwards with the additional continuous variables
    
    
    
    
    dens.grs <- density(grs.salia)
    
    # Simulate from the KDE
    g_main <- sample(dens.grs$x, N_samp, prob = dens.grs$y, replace = TRUE)
    
    emp.data <- cbind(pm25.salia,pm10.salia)
    kde_fit <- ks::kde(x = emp.data)
    apo <- ks::rkde(n = N_samp, fhat = kde_fit)
    # apo <- apo + matrix(
    #   rnorm(nrow(apo) * ncol(apo), mean = 0, sd = 1e-),
    #   nrow = nrow(apo)
    # )
    dens.no2 <- density(no2.salia)
    no2 <- sample(dens.no2$x, N_samp, prob = dens.no2$y, replace = TRUE)
    ENVdata <- cbind(apo,no2)
    colnames(ENVdata) <- c("ENV1","ENV2","ENV3")
    
    # environmental main effect
    #true_betas_env <- rnorm(3, 0, effect_scale_env)#runif(p_env,-log(1.01),log(1.01))
    true_betas_env <- c(0.3, 0.5, 0.1)  # Fixed, meaningful effects
    e_main <- ENVdata %*% true_betas_env #if (env_main) rnorm(1, 0, 1) * E else rep(0, N_samp)
    
    
    # print(var(g_main))
    # print(var(e_main))
    # print(var(g_main * e_main))
    
    # Center before scaling (avoid spurious interaction)
    g_main_c <- g_main - mean(g_main)
    e_main_c <- e_main - mean(e_main)
    
    # Now scale the centered versions
    g_main <- scale(g_main_c)[,1]
    e_main <- scale(e_main_c)[,1]
    
    # Create the interaction and scale it
    gxe_term <- g_main_c * e_main_c
    gxe_term_scaled <- scale(gxe_term)[,1]  # Now Var(gxe_term_scaled) = 1
    
    # Print diagnostics
    # cat("=== Marginal Statistics ===\n")
    # cat("g_main: mean =", mean(g_main), "var =", var(g_main), "\n")
    # cat("e_main: mean =", mean(e_main), "var =", var(e_main), "\n")
    # cat("gxe_term_scaled: mean =", mean(gxe_term_scaled), "var =", var(gxe_term_scaled), "\n")
    # 
    # # Check distributions
    # cat("\n=== Distribution Shapes ===\n")
    # cat("g_main: skewness =", moments::skewness(g_main), "\n")
    # cat("e_main: skewness =", moments::skewness(e_main), "\n")
    # 
    # # Check correlations
    # cat("\n=== Correlations ===\n")
    # cat("cor(g_main, e_main) =", cor(g_main, e_main), "\n")
    # cat("cor(g_main, gxe_term_scaled) =", cor(g_main, gxe_term_scaled), "\n")
    # cat("cor(e_main, gxe_term_scaled) =", cor(e_main, gxe_term_scaled), "\n")
    # 
    # 
    # # Calculate total signal variance
    # signal_var <- var(beta_GRS*g_main + beta_ERS*e_main + beta_GxE*gxe_term_scaled)
    
    # Set noise to achieve desired SNR (signal-to-noise ratio)
    # desired_snr <- 1  # adjust as needed
    # noise_sd <- sqrt(signal_var / desired_snr)
    noise <- rnorm(N_samp, 0, 1.2)
    
    # Build outcome: main effects use scaled versions, interaction uses centered originals
    y_cont <- beta_GRS*g_main + beta_ERS*e_main + beta_GxE*g_main_c*e_main_c + noise
    
    # Convert to probability WITHOUT scaling (removes the scaling artifact issue)
    p <- plogis(y_cont)
    y_bin <- rbinom(N_samp, size = 1, prob = p)
    cat("\n=== Outcome statistics ===\n")
    cat("var(y_cont):", var(y_cont), "\n")
    cat("range(y_cont):", range(y_cont), "\n")
    cat("mean(p):", mean(p), "\n")
    cat("var(p):", var(p), "\n")
    cat("range(p):", range(p), "\n")
    cat(table(y_bin))
    # cat("Var(beta_GRS*g_main) =", var(beta_GRS*g_main), "\n")
    # cat("Var(beta_ERS*e_main) =", var(beta_ERS*e_main), "\n")
    # cat("Var(beta_GxE*gxe_term_scaled) =", var(beta_GxE*gxe_term_scaled), "\n")
    # cat("Var(noise) =", var(noise), "\n")
    # cat("Total Var(y_cont) =", var(y_cont), "\n")
    
    real.lm <- glm(y_cont~g_main+e_main+g_main*e_main)
    real.glm <- glm(y_bin~g_main+e_main+g_main*e_main,family="binomial")
    pvals_realc <- c(summary(real.lm)$coefficients[,4]["e_main"],
                     summary(real.lm)$coefficients[,4]["g_main:e_main"],
                     summary(real.lm)$coefficients[,4]["g_main"])
    pvals_realb <- c(summary(real.glm)$coefficients[,4]["e_main"],
                     summary(real.glm)$coefficients[,4]["g_main:e_main"],
                     summary(real.glm)$coefficients[,4]["g_main"])
    
    
    
    train_idx <- sample(1:N_samp, round(N_samp*0.5))
    test_idx <- setdiff(1:N_samp, train_idx)
    
    subsets <- list(env=colnames(ENVdata))
    
    
    
    #  BINARYYY
    merged <- scale(cbind(ENVdata))
    dim(merged)
    
    weights0 <- numeric(0)
    weights0norm <- numeric(0)
    fit_rs_bin <- cv.glmnet(merged[train_idx, ], y_bin[train_idx], alpha=alpha,standardize=F,family="binomial")
    train_weights <- as.vector(coef(fit_rs_bin, s = "lambda.min"))[-1]
    weights0 <- rbind(weights0,train_weights)
    colnames(weights0) <- colnames(merged)
    weights0norm <- rbind(weights0norm, weights0)
    colnames(train_weights) <- names(coef(fit_rs_bin, s = "lambda.min"))[-1]
    scores0 <- matrix(NA,nrow=length(train_idx),ncol=length(subsets))
    for(j in seq_along(subsets)){
      # print(j)
      res <- riskscores_sub_abs(subsets[[j]], weights0[,], Xmat=merged[test_idx,])
      scores0[,j] <- res$risk
      weights0norm[,subsets[[j]]] <- res$weights.norm
    }
    
    for(j in seq_along(subsets)){
      if(IQR(scores0[,j])!=0){scores0[,j] <- scores0[,j]/IQR(scores0[,j])
      }else{print(paste("IQR ", subsets[[j]], " = 0"))}
    }
    colnames(scores0) <- c("ers")
    all0 <- apply(scores0,2,f<-function(x){all(x==0)})
    grs <- g_main
    newdata <- data.frame(y=y_bin[test_idx],grs=grs[test_idx],ers=scores0[,!all0])
    glm_result <- glm(y ~ grs + ers+grs*ers, data=newdata, family="binomial")
    # pvals_bin[l,"grs"] <- summary(glm_result)$coefficients[,4]["grs"]
    # pvals_bin[l,"ers"] <- summary(glm_result)$coefficients[,4]["ers"]
    # pvals_bin[l,"grs:ers"] <- summary(glm_result)$coefficients[,4]["grs:ers"]
    test.bin <- c(summary(glm_result)$coefficients[,4]["ers"],
                  summary(glm_result)$coefficients[,4]["grs:ers"],
                  summary(glm_result)$coefficients[,4]["grs"])
    clust.dat <- data.frame(ers=scores0[,!all0],grs=grs[test_idx],ENVdata[test_idx,])
    clus.pvals <- clus.inter(y_bin[test_idx],grs=grs,ers=ers,ersvars=colnames(ENVdata),
                             coefs=weights0norm[,grepl("ENV",colnames(weights0norm))],data=clust.dat,binary=T)
    outbin[[l]] <- data.frame(clus.pvals,test=test.bin,real=pvals_realb)
    
    
    ###
    
    #CONTINUOUS
    weights0 <- numeric(0)
    weights0norm <- numeric(0)
    fit_rs_cont <- cv.glmnet(merged[train_idx, ], y_cont[train_idx], alpha=alpha,standardize=F)
    train_weights <- as.vector(coef(fit_rs_cont, s = "lambda.min"))[-1]
    weights0 <- rbind(weights0,train_weights)
    colnames(weights0) <- colnames(merged)
    weights0norm <- rbind(weights0norm, weights0)
    colnames(train_weights) <- names(coef(fit_rs_cont, s = "lambda.min"))[-1]
    scores0 <- matrix(NA,nrow=length(train_idx),ncol=length(subsets))
    for(j in seq_along(subsets)){
      #print(j)
      res <- riskscores_sub_abs(subsets[[j]], weights0[,], Xmat=merged[test_idx,])
      scores0[,j] <- res$risk
      weights0norm[,subsets[[j]]] <- res$weights.norm
    }
    
    for(j in seq_along(subsets)){
      if(IQR(scores0[,j])!=0){scores0[,j] <- scores0[,j]/IQR(scores0[,j])
      }else{print(paste("IQR ", subsets[[j]], " = 0"))}
    }
    colnames(scores0) <- c("ers")
    all0 <- apply(scores0,2,f<-function(x){all(x==0)})
    newdata <- data.frame(y=y_cont[test_idx],grs=grs[test_idx],ers=scores0[,!all0])
    
    glm_result <- glm(y ~ grs + ers+grs*ers, data=newdata)
    # pvals_cont[l,"grs"] <- summary(glm_result)$coefficients[,4]["grs"]
    # pvals_cont[l,"ers"] <- summary(glm_result)$coefficients[,4]["ers"]
    # pvals_cont[l,"grs:ers"] <- summary(glm_result)$coefficients[,4]["grs:ers"]
    test.cont <- c(summary(glm_result)$coefficients[,4]["ers"],
                   summary(glm_result)$coefficients[,4]["grs:ers"],
                   summary(glm_result)$coefficients[,4]["grs"])
    clust.dat <- data.frame(ers=scores0[,!all0],grs=grs[test_idx],ENVdata[test_idx,])
    clus.pvals <- clus.inter(y_bin[test_idx],grs=grs[test_idx],ers=scores0[,!all0],ersvars=colnames(ENVdata),
                             coefs=weights0norm[,grepl("ENV",colnames(weights0norm))],data=clust.dat,binary=F)
    outcont[[l]] <- data.frame(clus.pvals,test=test.cont,real=pvals_realc)
    
    
  }
  summ.bin.GxE <- summ.bin <- summarize_pvalues(outbin,comparator="real")
  
  summ.cont.GxE <- summ.cont <- summarize_pvalues(outcont,comparator="real")
  
  output <- c(unlist(summ.bin$method_summary %>%   filter(variable == "gxe")%>%
                       select(fpr)),
              unlist(summ.cont$method_summary %>%   filter(variable == "gxe")%>%
                       select(fpr)))
  return(output)
}

#For sensitivity
interact.clust.test2 <- function(alpha,sd,beta_GxE,beta_GRS,beta_ERS,p_env){
  
  #Scenarios
  # noise_snp <- 0
  
  # p_env <- 5
  # beta_GRS <- 1.3
  # beta_ERS <- 1
  # beta_GxE <- 2.5
  # corr <- 0.5 #correlation within env vars
  # alpha <- 0
  # 
  effect_scale_env <- 0.05
  
  #CONSTANTS
  n_reps <- 100
  # maf_low <- 0.25
  # maf_high <- 0.45
  # effect_scale_snp <- 1
  # p_snp <- 50
  # p_causal <- p_snp-noise_snp
  N_samp <- 500
  # interaction_mode <- "all"
  
  #prevalence <- 0.5#c(0.115, 0.071, 0.047, 0.033, 0.024)
  pvals_bin <- matrix(NA,nrow=n_reps,ncol=3)
  pvals_cont <- matrix(NA,nrow=n_reps,ncol=3)
  colnames(pvals_bin) <- c("grs","ers","grs:ers")
  colnames(pvals_cont) <- c("grs","ers","grs:ers")
  
  load("salia_emp.Rdata")
  outbin <- list()
  outcont <- list()
  for(l in 1:n_reps) {
    if (l %% 10 == 0) {cat("data: ", l,"\n")}
    # Here, simulateSNPglm is only used for simulating the genotypes
    # The phenotype is simulated afterwards with the additional continuous variables
    
    
    
    
    dens.grs <- density(grs.salia)
    
    # Simulate from the KDE
    g_main <- sample(dens.grs$x, N_samp, prob = dens.grs$y, replace = TRUE)
    
    emp.data <- cbind(pm25.salia,pm10.salia)
    kde_fit <- ks::kde(x = emp.data)
    apo <- ks::rkde(n = N_samp, fhat = kde_fit)
    # apo <- apo + matrix(
    #   rnorm(nrow(apo) * ncol(apo), mean = 0, sd = 1e-),
    #   nrow = nrow(apo)
    # )
    dens.no2 <- density(no2.salia)
    no2 <- sample(dens.no2$x, N_samp, prob = dens.no2$y, replace = TRUE)
    ENVdata <- cbind(apo,no2)
    colnames(ENVdata) <- c("ENV1","ENV2","ENV3")
    
    # environmental main effect
    true_betas_env <- c(0.3, 0.5, 0.1) #true_betas_env <- rnorm(3, 0, effect_scale_env)#runif(p_env,-log(1.01),log(1.01))
    e_main <- ENVdata %*% true_betas_env #if (env_main) rnorm(1, 0, 1) * E else rep(0, N_samp)
    
    
    
    
    # Center before scaling (avoid spurious interaction)
    g_main_c <- g_main - mean(g_main)
    e_main_c <- e_main - mean(e_main)
    
    # Now scale the centered versions
    g_main <- scale(g_main_c)[,1]
    e_main <- scale(e_main_c)[,1]
    
    
    # Create the interaction and scale it
    gxe_term <- g_main_c * e_main_c
    gxe_term_scaled <- scale(gxe_term)[,1]  # Now Var(gxe_term_scaled) = 1
    
    
    # Print diagnostics
    # cat("=== Marginal Statistics ===\n")
    # cat("g_main: mean =", mean(g_main), "var =", var(g_main), "\n")
    # cat("e_main: mean =", mean(e_main), "var =", var(e_main), "\n")
    # cat("gxe_term_scaled: mean =", mean(gxe_term_scaled), "var =", var(gxe_term_scaled), "\n")
    
    # Calculate total signal variance
    # signal_var <- var(beta_GRS*g_main + beta_ERS*e_main + beta_GxE*g_main_c*e_main_c)
    
    # Set noise to achieve desired SNR (signal-to-noise ratio)
    # desired_snr <- 1  # adjust as needed
    # noise_sd <- sqrt(signal_var / desired_snr)
    noise <- rnorm(N_samp, 0, sd)
    
    
    # Build outcome: main effects use scaled versions, interaction uses centered originals
    y_cont <- beta_GRS*g_main + beta_ERS*e_main + beta_GxE*g_main_c*e_main_c + noise
    
    # Convert to probability WITHOUT scaling (removes the scaling artifact issue)
    p <- plogis(y_cont)
    y_bin <- rbinom(N_samp, size = 1, prob = p)
    
    # cat("Var(beta_GRS*g_main) =", var(beta_GRS*g_main), "\n")
    # cat("Var(beta_ERS*e_main) =", var(beta_ERS*e_main), "\n")
    # cat("Var(beta_GxE*gxe_term_scaled) =", var(beta_GxE*gxe_term_scaled), "\n")
    # cat("Var(noise) =", var(noise), "\n")
    # cat("Total Var(y_cont) =", var(y_cont), "\n")
    
    real.lm <- glm(y_cont~g_main+e_main+g_main*e_main)
    real.glm <- glm(y_bin~g_main+e_main+g_main*e_main,family="binomial")
    pvals_realc <- c(summary(real.lm)$coefficients[,4]["e_main"],
                     summary(real.lm)$coefficients[,4]["g_main:e_main"],
                     summary(real.lm)$coefficients[,4]["g_main"])
    pvals_realb <- c(summary(real.glm)$coefficients[,4]["e_main"],
                     summary(real.glm)$coefficients[,4]["g_main:e_main"],
                     summary(real.glm)$coefficients[,4]["g_main"])
    
    
    
    train_idx <- sample(1:N_samp, round(N_samp*0.5))
    test_idx <- setdiff(1:N_samp, train_idx)
    
    subsets <- list(env=colnames(ENVdata))
    
    
    
    #  BINARYYY
    merged <- scale(cbind(ENVdata))
    dim(merged)
    
    weights0 <- numeric(0)
    weights0norm <- numeric(0)
    fit_rs_bin <- cv.glmnet(merged[train_idx, ], y_bin[train_idx], alpha=alpha,standardize=F,family="binomial")
    train_weights <- as.vector(coef(fit_rs_bin, s = "lambda.min"))[-1]
    weights0 <- rbind(weights0,train_weights)
    colnames(weights0) <- colnames(merged)
    weights0norm <- rbind(weights0norm, weights0)
    colnames(train_weights) <- names(coef(fit_rs_bin, s = "lambda.min"))[-1]
    scores0 <- matrix(NA,nrow=length(train_idx),ncol=length(subsets))
    for(j in seq_along(subsets)){
      # print(j)
      res <- riskscores_sub_abs(subsets[[j]], weights0[,], Xmat=merged[test_idx,])
      scores0[,j] <- res$risk
      weights0norm[,subsets[[j]]] <- res$weights.norm
    }
    
    for(j in seq_along(subsets)){
      if(IQR(scores0[,j])!=0){scores0[,j] <- scores0[,j]/IQR(scores0[,j])
      }else{print(paste("IQR ", subsets[[j]], " = 0"))}
    }
    colnames(scores0) <- c("ers")
    all0 <- apply(scores0,2,f<-function(x){all(x==0)})
    grs <- g_main
    newdata <- data.frame(y=y_bin[test_idx],grs=grs[test_idx],ers=scores0[,!all0])
    glm_result <- glm(y ~ grs + ers+grs*ers, data=newdata, family="binomial")
    # pvals_bin[l,"grs"] <- summary(glm_result)$coefficients[,4]["grs"]
    # pvals_bin[l,"ers"] <- summary(glm_result)$coefficients[,4]["ers"]
    # pvals_bin[l,"grs:ers"] <- summary(glm_result)$coefficients[,4]["grs:ers"]
    test.bin <- c(summary(glm_result)$coefficients[,4]["ers"],
                  summary(glm_result)$coefficients[,4]["grs:ers"],
                  summary(glm_result)$coefficients[,4]["grs"])
    clust.dat <- data.frame(ers=scores0[,!all0],grs=grs[test_idx],ENVdata[test_idx,])
    clus.pvals <- clus.inter(y_bin[test_idx],grs=grs,ers=ers,ersvars=colnames(ENVdata),
                             coefs=weights0norm[,grepl("ENV",colnames(weights0norm))],data=clust.dat,binary=T)
    outbin[[l]] <- data.frame(clus.pvals,test=test.bin,real=pvals_realb)
    
    
    ###
    
    #CONTINUOUS
    weights0 <- numeric(0)
    weights0norm <- numeric(0)
    fit_rs_cont <- cv.glmnet(merged[train_idx, ], y_cont[train_idx], alpha=alpha,standardize=F)
    train_weights <- as.vector(coef(fit_rs_cont, s = "lambda.min"))[-1]
    weights0 <- rbind(weights0,train_weights)
    colnames(weights0) <- colnames(merged)
    weights0norm <- rbind(weights0norm, weights0)
    colnames(train_weights) <- names(coef(fit_rs_cont, s = "lambda.min"))[-1]
    scores0 <- matrix(NA,nrow=length(train_idx),ncol=length(subsets))
    for(j in seq_along(subsets)){
      #print(j)
      res <- riskscores_sub_abs(subsets[[j]], weights0[,], Xmat=merged[test_idx,])
      scores0[,j] <- res$risk
      weights0norm[,subsets[[j]]] <- res$weights.norm
    }
    
    for(j in seq_along(subsets)){
      if(IQR(scores0[,j])!=0){scores0[,j] <- scores0[,j]/IQR(scores0[,j])
      }else{print(paste("IQR ", subsets[[j]], " = 0"))}
    }
    colnames(scores0) <- c("ers")
    all0 <- apply(scores0,2,f<-function(x){all(x==0)})
    newdata <- data.frame(y=y_cont[test_idx],grs=grs[test_idx],ers=scores0[,!all0])
    
    glm_result <- glm(y ~ grs + ers+grs*ers, data=newdata)
    # pvals_cont[l,"grs"] <- summary(glm_result)$coefficients[,4]["grs"]
    # pvals_cont[l,"ers"] <- summary(glm_result)$coefficients[,4]["ers"]
    # pvals_cont[l,"grs:ers"] <- summary(glm_result)$coefficients[,4]["grs:ers"]
    test.cont <- c(summary(glm_result)$coefficients[,4]["ers"],
                   summary(glm_result)$coefficients[,4]["grs:ers"],
                   summary(glm_result)$coefficients[,4]["grs"])
    clust.dat <- data.frame(ers=scores0[,!all0],grs=grs[test_idx],ENVdata[test_idx,])
    clus.pvals <- clus.inter(y_bin[test_idx],grs=grs[test_idx],ers=scores0[,!all0],ersvars=colnames(ENVdata),
                             coefs=weights0norm[,grepl("ENV",colnames(weights0norm))],data=clust.dat,binary=F)
    outcont[[l]] <- data.frame(clus.pvals,test=test.cont,real=pvals_realc)
    
    
  }
  summ.bin.GxE <- summ.bin <- summarize_pvalues(outbin,comparator="real")
  
  summ.cont.GxE <- summ.cont <- summarize_pvalues(outcont,comparator="real")
  
  output <- c(unlist(summ.bin$method_summary %>%   filter(variable == "gxe")%>%
                       select(prop_sig)),
              unlist(summ.cont$method_summary %>%   filter(variable == "gxe")%>%
                       select(prop_sig)))
  return(output)
}
# set.seed(12345)
# test1sc <- interact.clust.test(alpha=0,beta_GxE=1e-16,beta_GRS=0.3,beta_ERS=0.5,p_env=3)
# 
# sc.ex <- expand.grid(alpha=c(0,0.5,0.99),beta_GxE=c(1e-16), beta_GRS=c(-0.3,1e-16,0.3),
#                      beta_ERS=c(-0.3,1e-16,0.3))
# 
# Vec.GxE.Sim<- Vectorize(interact.clust.test,vectorize.args=c("alpha","beta_GxE","beta_GRS",
#                                                              "beta_ERS"))
# #trace("Vec.GxE.Sim", tracer = quote(print(list(...))), print = FALSE)
# try <- Vec.GxE.Sim(alpha=sc.ex$alpha,beta_GxE=sc.ex$beta_GxE,
#                    beta_GRS=sc.ex$beta_GRS,beta_ERS=sc.ex$beta_ERS)
# #GxE exists - get sensitivity
# 
# #GxE non-existent - get fp
# binres <- try[1:10,]
# contres <- try[11:20,]
# rownames(binres) <- rownames(contres) <- c("clus","gmm","median","mid80","mid90","quartile","real","tertile","test","wclus")
# 
# out.bin.final <- cbind(sc.ex,t(binres))
# out.cont.final <- cbind(sc.ex,t(contres))
# save(out.bin.final,file=paste("/home/wincy.reyes@IUF.LAN/epishare/Transfer/Wincy Reyes/sim-fpr-bin2_eff.Rdata"))
# save(out.cont.final,file=paste("/home/wincy.reyes@IUF.LAN/epishare/Transfer/Wincy Reyes/sim-fpr-cont2_eff.Rdata"))
# save(out.bin.final,file=paste("/sim-fpr-bin2_eff-gmm-grssim.Rdata"))
# save(out.cont.final,file=paste("/sim-fpr-cont2_eff-gmm-grssim.Rdata"))


#Only Positives
####SENSITIVITY (true GxE exists)
# 
# test1sn <- interact.clust.test2(alpha=0,beta_GxE=0.1,beta_GRS=0.3,beta_ERS=0.2)
# sn.ex <- expand.grid(alpha=c(0,0.5,0.99),beta_GxE=c(0.5,1,1.5), beta_GRS=c(0.5,1),
# beta_ERS=c(1e-16,0.5))
sn.ex <- expand.grid(alpha=c(0,0.5,0.99),beta_GxE=c(0.1,0.2,0.3,0.4), beta_GRS=c(0,0.2),
                     beta_ERS=c(0,0.2),sd=c(0.8))
#sn.ex
Vec.GxE.Sim2<- Vectorize(interact.clust.test2,vectorize.args=c("alpha","beta_GxE","beta_GRS",
                                                               "beta_ERS","sd"))
#trace("Vec.GxE.Sim", tracer = quote(print(list(...))), print = FALSE)
try_sensi <- Vec.GxE.Sim2(alpha=sn.ex$alpha,beta_GxE=sn.ex$beta_GxE,
                          beta_GRS=sn.ex$beta_GRS,beta_ERS=sn.ex$beta_ERS,sn.ex$sd)

binres2 <- try_sensi[1:10,]
contres2 <- try_sensi[11:20,]
rownames(binres2) <- rownames(contres2) <- c("clus","gmm","median","mid80","mid90","quartile","real","tertile","test","wclus")
out.bin.final2 <- cbind(sn.ex,t(binres2))
out.cont.final2 <- cbind(sn.ex,t(contres2))
save(out.bin.final2,file=paste("/home/wincy.reyes@IUF.LAN/epishare/Transfer/Wincy Reyes/sim-sn-bin2_pos1_500.Rdata"))
save(out.cont.final2,file=paste("/home/wincy.reyes@IUF.LAN/epishare/Transfer/Wincy Reyes/sim-sn-cont2_pos1_500.Rdata"))