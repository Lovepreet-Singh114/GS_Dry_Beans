## 3 approaches
library(tidyverse)
library(rrBLUP)
library(BGLR)
library(metan)

# =============================================================================
# GET REP ID FROM COMMAND LINE
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
rep_id <- 1L

cat("=== Rep:", rep_id, "===\n")

# =============================================================================
# READ DATA
# =============================================================================
blue_path <- if (length(args) >= 2) args[2] else "C:/Users/singhlo2/OneDrive - Michigan State University/Documents/Projects/GS_NIFA/Parametric_Analysis/14. GGE_Informed_NIFA/blues_environment.csv"
blues <- read.csv(blue_path) %>%
  rename(GEN = Name, YIELD = BLUE, ENV = Group) %>%
  select(GEN, ENV, YIELD) %>%
  mutate(GEN = as.character(GEN), ENV = as.character(ENV))

geno_raw <- read.csv("C:/Users/singhlo2/OneDrive - Michigan State University/Documents/Projects/GS_NIFA/Parametric_Analysis/14. GGE_Informed_NIFA/geno_f.csv")
geno_names <- geno_raw[[1]]
M <- as.matrix(geno_raw[, -1])
rownames(M) <- geno_names

M_centered <- M - 1
K <- A.mat(M_centered)

align_geno <- function(blues_sub, K) {
  common <- intersect(unique(blues_sub$GEN), rownames(K))
  blues_sub <- blues_sub %>% filter(GEN %in% common)
  K_sub <- K[common, common]
  list(blues = blues_sub, K = K_sub, geno = common)
}

run_rrblup <- function(y_train, M_train, M_test) {
  fit <- mixed.solve(y = y_train, Z = M_train)
  y_hat <- M_test %*% fit$u + as.numeric(fit$beta)
  as.numeric(y_hat)
}

# =============================================================================
# SUBSAMPLE FRACTION
# =============================================================================
subsample_frac <- 0.8

# =============================================================================
# STRATEGY 0: BASELINE rrBLUP
# =============================================================================
run_baseline_loeo <- function(blues, M_centered, K) {
  envs <- unique(blues$ENV)
  results <- data.frame()
  
  for (test_env in envs) {
    train <- blues %>% filter(ENV != test_env)
    test  <- blues %>% filter(ENV == test_env)
    
    train_means <- train %>%
      group_by(GEN) %>%
      summarise(YIELD = mean(YIELD, na.rm = TRUE), .groups = "drop")
    
    aligned <- align_geno(train_means, K)
    
    # Subsample training genotypes
    n_train <- length(aligned$geno)
    samp_idx <- sort(sample(n_train, round(n_train * subsample_frac)))
    samp_geno <- aligned$geno[samp_idx]
    
    test_geno <- intersect(test$GEN, samp_geno)
    test_sub <- test %>% filter(GEN %in% test_geno)
    
    M_tr <- M_centered[samp_geno, ]
    M_te <- M_centered[test_geno, ]
    y_tr <- aligned$blues$YIELD[match(samp_geno, aligned$blues$GEN)]
    
    pred_rr <- run_rrblup(y_tr, M_tr, M_te)
    obs <- test_sub$YIELD[match(test_geno, test_sub$GEN)]
    acc_rr <- cor(pred_rr, obs, use = "complete.obs")
    
    results <- rbind(results, data.frame(
      Strategy = "Baseline",
      Test_Env = test_env,
      Model = "rrBLUP",
      Accuracy = acc_rr
    ))
  }
  results
}

# =============================================================================
# STRATEGY 1: REACTION NORM rrBLUP
# =============================================================================
run_reaction_norm_loeo <- function(blues, M_centered, K) {
  envs <- unique(blues$ENV)
  results <- data.frame()
  
  for (test_env in envs) {
    train <- blues %>% filter(ENV != test_env)
    test  <- blues %>% filter(ENV == test_env)
    
    common <- intersect(unique(train$GEN), rownames(K))
    common <- intersect(common, unique(test$GEN))
    train <- train %>% filter(GEN %in% common)
    test  <- test %>% filter(GEN %in% common)
    
    # Subsample training genotypes
    train_geno <- unique(train$GEN)
    n_train <- length(train_geno)
    samp_geno <- sort(sample(train_geno, round(n_train * subsample_frac)))
    
    train <- train %>% filter(GEN %in% samp_geno)
    
    # All phenotypic preprocessing is restricted to this sampled training fold.
    train$ENV <- droplevels(factor(train$ENV))
    train$GEN <- droplevels(factor(train$GEN))
    stopifnot(!as.character(test_env) %in% as.character(train$ENV))
    gge_fold <- metan::gge(train, ENV, GEN, YIELD,
                          centering = "environment", scaling = FALSE,
                          svp = "environment")$YIELD
    scores <- as.data.frame(gge_fold$coordenv[, 1:2])
    names(scores) <- c("E_PC1", "E_PC2")
    scores$ENV <- gge_fold$labelenv
    train <- left_join(train, scores, by = "ENV")
    stopifnot(all(is.finite(train$E_PC1)), all(is.finite(train$E_PC2)))
    env_fit <- lm(YIELD ~ E_PC1 + E_PC2, data = train)
    y_resid_train <- residuals(env_fit)
    
    train_resid <- data.frame(GEN = train$GEN, resid = y_resid_train) %>%
      group_by(GEN) %>%
      summarise(resid = mean(resid), .groups = "drop")
    
    geno_order <- as.character(train_resid$GEN)
    M_tr <- M_centered[geno_order, ]
    stopifnot(identical(rownames(M_tr), geno_order))
    
    fit_rr <- mixed.solve(y = train_resid$resid, Z = M_tr)
    
    test_unique <- test %>%
      filter(GEN %in% samp_geno) %>%
      distinct(GEN, .keep_all = TRUE)
    
    marker_pred <- M_centered[as.character(test_unique$GEN), ] %*% fit_rr$u + as.numeric(fit_rr$beta)
    pred_rr <- as.numeric(marker_pred)
    obs <- test_unique$YIELD
    acc_rr <- cor(pred_rr, obs, use = "complete.obs")
    
    results <- rbind(results, data.frame(
      Strategy = "GGE training only",
      Test_Env = test_env,
      Model = "rrBLUP",
      Accuracy = acc_rr
    ))
  }
  results
}


stopifnot(!anyDuplicated(blues[c("GEN", "ENV")]), all(is.finite(blues$YIELD)))
stopifnot(!anyDuplicated(rownames(M_centered)), all(is.finite(M_centered)))
n_reps <- if (length(args)) as.integer(args[1]) else 50L
out <- if (length(args) >= 3) args[3] else "outputs/comment08"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
all_results <- list()
for (rep_id in seq_len(n_reps)) {
  set.seed(rep_id)
  # Preserve original baseline RNG consumption before the corrected strategy.
  baseline <- run_baseline_loeo(blues, M_centered, K)
  corrected <- run_reaction_norm_loeo(blues, M_centered, K)
  result <- bind_rows(baseline, corrected)
  result$Rep <- rep_id
  stopifnot(all(is.finite(result$Accuracy)))
  write.csv(result, file.path(out, paste0("rep_", rep_id, ".csv")), row.names=FALSE)
  all_results[[rep_id]] <- result
  cat("Completed replicate", rep_id, "of", n_reps, "\n")
}
combined <- bind_rows(all_results)
write.csv(combined, file.path(out, "all_results.csv"), row.names=FALSE)
summary <- combined %>% group_by(Strategy, Test_Env) %>%
 summarise(mean_r=mean(Accuracy), sd_r=sd(Accuracy), n_reps=n(), .groups="drop")
write.csv(summary, file.path(out, "summary.csv"), row.names=FALSE)
print(summary)
