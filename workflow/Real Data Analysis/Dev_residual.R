############################################################
# ACJ model deviance analysis (Table S9 + time plots)
#
# Models:
#   ext = 0, gamma_fix = 0 -> "Prop"
#   ext = 1, gamma_fix = 0 -> "Prop (Sect. 2.3.2)"
#   ext = 0, gamma_fix = 1 -> "Prop (Sect. 2.3.3)"
#   ext = 1, gamma_fix = 1 -> "Prop (Sect. 2.3.2 and 2.3.3)"
############################################################

############################
# 1. Helper functions
############################

xi.func <- function(x) {
  expx <- exp(x)
  out  <- expx / (1 + expx)
  return(out)
}

# Convert J × (Tp+1) deviance matrix to long format
# One row = (item, time, J, gamma_fix, ext, model, deviance)
make_dev_long <- function(M, J, gamma_fix, ext) {
  df <- as.data.frame(M)
  df$item <- 1:J
  
  df <- tidyr::pivot_longer(
    df,
    cols      = -item,
    names_to  = "time",
    values_to = "deviance"
  )
  df$time <- as.numeric(gsub("V", "", df$time))
  
  df$J <- J
  
  gamma_fix_val <- gamma_fix
  ext_val       <- ext
  
  df$gamma_fix <- factor(
    gamma_fix_val,
    levels = c(0, 1),
    labels = c("gamma_time", "gamma_fixed")
  )
  df$ext <- factor(
    ext_val,
    levels = c(0, 1),
    labels = c("no_ext", "ext")
  )
  
  # Manuscript model labels
  df$model <- dplyr::case_when(
    ext_val == 0 & gamma_fix_val == 0 ~ "Prop",
    ext_val == 1 & gamma_fix_val == 0 ~ "Prop (Sect. 2.3.2)",
    ext_val == 0 & gamma_fix_val == 1 ~ "Prop (Sect. 2.3.3)",
    ext_val == 1 & gamma_fix_val == 1 ~ "Prop (Sect. 2.3.2 and 2.3.3)"
  )
  
  df
}

############################
# 2. Libraries
############################

library(readr)
library(lmtest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)

############################
# 3. Global settings
############################

L  <- 4
Tp <- 24

Jlong <- (1:L) * 100        # 100, 200, 300, 400
Nlong <- rep(800, L)

demo       <- "HI"          # "HI": Household and Income
gamma_vals <- c(0, 1)
ext_vals   <- c(0, 1)

# Collect deviance for all (J, gamma_fix, ext) in long format
all_dev_long <- list()
idx <- 1

############################
# 4. Main loop over models and J
############################

for (gamma_fix in gamma_vals) {
  for (ext in ext_vals) {
    for (l in 1:L) {
      set.seed(123)
      J <- Jlong[l]
      N <- Nlong[l]
      
      #---------------------------------------------------------
      # Read data
      #---------------------------------------------------------
      transaction_data <- read_csv("Complete journey/transaction_data.csv")
      demographic_data <- read_csv("Complete journey/hh_demographic_extend.csv")
      
      Last <- 24 * 4
      
      #Load data
      load(sprintf( "FCJ_fixg%dext%dJ%dN%dTp%d.Rdata",gamma_fix*1,ext*1,J,N,Tp ) )
      
      px <- ncol(X)
      Z  <- array(0, dim = c(0, 0, 0))   # not used but kept for completeness
      
      #---------------------------------------------------------
      # Build evaluation data (out-of-sample)
      #---------------------------------------------------------
      indices_eval <- which(
        transaction_data$WEEK_NO %in% (Last + 1):(Last + 4) &
          transaction_data$PRODUCT_ID %in% event_mapping$PRODUCT_ID &
          transaction_data$household_key %in% mapping$household_key
      )
      
      transaction_data_eval <- transaction_data[indices_eval, ]
      transaction_data_eval <- merge(transaction_data_eval, mapping,
                                     by = "household_key", all.x = TRUE)
      transaction_data_eval <- merge(transaction_data_eval, event_mapping,
                                     by = "PRODUCT_ID", all.x = TRUE)
      
      # For R_eval: which households are observed in evaluation window
      indices_eval_R <- which(
        transaction_data$WEEK_NO %in% (Last + 1):(Last + 4) &
          transaction_data$household_key %in% mapping$household_key
      )
      
      transaction_data_eval_R <- transaction_data[indices_eval_R, ]
      transaction_data_eval_R <- merge(transaction_data_eval_R, mapping,
                                       by = "household_key", all.x = TRUE)
      
      Obs_id <- unique(transaction_data_eval_R$id)
      R_eval <- matrix(0, nrow = N, ncol = 1)
      R_eval[Obs_id, 1] <- 1
      
      # Y: transaction record to be predicted (out-of-sample)
      Y_eval <- matrix(0, nrow = N, ncol = J)
      toeval <- transaction_data_eval[, c("id", "event_id")]
      
      for (i in 1:N) {
        datai    <- toeval[toeval$id == i, ]
        obsevent <- unique(datai$event_id)
        Y_eval[i, obsevent] <- 1
      }
      
      #---------------------------------------------------------
      # Extract fitted parameters
      #---------------------------------------------------------
      cand  <- object$all_output
      Khat  <- object$Khat
      g_len <- if (gamma_fix) 1 else Tp
      
      # Deviance: in-sample (1:Tp) + out-of-sample (Tp+1)
      Deviance <- matrix(0, nrow = J, ncol = Tp + 1)
      
      Ahat     <- object$Ahat
      Betahat  <- object$Betahat
      Gammahat <- object$Gammahat
      Thetahat <- object$Thetahat
      
      #---------------------------------------------------------
      # In-sample deviance (t = 1..Tp)
      #---------------------------------------------------------
      for (t in 1:Tp) {
        if (!gamma_fix) {
          Gamma_part <- matrix(Gammahat[, t], nrow = N, ncol = J, byrow = TRUE)
        } else {
          Gamma_part <- (t) * matrix(Gammahat[, 1], nrow = N, ncol = J, byrow = TRUE)
        }
        
        if (!ext) {
          Betapart <- X %*% t(Betahat)
          Apart    <- Thetahat %*% t(Ahat)
        } else {
          Betat    <- Betahat[, ((t - 1) * px + 1):(t * px)]
          Betapart <- X %*% t(Betat)
          
          At    <- Ahat[, ((t - 1) * Khat + 1):(t * Khat)]
          Apart <- Thetahat %*% t(At)
        }
        
        pred_prob <- xi.func(Gamma_part + Betapart + Apart)
        
        # Deviance contribution at (i, j, t)
        Yt      <- Y[, , t]
        dev_mat <- -2 * (Yt * log(pred_prob) + (1 - Yt) * log(1 - pred_prob))
        
        # Sum over i → item-level deviance at time t
        Deviance[, t] <- colSums(dev_mat, na.rm = TRUE)
      }
      
      #---------------------------------------------------------
      # Out-of-sample deviance (t = Tp+1)
      #---------------------------------------------------------
      t <- Tp + 1
      if (gamma_fix) {
        Gamma_part <- (Tp + 1) * matrix(Gammahat[, 1], nrow = N, ncol = J, byrow = TRUE)
      }
      # For gamma_time models we reuse Betapart/Apart from t = Tp
      
      pred_prob <- xi.func(Gamma_part + Betapart + Apart)
      Yt        <- Y_eval
      
      dev_mat <- -2 * (Yt * log(pred_prob) + (1 - Yt) * log(1 - pred_prob))
      Deviance[, t] <- colSums(dev_mat, na.rm = TRUE)
      
      #---------------------------------------------------------
      # Store deviance for this (J, gamma_fix, ext) in long format
      #---------------------------------------------------------
      data_dev <- make_dev_long(Deviance, J, gamma_fix, ext)
      all_dev_long[[idx]] <- data_dev
      idx <- idx + 1
    }
  }
}

############################
# 5. Bind all deviance data
############################

heatmap_df <- dplyr::bind_rows(all_dev_long)

# Consistent model order
heatmap_df$model <- factor(
  heatmap_df$model,
  levels = c("Prop",
             "Prop (Sect. 2.3.2)",
             "Prop (Sect. 2.3.3)",
             "Prop (Sect. 2.3.2 and 2.3.3)")
)

Tp <- 24  # in-sample horizon

############################
# 6. Deviance over time
############################

# Sum deviance over items to get total deviance at each time point
dev_time_df <- heatmap_df %>%
  group_by(J, gamma_fix, ext, model, time) %>%
  summarise(
    dev_total = sum(deviance, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(
    sample_type = ifelse(time <= Tp, "in-sample", "out-of-sample")
  )

############################
# 7. Total deviance summary 
############################

dev_summary <- dev_time_df %>%
  group_by(J, model, sample_type) %>%
  summarise(
    dev_sum = sum(dev_total, na.rm = TRUE),
    .groups = "drop"
  )

dev_table_all <- dev_summary %>%
  arrange(model, J, sample_type) %>%
  tidyr::pivot_wider(
    names_from  = c(J, sample_type),
    names_glue  = "J{J}_{sample_type}",
    values_from = dev_sum
  )

dev_table_all 

############################
# 8. Line plot – deviance over time
############################

line_plot <- ggplot(
  dev_time_df,
  aes(x = time, y = dev_total, color = model, group = model)
) +
  geom_line(size = 0.6) +
  geom_point(size = 0.9) +
  geom_vline(
    xintercept = Tp + 0.5,
    linetype   = "dashed"
  ) +
  facet_wrap(~ sprintf("J=%d",J), scales = "free_y") +
  labs(
    title = "Residual Deviance Over Time",
    x     = "Time point",
    y     = "Residual deviance (summed over items)",
    color = "Model"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 9)
  )