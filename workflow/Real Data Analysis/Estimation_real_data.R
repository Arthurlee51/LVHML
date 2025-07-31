#General script to handle covariates. 
# Load required libraries
library(readr)

# Loop over a specific range; 
for (J in c(100)) {
  # Read transaction and demographic data
  transaction_data <- read_csv("Complete journey/transaction_data.csv")
  demographic_data <- read_csv("Complete journey/hh_demographic_extend.csv")
  
  # Set time period constants
  Tp <- 24 # Total number of 4-week periods considered
  Last <- 24 * 4 # Last week number to consider
  
  # Trim transaction data to only include the last Tp periods
  indices <- which(transaction_data$WEEK_NO %in% (Last - Tp * 4 + 1):Last)
  transaction_data_trimmed <- transaction_data[indices, ]
  rm(transaction_data)
  
  # Further trim transaction data to include only households with available demographic information
  indices <- which(transaction_data_trimmed$household_key %in% demographic_data$household_key)
  transaction_data_trimmed <- transaction_data_trimmed[indices, ]
  
  # Group each 4-week period into one timeslot
  breaks <- seq((Last - Tp * 4), Last, by = 4)
  labels <- 1:Tp
  fourweeks <- cut(transaction_data_trimmed$WEEK_NO, breaks = breaks, labels = labels, include.lowest = TRUE)
  fourweeks <- as.numeric(as.character(fourweeks))
  transaction_data_trimmed <- cbind(transaction_data_trimmed, fourweeks)
  
  
  # Identify and sort the most popular J items within the selected period
  product_id_freq <- table(transaction_data_trimmed$PRODUCT_ID)
  sorted_freq <- sort(product_id_freq, decreasing = TRUE)
  
  # Select top J items with highest frequency 
  top <- head(sorted_freq, J)
  top_numbers <- as.numeric(names(top))
  indices <- which(transaction_data_trimmed$PRODUCT_ID %in% top_numbers)
  toanal <- transaction_data_trimmed[indices, c("household_key", "PRODUCT_ID", "fourweeks")]
  
  # Assign a unique numeric ID to each household
  id <- as.integer(factor(toanal$household_key))
  N <- length(unique(id))
  toanal <- cbind(toanal, id)
  toanal <- toanal[order(toanal$id), ]
  
  # Assign a unique numeric event ID for each product
  event_id <- as.integer(factor(toanal$PRODUCT_ID))
  toanal <- cbind(toanal, event_id)
  
  # Keep only the data for households considered in 'toanal'
  transaction_data_trimmed <- transaction_data_trimmed[which(transaction_data_trimmed$household_key %in% toanal$household_key), ]
  
  # Prepare matrices for analysis
  allevents <- 1:J
  to_transtime <- toanal[, c("id", "fourweeks")]
  R <- matrix(0, nrow = N, ncol = Tp)
  
  for (i in 1:N) {
    datai <- to_transtime[which(to_transtime$id == i), ]
    times <- sort(unique(c(datai$fourweeks)))
    R[i, times] <- 1
  }
  
  # Initialize Y, a 3D array to track transactions across individuals, items, and time
  Y <- array(0, dim = c(N, J, Tp))
  
  for (i in 1:N) {
    datai <- toanal[which(toanal$id == i), ]
    obsevent <- unique(datai$event_id)
    
    for (j in obsevent) {
      Y[i, j, datai$fourweeks[datai$event_id == j]] <- 1
    }
    
    Y[i, , which(R[i, ] == 0)] <- NA
  }
  
  # Map unique IDs to the transaction data
  mapping <- unique(toanal[, c("household_key", "id")])
  transaction_data_trimmed <- merge(transaction_data_trimmed, mapping, by = "household_key", all.x = TRUE)
  
  # Trim demographic data to include only relevant households
  indices <- which(demographic_data$household_key %in% transaction_data_trimmed$household_key)
  demographic_data_trimmed <- demographic_data[indices, ]
  demographic_data_trimmed <- merge(demographic_data_trimmed, mapping, by = "household_key", all.x = TRUE)
  
  #-----------------------------------------------------------------------------------------------------------
  #Prepare covariates
  # Initialize and populate the covariate matrix
    px <- 4
 
  #Initialise the covariate matrix X
    X <- matrix(0, nrow = N, ncol = px)
    Income_groups<- demographic_data_trimmed$Income_group
    dummy_vars <- model.matrix(~Income_groups - 1)
    #Middle
    X[,1]<- dummy_vars[,3]
    #High 
    X[,2]<- dummy_vars[,1]

    Household_groups<- demographic_data_trimmed$Household_group
    dummy_vars <- model.matrix(~Household_groups - 1)
    #Middle(2)
    X[,3]<- dummy_vars[,2]
    #High(3+) 
    X[,4]<- dummy_vars[,3]
  
  #----------------------------------------------------------------------------------------------------------
  # Load libraries and source functions for analysis
  library(BB)
  library(expm)
  library(doParallel)
  library(Rcpp)
  library(truncnorm)
  library(glmnet)
  library(Rcpp)
  library(LVHML)
  
  # Set parameters for analysis
  set.seed(123)
  Kset <- 1:15#Candidate set. Use 1:15 for real data analysis.
  par <- 0
  n.cores <- 1
  Z <- array(0, dim = c(0, 0, 0))
  ext = FALSE
  gamma_fix = FALSE
  event_mapping <- unique(toanal[, c("PRODUCT_ID", "event_id")])
  # Perform estimation
  object <- lvhml_est(Y, R,  X , Kset , par, n.cores, Asymp = FALSE, Silent = FALSE, Z ,full = TRUE, ext,gamma_fix)
  # Create and save event mapping for later evaluation
  save(object, N, J, R, X, Tp, mapping,covariates, event_mapping, Y, file = sprintf("FCJ_fixg%dext%dJ%dN%dTp%d.Rdata", gamma_fix,ext,J, N, Tp))
}
