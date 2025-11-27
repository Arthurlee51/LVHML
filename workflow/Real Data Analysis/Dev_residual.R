#Script to compute deviance residual, presented as Table S9 in the online supplementary material. 
#Plots of heatmaps of deviance residuals in various models across different settings (4 models and J = 100,200,300,400). Expect a 4 times 4 heatmap. 
xi.func = function(x){
  expx = exp(x)
  out=expx/(1 + expx)
  return(out)
}

make_dev_long <- function(M, J, gamma_fix, ext) {
  df <- as.data.frame(M)
  df$item <- 1:J
  
  df <- pivot_longer(
    df,
    cols      = -item,
    names_to  = "time",
    values_to = "deviance"
  )
  df$time <- as.numeric(gsub("V", "", df$time))
  
  df$J <- J
  df$gamma_fix <- factor(gamma_fix, levels = c(0, 1),
                         labels = c("gamma_time", "gamma_fixed"))
  df$ext <- factor(ext, levels = c(0, 1),
                   labels = c("no_ext", "ext"))
  df$model <- interaction(df$gamma_fix, df$ext, sep = " + ")

  df
}


#This script is to evaluate the predictive ability of the ACJ models under different settings.
library(readr)
library(lmtest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)
L=4
Tp=24
#Result to be copied
#Load the data first. 
Jlong = (1:L)*100
Nlong = rep(800,L)


#"HI": Household and Income
demo =  "HI"
gamma_vals <- c(0, 1)
ext_vals   <- c(0, 1)
Deviance_all <- data.frame(matrix(0,nrow = L*4*2, ncol = 4))
colnames(Deviance_all) <- c("Deviance","J","method","type")
idx <- 1
for (gamma_fix in gamma_vals) {
  for (ext in ext_vals) {
    for ( l in 1:L){
      set.seed(123)
      J = Jlong[l]
      N = Nlong[l]
      
      #Read data
      transaction_data <- read_csv("Complete journey/transaction_data.csv")
      demographic_data <- read_csv("Complete journey/hh_demographic_extend.csv")
      
      Last = 24*4
      #Load data
      if(!ext){
        if(!gamma_fix){
          load(sprintf(   "~/Library/CloudStorage/OneDrive-LondonSchoolofEconomics/codes/MTDRC/Code/Discrete survival/v12/Data analysis/FCJ_fixg%dext%dJ%dN%dTp%d.Rdata",gamma_fix*1,ext*1,J,N,Tp ) )
        }else{
          load(sprintf(   "~/Library/CloudStorage/OneDrive-LondonSchoolofEconomics/codes/MTDRC/Code/Discrete survival/v12/Data analysis/FCJ_fixg%dext%dJ%dN%dTp%d_FULL%s.Rdata",gamma_fix*1,ext*1,J,N,Tp,demo ) )
        }
        
      }else{
        load(sprintf(   "~/Library/CloudStorage/OneDrive-LondonSchoolofEconomics/codes/MTDRC/Code/Discrete survival/v12/Data analysis/altICFCJ_fixg%dext%dJ%dN%dTp%d_FULL%s.Rdata",gamma_fix*1,ext*1,J,N,Tp,demo ) )
      }
      
      px = ncol(X)
      Z = array(0, dim=c(0,0,0))
      
      #Indicies of dataset for evaluation:
      indices_eval <- which(transaction_data$WEEK_NO %in% (Last+1):(Last+4) & transaction_data$PRODUCT_ID %in% event_mapping$PRODUCT_ID &transaction_data$household_key %in% mapping$household_key)
      transaction_data_eval <- transaction_data[indices_eval,]
      #Add id to dtransaction_data_eval through the mapping.
      transaction_data_eval <- merge(transaction_data_eval, mapping, by = "household_key", all.x = TRUE)
      #Add event_id to transaction_data_eval through the mapping.
      transaction_data_eval <- merge(transaction_data_eval,event_mapping, by = "PRODUCT_ID", all.x = TRUE)
      
      #Indices of dataset for computation of R: 
      indices_eval_R <- which(transaction_data$WEEK_NO %in% (Last+1):(Last+4)  &transaction_data$household_key %in% mapping$household_key)
      transaction_data_eval_R <- transaction_data[indices_eval_R,]
      #Add id to dtransaction_data_eval through the mapping.
      transaction_data_eval_R <- merge(transaction_data_eval_R, mapping, by = "household_key", all.x = TRUE)
      #Get unique householdkeys
      Obs_id <- unique(transaction_data_eval_R$id)
      R_eval <- matrix(0, nrow = N,ncol = 1)
      R_eval[Obs_id,1] <- 1
      
      #Y: Transaction record to be predicted (out-sample)
      Y_eval = matrix(0, nrow=N, ncol=J)  
      toeval <- transaction_data_eval[,c("id","event_id")  ]
      for ( i in 1:N){
        #Extract the part with the i-th id
        datai<- toeval[which(toeval$id==i),]
        #id for observed events 
        obsevent= unique(datai$event_id)
        #Assign value to yrijs and rijs for obsevent
        Y_eval[i,obsevent]=1
      }
      
      #Get R for evalutaion data as well.
      
      #Now do prediction based on the estimates
      #Get parameters
      cand = object$all_output
      #GET KHAT
      Khat = object$Khat
      #Old version:
      #Khat = object$Khat2
      
      g_len <- if(gamma_fix) 1 else Tp
      #Initialise deviance (In-sample (1:24) + Out-sample(25))
      Deviance <- matrix(0, nrow = J, ncol = Tp+1)
      Ahat <- object$Ahat
      Betahat <- object$Betahat
      Gammahat <- object$Gammahat
      Thetahat <- object$Thetahat
      #==============================================================================================================
      #In-sample deviance and likelihood.
      #Get estimates and 
      Deviance <- 0
      for( t in 1:Tp){
        if(!gamma_fix){
          Gamma_part <- matrix(  Gammahat[,t], nrow=N,ncol=J, byrow=TRUE)
        }else{
          Gamma_part <- (t)*matrix( Gammahat[,1], nrow=N,ncol=J, byrow=TRUE)
        }
        
        if(!ext){
          Betapart <- X%*%t(Betahat)
          Apart <- Thetahat%*%t(Ahat)
        }else{
          Betat <- Betahat[,((t-1)*(px)+1):((t)*(px))]
          Betapart <-  X%*%t(Betat )
          
          At <- Ahat[,((t-1)*(Khat)+1):((t)*(Khat))]
          Apart <- Thetahat%*%t(At)
        }
        pred_prob <- xi.func(Gamma_part + Betapart+ Apart)
        
        #Get Yt and compute variance.
        Yt <- Y[,,t]
        dev_mat <- -2*(Yt*log(pred_prob) + (1-Yt)*log(1 - pred_prob) )
        Deviance <-  Deviance + sum(dev_mat, na.rm = TRUE)
      }
      Deviance_all[idx,] <- c(Deviance,J,sprintf("ext%dgamma_fix%d",ext,gamma_fix),"in-sample")
      idx <- idx + 1
      #==============================================================================================================
      #Out-sample predictions. No need to compute again, except slightly adjust for gamma.
      t <- Tp+1
      if(gamma_fix){ Gamma_part <- (Tp+1)*matrix( Gammahat[,1], nrow=N,ncol=J, byrow=TRUE)}
      pred_prob <- xi.func(Gamma_part + Betapart+ Apart)
      Yt <- Y_eval
      dev_mat <- -2*(Yt*log(pred_prob) + (1-Yt)*log(1 - pred_prob) )
      Deviance <-  sum(dev_mat, na.rm = TRUE)
      Deviance_all[idx,] <- c(Deviance,J,sprintf("ext%dgamma_fix%d",ext,gamma_fix),"out-sample")
      idx <- idx + 1
      #---------------------------------------------------------------------------------------------------------
      #========================================================
      #========================================================
    }
  }
}


Deviance_tab <- Deviance_all %>%
  mutate(
    # nice model labels for the rows
    model = recode(method,
                   "ext0gamma_fix0" = "Prop",
                   "ext1gamma_fix0" = "Prop (sect 2.3.2)",
                   "ext0gamma_fix1" = "Prop (sect 2.3.3)",
                   "ext1gamma_fix1" = "Prop (sect 2.3.2 & 2.3.3)"
    ),
    J    = factor(J, levels = c(100, 200, 300, 400)),
    type = factor(type, levels = c("in-sample", "out-sample"))
  )

## In-sample block (4 cols: J = 100, 200, 300, 400)
tab_in <- Deviance_tab %>%
  filter(type == "in-sample") %>%
  select(model, J, Deviance) %>%
  pivot_wider(
    names_from  = J,
    values_from = Deviance,
    names_prefix = "J = "
  )

## Out-of-sample block (4 cols: J = 100, 200, 300, 400)
tab_out <- Deviance_tab %>%
  filter(type == "out-sample") %>%
  select(model, J, Deviance) %>%
  pivot_wider(
    names_from  = J,
    values_from = Deviance,
    names_prefix = "J = "
  )



