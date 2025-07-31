#This script is to evaluate the recommendation performances of various models based on sensitivity, as in Table 2. 
#Sensitivity of the methods are stored in Sen_result.
#================================================================================================
xi.func = function(x){
  expx = exp(x)
  out=expx/(1 + expx)
  return(out)
}

get_top.func <- function(x,no) {
  sorted_x <- sort(x, decreasing = TRUE)
  top<- head(sorted_x, no)
  return(top)
}

#Function for classification and printing of results
class.func <- function(trained, Y){
  train_con_mat = confusionMatrix(data = factor(trained), reference =  factor(Y), positive = "1" )
  print(c(train_con_mat$overall["Accuracy"], 
          train_con_mat$byClass["Sensitivity"], 
          train_con_mat$byClass["Specificity"]))
}
#================================================================================================
library(readr)
library(lmtest)
library(caret)
library(pROC)
L=4
toexcel = matrix(0, nrow=L, ncol=16)
#Number of recommendations made
norecom <- 10
#Class_result: matrix to store classification result of different recommendation methods give the number of recommendations. For comparison only. Not the main focus.
Class_result <- matrix(NA, nrow = 4, ncol = L*4)
Tp=24
#Result to be copied
#Load the data first. 
Jlong = (1:L)*100
Nlong = rep(800,L)
gamma_fix <- FALSE
ext <- FALSE

for ( l in 1:L){
  #Prepare data for evaluation
  set.seed(123)
  J = Jlong[l]
  N = Nlong[l]
  
  #Read data
  transaction_data <- read_csv("Complete journey/transaction_data.csv")
  demographic_data <- read_csv("Complete journey/hh_demographic_extend.csv")
  
  Last = 24*4
  load(sprintf("FCJ_fixg%dext%dJ%dN%dTp%d.Rdata",gamma_fix*1,ext*1,J,N,Tp))
 
  px = ncol(X)
  Z = array(0, dim=c(0,0,0))
  
  #Indicies of dataset for evaluation:
  indices_eval <- which(transaction_data$WEEK_NO %in% (Last+1):(Last+4) & transaction_data$PRODUCT_ID %in% event_mapping$PRODUCT_ID &transaction_data$household_key %in% mapping$household_key)
  transaction_data_eval <- transaction_data[indices_eval,]
  #Add id to dtransaction_data_eval through the mapping.
  transaction_data_eval <- merge(transaction_data_eval, mapping, by = "household_key", all.x = TRUE)
  #Add event_id to transaction_data_eval through the mapping.
  transaction_data_eval <- merge(transaction_data_eval,event_mapping, by = "PRODUCT_ID", all.x = TRUE)
  
  #Y: Transaction record to be predicted
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
  
  #Now do prediction based on the estimates
  #GET KHAT
  Khat = object$Khat
  Thetahat = object$Thetahat
  g_len <- if(gamma_fix) 1 else Tp
  
  
  #Get parameters and predict probabilities
  Gammahat <- object$Gammahat
  LastGamma = Gammahat[,g_len]
  Ahat <- object$Ahat
  Betahat <- object$Betahat
  if(!gamma_fix){
   Gamma_part <- matrix(LastGamma, nrow=N,ncol=J, byrow=TRUE)
  }else{
    Gamma_part <- (Tp+1)*matrix(LastGamma, nrow=N,ncol=J, byrow=TRUE)
  }
  if(!ext){
    Betapart <- X%*%t(Betahat)
    Apart <- Thetahat%*%t(Ahat)
  }else{
    LastBeta <- Betahat[,((Tp-1)*(px)+1):ncol(Betahat)]
    Betapart <-  X%*%t(LastBeta)
    
    LastA <- Ahat[,((Tp-1)*(Khat)+1):ncol(Ahat)]
    Apart <- Thetahat%*%t(LastA)
  }
  pred_prob <- xi.func(Gamma_part + Betapart+ Apart)
  
  #---------------------------------------------------------------------------------------------------------
  
  
  #Get historical history by summing up
  Y[is.na(Y)] = 0
  histfreq <- apply(Y,c(1,2),sum)
  #Recommendation using history:
  hist_Recommendation = matrix(0, nrow=N,ncol=J)
  for ( i in 1:N){
    histi <- histfreq[i,]
    top_indices <- order(histi, decreasing = TRUE)[1:norecom]
    #Find the item that gets the last value, randomly choose those with ties. 
    last_value = histi[top_indices[norecom]]
    toreplace<-norecom-sum(histi> last_value )
    top_indices[(norecom-toreplace +1):norecom ] <- sample(which(histi==last_value),size = toreplace)
    hist_Recommendation[i,top_indices]<-1
  }
  Class_result[1,(4*l-3):(4*l-1)] <- class.func(hist_Recommendation,Y_eval)
  
  
  #Now consider Recommendation: recommend the items with highest probabilities
  Recommendation = matrix(0, nrow=N,ncol=J)
  for ( i in 1:N){
    # Find the indices of the 20 highest probabilities
    top_indices <- order(pred_prob[i, ], decreasing = TRUE)[1:norecom]
    Recommendation[i, top_indices]<- 1
  }
  # 
  Class_result[2,(4*l-3):(4*l-1)] <-class.func(Recommendation,Y_eval)
  

  
  #Now Hist+popupar_ties(Or Hist-hist in manuscript): hist first, complemented by the most popular items for all ties.
  #histfreq_J: sum up individuals
  histfreq_J = colSums(histfreq)
  #Now fit the prediction using ordered_histfreq_J for individuals with no record.
  Hist_pop_ties_Recommendation = matrix(0, nrow=N, ncol=J)
  for ( i in 1:N){
    histi <- histfreq[i,]
    top_indices <- order(histi, decreasing = TRUE)[1:norecom]
    #Find the item that gets the last value, randomly choose those with ties. 
    last_value = histi[top_indices[norecom]]
    toreplace<-norecom-sum(histi> last_value )
    pred_set <- histfreq_J
    pred_set[which(histi!=last_value)]<-NA
    top_indices[(norecom-toreplace +1):norecom ] <-  order( pred_set, decreasing = TRUE)[1:toreplace]
    Hist_pop_ties_Recommendation[i,top_indices]<-1
  }
  Class_result[3,(4*l-3):(4*l-1)] <- class.func(Hist_pop_ties_Recommendation,Y_eval)

  
  #Comb_ties(or Hist-Prop in manuscript): Recommendation using history supported by latent factor model, now not only for individuals with no purchasing record, but also for bottom ties. 
  comb_ties_Recommendation = matrix(0, nrow=N, ncol=J)
  #pred_prob_comb: Predicted probability for comb
  pred_prob_comb_ties <- pred_prob
  count=0#count number of individuals that used the combination.
  for ( i in 1:N){
    histi <- histfreq[i,]
    top_indices <- order(histi, decreasing = TRUE)[1:norecom]
    #Find the item that gets the last value, randomly choose those with ties. 
    last_value = histi[top_indices[norecom]]
    toreplace<-norecom-sum(histi> last_value )
    pred_prob_comb_ties[i,which(histi!=last_value)]<- NA
    top_indices[(norecom-toreplace +1):norecom ] <-  order( pred_prob_comb_ties[i, ], decreasing = TRUE)[1:toreplace]
    comb_ties_Recommendation[i,top_indices]<-1
  }
  Class_result[4,(4*l-3):(4*l-1)] <- class.func(comb_ties_Recommendation,Y_eval)
}

 #Sensitivity result. Rows 1 to 4 correspond to methods in the order Hist, Prop(depends on ext and gamma_fix), Hist-Hist and Hist_Prop.
 Sen_result <- Class_result[,(1:L)*4-2]
   