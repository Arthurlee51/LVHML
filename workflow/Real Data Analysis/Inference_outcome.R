#Script to produce inference outcome for real data analysis. Specifically, the object product_data_size corresponds to the results in Table S6, and product_data_income corresponds to the results in Table S7. 
library(readr)
library(dplyr)
library(LVHML)

source("functions_data_analysis.R")
J=100
ext = FALSE
gamma_fix = FALSE
#Load data
load(sprintf("FCJ_fixg0ext0J%dN800Tp24.Rdata",J))

#Get estimate based on the Khat chosen by IC.
out <- object$all_output[[object$Khat]]


#Extract Betahat and Gammahats
Z = array(0, dim = c(0,0,0))
px <- ncol(X) 
Betahat<- out$Uhat[,(Tp+1):(Tp+px)]
Gammahat <- out$Uhat[, 1:Tp]
#For reference: 
#X1:Middle income
#X2:High income
#X3:Middle household size
#X4:Large household size

#load transaction data first 
transaction_data <- read_csv("Complete journey/transaction_data.csv")
Last = 24*4

indices <- which(transaction_data$WEEK_NO %in% (Last - Tp * 4 + 1):Last)
transaction_data_trimmed <- transaction_data[indices, ]
rm(transaction_data)

# Further trim transaction data to include only households with available demographic information
indices <- which(transaction_data_trimmed$household_key %in% mapping$household_key)
transaction_data_trimmed <- transaction_data_trimmed[indices, ]

indices <- which(transaction_data_trimmed$PRODUCT_ID %in% event_mapping$PRODUCT_ID)
transaction_data_trimmed <- transaction_data_trimmed[indices, ]

# Group each 4-week period into one timeslot
breaks <- seq((Last - Tp * 4), Last, by = 4)
labels <- 1:Tp
fourweeks <- cut(transaction_data_trimmed$WEEK_NO, breaks = breaks, labels = labels, include.lowest = TRUE)
fourweeks <- as.numeric(as.character(fourweeks))
transaction_data_trimmed <- cbind(transaction_data_trimmed, fourweeks)
transaction_data_trimmed <- merge(transaction_data_trimmed,event_mapping, by = "PRODUCT_ID", all.x = TRUE)


#Get average price of the items in each time periods for future comparisons.
Average_Price<- Get_avg_price.func(transaction_data_trimmed)
#==========================================================================
#Find products associated with different covariates based on BY procedure.
#Prepare product data first for later comparisons. 
product_data <- read_csv("Complete journey/product.csv")
product_data_trimmed<- product_data[product_data$PRODUCT_ID %in% event_mapping$PRODUCT_ID, ]
product_data_trimmed <- merge(product_data_trimmed, event_mapping , by = "PRODUCT_ID", all.x = TRUE)



#FDR_control
Thetahat <- out$Thetahat
Uhat <- out$Uhat
#Get asymptotic variances
AsymVarhatstore <-LVHML:::calasympvar.func(X, Z, Thetahat, Uhat, R, J, N, Tp, ext,gamma_fix)


#Wald test version: Hypothesis:
#I %*% (beta1, beta2)^{\top} = (0,0)^{\top} (Or beta 3 , beta 4)
#Compute test statistics
#Get Betahat_var first
Beta_var_store <- AsymVarhatstore[(Tp+1):(Tp+px),(Tp+1):(Tp+px),]
Test_stat_wald<- matrix(0, nrow = J, ncol = 2)
for ( j in 1:J){
  Beta_var <- Beta_var_store[,,j]
  #Test statistic for beta1, beta2
  Test_stat_wald[j,1]<- Betahat[j,1:2]%*%solve(Beta_var[1:2,1:2]/N)%*%Betahat[j,1:2]
  #Test statistic for beta3, beta4
  Test_stat_wald[j,2]<- Betahat[j,3:4]%*%solve(Beta_var[3:4,3:4]/N)%*%Betahat[j,3:4]
}
#Get p-values by comparing to the relevant chi^2 distribution, calculate them separetely. 
p_values<- 1- pchisq(Test_stat_wald ,df=2)
#Adjust each of the p-values using FDR, BH procedure
#Adjust p_valus separately. 
adjusted_p_values_income <- p.adjust(p_values[,1], method = "BY") 
adjusted_p_values_size <- p.adjust(p_values[,2], method = "BY") 
original_p_values_income <- p_values[,1]
original_p_values_size <- p_values[,2]

#Get Significance matrix:
  alpha <- 0.05
  Sig_matrix <- matrix(0, nrow = J, ncol = px)
  Sig_matrix[,1:2]<- rep(adjusted_p_values_income<= alpha ,2)
  Sig_matrix[,3:4]<- rep(adjusted_p_values_size<= alpha ,2)
Betahat_show <- Betahat
Betahat_show[!Sig_matrix]<-NA

#Compute and add overall_average_price to product_data_trimmed
Overall_Average_Price <-rowMeans(Average_Price, na.rm = TRUE)
Overall_var_Price <- apply(Average_Price,c(1),function(x) var(x, na.rm = TRUE) )
Price_mapping <- data.frame(cbind(Overall_Average_Price,Overall_var_Price,Betahat_show,1:J,adjusted_p_values_income,adjusted_p_values_size,original_p_values_income,original_p_values_size ) )

names(Price_mapping)[3:7] <- c("b1:in_mid","b2:in_high","b3:hou_2","b4:hou_3+","event_id")
names(Price_mapping)[8:9] <- c("p_in","p_size")
names(Price_mapping)[10:11] <- c("ori_p_in","ori_p_size")
product_data_trimmed <-merge(product_data_trimmed, Price_mapping , by = "event_id", all.x = TRUE)



# Reorder the dataset by frequency of COMMODITY_DESC and then alphabetically by COMMODITY_DESC
product_data_ordered <- product_data_trimmed %>%
  group_by(COMMODITY_DESC) %>%
  mutate(freq = n()) %>%
  ungroup() %>%
  arrange(desc(freq), COMMODITY_DESC,SUB_COMMODITY_DESC ) # Now also sorting alphabetically by COMMODITY_DESC



#Create bigger categories for analysis
product_data_ordered <- product_data_ordered %>%
  mutate(Broad_Category = case_when(
    COMMODITY_DESC %in% c("APPLES", "BERRIES", "CITRUS", "GRAPES", "MELONS", "TROPICAL FRUIT") ~ "Fruits",
    COMMODITY_DESC %in% c("BROCCOLI/CAULIFLOWER", "CARROTS", "CORN", "MUSHROOMS", "ONIONS", "PEPPERS-ALL", "POTATOES", "SALAD MIX", "SQUASH", "TOMATOES", "VEGETABLES - ALL OTHERS", "VEGETABLES SALAD","SALAD BAR","VEGETABLES - SHELF STABLE") ~ "Vegetables",
    COMMODITY_DESC %in% c("FLUID MILK PRODUCTS", "MILK BY-PRODUCTS", "CHEESE", "EGGS") ~ "Dairy and Eggs",
    COMMODITY_DESC %in% c("BEEF", "CHICKEN", "DELI MEATS", "HOT DOGS") ~ "Meats",
    COMMODITY_DESC %in% c("BAKED BREAD/BUNS/ROLLS", "BREAKFAST SWEETS") ~ "Breakfast",
    COMMODITY_DESC %in% c("REFRGRATD JUICES/DRNKS", "SOFT DRINKS") ~ "Beverages",
    COMMODITY_DESC %in% c("BAG SNACKS", "CANDY - CHECKLANE", "CRACKERS/MISC BKD FD") ~ "Miscellaneous",
    COMMODITY_DESC %in% c("BATH TISSUES", "NEWSPAPER", "PAPER TOWELS") ~ "Miscellaneous",
    TRUE ~ "Miscellaneous" # Catch-all for any categories not explicitly listed above
  ))

# Desired order for the 'category' column
desired_order <- c("Beverages", "Breakfast", "Dairy and Eggs", "Fruits", "Vegetables")

all_categories <- unique(product_data_ordered$Broad_Category)
additional_categories <- setdiff(all_categories, desired_order)
full_order <- c(desired_order, additional_categories)


#Get product ordered for income and household sizes respectively.
product_data_income <- product_data_ordered[c("Broad_Category","SUB_COMMODITY_DESC","Overall_Average_Price","CURR_SIZE_OF_PRODUCT","b1:in_mid","b2:in_high","p_in","ori_p_in")]
#Trim those with p-value larger than 0.05
product_data_income <- product_data_income[which(product_data_income$p_in <=alpha),]


product_data_income$Broad_Category <- factor(product_data_income$Broad_Category, levels = full_order, ordered = TRUE)
#Now order it accoridng to category 
product_data_income <- product_data_income %>%
  arrange(Broad_Category, SUB_COMMODITY_DESC) 


#Get product ordered for income and household sizes respectively.
product_data_size <- product_data_ordered[c("Broad_Category","SUB_COMMODITY_DESC","Overall_Average_Price","CURR_SIZE_OF_PRODUCT","b3:hou_2","b4:hou_3+","p_size","ori_p_size")]
#Trim those with p-value larger than 0.05
product_data_size <- product_data_size[which(product_data_size$p_size <=alpha),]

product_data_size$Broad_Category <- factor(product_data_size$Broad_Category, levels = full_order, ordered = TRUE)
#Now order it accoridng to category 
product_data_size <- product_data_size %>%
  arrange(Broad_Category, SUB_COMMODITY_DESC)
  

