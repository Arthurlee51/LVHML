#Functions used only in data analysis.
Get_avg_price.func = function(transaction_data_trimmed){
  Price<- transaction_data_trimmed$SALES_VALUE/transaction_data_trimmed$QUANTITY
  #Set Items where quantity/sales values 0 as NA
  Price[transaction_data_trimmed$QUANTITY==0 ] = NA
  Price[transaction_data_trimmed$SALES_VALUE==0 ] = NA
  transaction_data_trimmed<- cbind(transaction_data_trimmed , Price)
  Average_Price = matrix(0, nrow = J, ncol =Tp)
  #Calculate average price in each period
  for ( t in 1:Tp){
    data_t <- transaction_data_trimmed[transaction_data_trimmed$fourweeks == t ,]
    for( j in 1:J){
      Average_Price[j,t]<- mean(data_t$Price[data_t$event_id==j], na.rm = TRUE)
    }
  }
  return(Average_Price)
}

