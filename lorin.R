lorin = function(x,y,z){
  
  ### Checking for non-AP ###
  if(x < 1 || y < 1){
    return(0)
  }
  
  ### make sure x < y ###
  maxer = max(x,y)
  miner = min(x,y)
  
  x = miner
  y = maxer
  
  ### Parameters Set Up ###
  n = x + y + z
  
  ### Vectors to Save Results ###
  Prob = c();
  E_k = c()
  
  ### Compute the \sum_k k*P[D = k] ###
  for(k in 2:(z + y + 1)){
    
    sum_j = 0
    
    for(j in 1:(k-1)){
      term1 = ifelse(z>=(k-1-j),choose(k-1,j)*perm(z,k-1-j),0)
      term2 = ifelse(x>=j,y*perm(x,j),0) + ifelse(y>=j,x*perm(y,j),0)
      
      sum_j = sum_j + (term1*term2)
      
      # cat("\n\n for j = ", j,
      #     "\n term1 = ", term1,
      #     "\n term2 = ", term2,
      #     "\n sum_j = ", sum_j, sep = "")
    }
    Prob = c(Prob,sum_j/perm(n,k)) #Probability for Draw k
    E_k = c(E_k,k*Prob[k-1]) #Value for Draw k 
    
    # cat("\n\nfor k = ", k, ", Tracked Expectation = ", sum(E_k), sep = "")
    
  }
  
  ### Compute the Expectation ###
  Exp = sum(E_k) 
  return(Exp)
}