ziwei = function(x,y,z){
  
  if(x < 1 || y < 1){
    return(0)
  }
  
  term1.1 = (y * factorial(x + z))/factorial(x + y + z)
  sum1 = 0
  for(k in 2:(x + z + 1)){
    temp1 = (k * factorial(x + y + z - k))/factorial(x + z - k + 1)
    sum1 = sum1 + temp1
  }
  
  term2.1 = (x * factorial(y + z))/factorial(x + y + z)
  sum2 = 0
  for(k in 2:(y + z + 1)){
    temp2 = (k * factorial(x + y + z - k))/factorial(y + z - k + 1)
    sum2 = sum2 + temp2
  }
  
  term3.1 = ((x + y) * factorial(z))/factorial(x + y + z)
  sum3 = 0
  
  if((z + 1) >= 2){
    for(k in 2:(z + 1)){
      temp3 = (k * factorial(x + y + z - k))/factorial(z - k + 1)
      sum3 = sum3 + temp3
    }
  } else{
    sum3 = 0
  }
  
  solution = term1.1*sum1 + term2.1*sum2 - term3.1*sum3
  
  return(solution)
}