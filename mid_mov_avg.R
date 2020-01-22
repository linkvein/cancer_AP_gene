mid_mov_avg = function(input_table, percent){
  
  ########################################
  # moving avg analysis
  ########################################
  
  # this needs to take in a table of genes(rows) vs. conditions(columns)
  # and produce a table of genes vs. s/i/r
  
  cell_line_vec = c()
  min_point_vec = c()
  max_point_vec = c()
  min_per_vec = c()
  max_per_vec = c()
  
  essential_genes = list()  # this will be a list of all essential genes by cell line
  permissive_genes = list() # " permissive "
  
  for(i in 1:length(input_table[1,])){
    
    temp_data = input_table[,i]
    
    temp_y = temp_data[order(temp_data)]
    temp_x = seq(1,length(temp_y))
    
    running_slope = c()
    for(j in 3:(length(temp_x)-2)){
      
      temp_start = j - 2
      temp_finish = j + 2
      temp_slope = (temp_y[temp_finish] - temp_y[temp_start]) / 4
      
      running_slope[j-2] <- temp_slope
    }
    
    running_slope_avg = c()
    for(k in 5:(length(running_slope)-4)){
      temp_avg = mean(running_slope[(k-4):(k+4)])
      running_slope_avg[k-4] <- temp_avg
    }
    
    #let's look at slope of middle percentage-fraction
    temp_fraction = round(length(temp_x)* ((50 - percent/2)*0.01))
    middle_slope = (temp_y[length(temp_x) - temp_fraction] - temp_y[temp_fraction]) / ((length(temp_x) - temp_fraction) - temp_fraction)
    
    min_point = min(which(running_slope_avg < middle_slope))
    max_point = max(which(running_slope_avg < middle_slope))
    
    total_genes = length(temp_data)
    
    min_percent = min_point/total_genes
    max_percent = max_point/total_genes
    
    # recording vector data
    cell_line_vec[length(cell_line_vec)+1] <- colnames(input_table)[i]
    min_point_vec[length(min_point_vec)+1] <- min_point
    max_point_vec[length(max_point_vec)+1] <- max_point
    min_per_vec[length(min_per_vec)+1] <- min_percent
    max_per_vec[length(max_per_vec)+1] <- max_percent
    
    # recording essential / permissive gene data
    reordered_genes = rownames(input_table)[order(temp_data)]
    temp_table = cbind(reordered_genes, temp_y)
    
    essential_genes[[i]] <- temp_table[1:min_point,1]
    names(essential_genes)[i] <- colnames(input_table)[i]
    permissive_genes[[i]] <- temp_table[max_point:length(temp_table[,1]),1]
    names(permissive_genes)[i] <- colnames(input_table)[i]
    
  }
  
  return(c(min_percent, max_percent))
}