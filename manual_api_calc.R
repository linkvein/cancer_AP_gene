# version of API_calculator where thresholds are manually defined

manual_api_calc = function(input_table, input_min, input_max){
  
  ########################################
  # this needs to take in a table of genes(rows) vs. conditions(columns)
  # and produce a table of genes vs. s/i/r
  ########################################
  
  essential_genes = list()  # this will be a list of all essential genes by cell line
  permissive_genes = list() # " permissive "
  
  for(i in 1:length(input_table[1,])){
    
    temp_data = input_table[,i]
    
    temp_y = temp_data[order(temp_data)]
    temp_x = seq(1,length(temp_y))
    
    total_genes = length(temp_data)
    
    min_point = round(input_min * 0.01 * total_genes)
    max_point = round(input_max * 0.01 * total_genes)
    
    # recording essential / permissive gene data
    reordered_genes = rownames(input_table)[order(temp_data)]
    temp_table = cbind(reordered_genes, temp_y)
    
    essential_genes[[i]] <- temp_table[1:min_point,1]
    names(essential_genes)[i] <- colnames(input_table)[i]
    permissive_genes[[i]] <- temp_table[max_point:length(temp_table[,1]),1]
    names(permissive_genes)[i] <- colnames(input_table)[i]
    
  }
  
  ########################################
  # calculate the APS for each genes
  ########################################
  # goal is to identify the number of s/i/r cell lines for each gene 
  # then calculate the overall APS and lineage-specific APS scores
  
  # 1. go through each of the conditions and check each of the genes for ess or per
  # 2. annotate any genes that aren't ess or per as in(ert)
  # 3. calculate APS for each gene
  
  # initializing matrix to store s/i/r designation
  APS_matrix = matrix(, nrow = length(input_table[,1]), ncol = length(input_table[1,]))
  colnames(APS_matrix) <- colnames(input_table)
  rownames(APS_matrix) <- rownames(input_table)
  
  # populating s/i/r matrix
  for(i in 1:length(input_table[1,])){
    temp_line = colnames(input_table)[i]
    
    for(j in 1:length(input_table[,1])){
      temp_gene = rownames(input_table)[j]
      index = "NA"
      
      if(temp_gene %in% essential_genes[[i]]){index = "s"
      } else if(temp_gene %in% permissive_genes[[i]]){index = "r"
      } else{index = "i"}
      
      APS_matrix[j,i] <- index
    }
    # print(i)
  }
  
  APS_vector = c()
  # calculating ziwei score for each of the genes
  for(i in 1:length(APS_matrix[,1])){
    temp_ziwei = ziwei(x = length(which(APS_matrix[i,] == "r")),
                       y = length(which(APS_matrix[i,] == "s")),
                       z = length(which(APS_matrix[i,] == "i")))
    
    APS_vector[i] <- temp_ziwei
    
    # print(i)
  }
  names(APS_vector) <- rownames(APS_matrix)
  
  APS_vector[which(APS_vector == 0)] <- NaN
  
  return(APS_vector)
}