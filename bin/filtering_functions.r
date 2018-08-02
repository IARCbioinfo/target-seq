# function to compute, given a variant (dataframe + row number), the minimum distance to another one
min_distance <-function(all_my_data, row_number_variant, max_af=0.1){
  cur_data = all_my_data[-row_number_variant,] #remove current variant from distances computations
  cur_data = cur_data[which(cur_data$Chr==all_my_data[row_number_variant,"Chr"] &
                              cur_data$old_SM==all_my_data[row_number_variant,"old_SM"] &
                              cur_data$AF>=max_af) ,]
  indels = which(as.numeric(cur_data$End) > as.numeric(cur_data$Start))
  min( unlist(lapply(all_my_data[row_number_variant,"Start"]:all_my_data[row_number_variant,"End"], function(var_pos) {
    #create a vector of all position where we observed a variant
    all_pos_sort = sort(as.numeric(unique(c(cur_data$Start, unlist(lapply(indels, function(i) cur_data[i,"Start"]:cur_data[i,"End"]))))))
    min(var_pos - all_pos_sort[which(all_pos_sort<=var_pos)[length(which(all_pos_sort<=var_pos))]] ,
        all_pos_sort[which(all_pos_sort>=var_pos)[1]] - var_pos, na.rm = T) #na if only one variant compared
  })) ) #min to minimum distance over all bases of the variant
}

# function giving pvalue from a mutation barcode, a data frame and a number of CIR
pvalue_from_mut <- function(mut, data, CIR){
  mut_data = data[data$mutation_bc==mut,]
  nb_var = dim(mut_data)[1]
  nb_pairs = sum(table(mut_data$SM)==2)
  nb_indiv = as.numeric(nb_indivs[as.character(CIR)])/2
  pval = pvalue_pair(p = nb_pairs, N = nb_indiv, k = nb_var)
  return(pval)
}

# function to compute pvalue from Matthieu's formula for each row of the data matrix (unique pvalue by mutation_bc but duplicated by rows)
get_pvalue_pairs <- function(data, CIR){
  data$mutation_bc=paste(data$Chr,data$Start,data$End,data$Ref,data$Alt,sep="_")
  mut_bc=unique(paste(data$Chr,data$Start,data$End,data$Ref,data$Alt,sep="_"))
  pvals = unlist(lapply(mut_bc, pvalue_from_mut, data, CIR)) #create a vector of pvalues for each mutation_bc (names)
  names(pvals) = mut_bc
  unlist(lapply(1:dim(data)[1], function(i, data){ # for each line of the data frame return the pvalue of the mutation
    pvals[paste(data[i,"Chr"],data[i,"Start"],data[i,"End"],data[i,"Ref"],data[i,"Alt"],sep="_")]
  }, data))
}

# function to get minimum RVSB from paired mutations (in two libraries sharing a same SM)

get_minRVSB <- function(data){
  bc_data = paste(data$SM, data$Chr ,data$Start, data$Ref, data$Alt)
  lib1 = which(duplicated(bc_data))
  lib2 = unlist(lapply(lib1, function(i){
    ids = which(bc_data[i] == bc_data)
    ids[which(ids!=i)]
  }))
  data[lib1,"minRVSB"] = data[lib2,"minRVSB"] = pmin(data[lib1,"RVSB"], data[lib2,"RVSB"])
  data$minRVSB
}
