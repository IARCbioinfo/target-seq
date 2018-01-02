setwd(QC3_output_folder)

summary=read.table("bamResult/bamSummary.txt", stringsAsFactors = F, h=T, sep="\t")

#here need to get true sample name from barcodes.
#sample name are identical for both libraries.
#in this example, barcodes are in the form IonXpress_000-RSXXX and true sample name are in the form RSXXX
summary$SM = unlist(lapply(summary$SM, function(x) unlist(strsplit(x,"-"))[2]))

 #the following code order the summary dataframe by sample name, i.e. take lines 2 by 2 form library pairs.
 samples = summary$Sample
 names(samples)=rank(summary$SM,ties.method = "first") # name the sample by its rank
 samples = samples[as.character(sort(as.numeric(names(samples))))] #sort the samples by the rank, so library pairs will be in chain
 row.names(summary)=summary$Sample
 summary=summary[samples,] #to order the dataframe

#here a dataframe containing the median depth for the libraries of each sample is computed.
data_depth = data.frame(cbind("SM"=summary[c(T,F),"SM"],
                                       "median.lib1"=as.numeric(summary[c(T,F),"On.target.Median.Depth."]),
                                       "median.lib2"=as.numeric(summary[c(F,T),"On.target.Median.Depth."])), stringsAsFactors = F)

#finally, the sample with at least one library with a median coverage below "depth_thr" are extracted.
depth_thr = 1000
data_depth[as.numeric(data_depth$median.lib1)<depth_thr | as.numeric(data_depth$median.lib2) < depth_thr,]
