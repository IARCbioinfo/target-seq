library(xlsx)
source("https://gist.githubusercontent.com/mfoll/a4dfbb92068dc559f130/raw/714dc8c2e97987fd4385dcef2722b3ef986d38d6/get_vcf_data.r")
source("filtering_functions.r")
source("prob_pairs.r")
source("nb_muts.r")

interesting_cols=c("SM","old_SM","Chr","Start","End","Ref","Alt","Gene.refGene","TYPE","Func.refGene","GeneDetail.refGene",
                   "ExonicFunc.refGene","AAChange.refGene","AAChange.refGene_NM_546","PopFreqMax","ESP6500siv2_ALL",
                   "X1000G_ALL","ExAC_ALL","ExAC_nontcga_ALL","avsnp147","SIFT_pred","Polyphen2_HDIV_pred",
                   "clinvar_20150330","cosmic79","MCAP","REVEL","RVSB","QVAL","AO","DP","AF","CONT",
                   "pvalue_pairs50") # ,"AAChange.refGene_NM_546"

nb_indivs = 30 #number of total libraries (wc -l *.bam)

# reading of the output table from annovar
data_annotated=read.table("variants_annotated.txt",quote="\"",stringsAsFactors=F,sep="\t",header=T)
# modifying sample name (SM), to have same SM for both libraries
data_annotated$old_SM=data_annotated$SM
data_annotated$SM=gsub(".*-(.*)$","\\1",data_annotated$SM)

# reformat some annovar annotations
data_annotated$ExonicFunc.refGene[which(data_annotated$Func.refGene=="splicing")]="splicing"
data_annotated$AAChange.refGene[which(data_annotated$Func.refGene=="splicing")]=data_annotated$GeneDetail.refGene[which(data_annotated$Func.refGene=="splicing")]
NM_546=lapply(strsplit(data_annotated$AAChange.refGene,","),function(test){gsub("TP53:|NM_000546:","",test[grep("NM_000546",test)])})
NM_546[sapply(NM_546,length)==0]<-NA
data_annotated$AAChange.refGene_NM_546=unlist(NM_546)

# get usefull VCF fields
data_annotated$AF=get_genotype(data_annotated$GENOTYPE,data_annotated$FORMAT[1],"AF")
data_annotated$AO=get_genotype(data_annotated$GENOTYPE,data_annotated$FORMAT[1],"AO")
data_annotated$DP=get_genotype(data_annotated$GENOTYPE,data_annotated$FORMAT[1],"DP")
data_annotated$RVSB=get_genotype(data_annotated$GENOTYPE,data_annotated$FORMAT[1],"RVSB")
data_annotated$QVAL=get_genotype(data_annotated$GENOTYPE,data_annotated$FORMAT[1],"QVAL")
data_annotated$TYPE=get_info(data_annotated$INFO,"TYPE",num=F)
data_annotated$ERR=get_info(data_annotated$INFO,"ERR")
data_annotated$CONT=get_info(data_annotated$INFO,"CONT",num=F)

# filter on QVAL if necessary
data_annotated=data_annotated[data_annotated$QVAL>=50,]

# compute for each mutation, the probability that the position is randomly noisy by errors
# we compute this value by counting, amoung the number of positive libraries (QVAL>50), the number of pairs
# and we compare this to what is expected by a random sampling
data_annotated$pvalue_pairs50 = NA
data_annotated$pvalue_pairs50[which(data_annotated$QVAL>=50)] = get_pvalue_pairs(data = data_annotated[which(data_annotated$QVAL>=50),], CIR = run)

# filtering on columns
data_annotated=data_annotated[,interesting_cols]

#number of mutations without any filters, before and after library dup check
dat = data.frame("SM"=names(table(data_annotated$old_SM)),
                   "mutations"=as.numeric(table(data_annotated$old_SM)))
#compute the threshold of number of mutations per library, and the corresponding samples to exclude
thr = nb_mut_thr(dat$mutations)
excluded_libraries = as.character(dat[which(dat$mutations>=thr),"SM"])
excluded_samples = unlist(lapply(excluded_libraries, function(x) unlist(strsplit(x,"-"))[2]))
data_annotated = data_annotated[which(!data_annotated$SM %in% excluded_samples),]

# compute mutations barcodes, that will be used to control for presence in the two libraries
# mutations barcodes are in format: SM_CHR_START_END_REF_ALT
mutation_sample_bc=paste(data_annotated$SM,data_annotated$Chr,data_annotated$Start,data_annotated$End,data_annotated$Ref,data_annotated$Alt,sep="_")
# keep a mutation only if its barcode is present twice
data_annotated=data_annotated[which(duplicated(mutation_sample_bc) | duplicated(mutation_sample_bc, fromLast = TRUE) ),]

#count the minimum distance from another variant, to remove alignment errors
data_annotated$MIN_DIST = unlist(lapply(1:nrow(data_annotated), function(i) min_distance(data_annotated, i) ))
#filter
data_annotated = data_annotated[which(data_annotated$MIN_DIST>5 | is.na(data_annotated$MIN_DIST)),]

write.xlsx(data_annotated,file="all_data_annotated_filtered.xlsx",col.names=T,row.names=F)
