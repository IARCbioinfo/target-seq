#! /usr/bin/env nextflow
// run using for ex.: nextflow run vcf_annovar.nf --vcf test.vcf

vcf = file( params.vcf )

// convert vcf to annovar format by splitting by samples
process convert2annovar {

	input:
	file vcf

	output:
	file "all_variants.*.avinput" into AVI mode flatten

	script:
	"""
	convert2annovar.pl $vcf -format vcf4 -includeinfo -allsample -outfile all_variants
	"""
}

// annovar annotation
process annovar {

	errorStrategy 'ignore'

  tag { SM }

	input:
	file AVI

	output:
	file "${AVI}.hg19_multianno.txt" into AV

	script:
	out_prefix = AVI.baseName
	(full,SM) = (out_prefix =~ /all_variants\.(.+)/)[0]
	"""
	# test emtpy sample annovar input
	nb_var=\$(wc -l < $AVI)
	if [ \$nb_var -gt 0 ]; then
		# add additional columns to the annovar file (bam file name, sample name, barcode, run number
		sed -i "s/\$/\t${SM}/" $AVI
		table_annovar.pl -nastring NA -buildver hg19 -remove -protocol refGene,knownGene,ensGene,cytoBand,genomicSuperDups,tfbsConsSites,gwasCatalog,avsnp147,popfreq_all_20150413,exac03nontcga,kaviar_20150923,cosmic79,icgc21,clinvar_20150330,mcap,revel,dbnsfp30a,dbnsfp31a_interpro,dbscsnv11 -operation g,g,g,r,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f -otherinfo $AVI '${params.avdb}'
		# put back the columns names from the VCF file that annovar replaces with Otherinfo
		sed -i '1s/Otherinfo/CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGENOTYPE\tSM/' ${AVI}.hg19_multianno.txt
	fi
	"""
}

// merge all annotated files in one big file
process collect_result {

	storeDir { vcf.getParent() }

	input:
	file '*.hg19_multianno.txt' from AV.toList()

	output:
	file "${out_prefix}_annotated.txt" into AVO

	script:
	out_prefix = vcf.baseName
	"""
	# keep the header (column names) only from the first file
	awk 'FNR==1 && NR!=1{next;}{print}' *.hg19_multianno.txt > ${out_prefix}_annotated.txt
	"""
}
