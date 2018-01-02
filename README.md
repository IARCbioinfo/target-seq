# Calling mutations with needlestack on targeted sequencing data

This is the description of the IARC workflow to call somatic mutations with needlestack on targeted sequencing data.  
Usually, IARC data are generating with an IonTorrent Proton sequencer.  
Here samples are sequenced twice to remove library-preparation errors.  
This workflow contains 5 disctinct steps:

## STEP 1: local re-alignment with abra2

[abra2](https://github.com/mozack/abra2) is a NGS realigner that uses localized assembly and global realignment to improve detection of complex variant such as insertion or deletion.  
IARC bioinformatic platform has developed an [abra nextflow pipeline](https://github.com/IARCbioinfo/abra-nf) for a user-friendly usage of the software.

Command line example:

```
nextflow run iarcbioinfo/abra-nf --bed myBedFile.bed --ref genome.fasta --abra_path pathToAbra.jar --single --bam_folder BAM/ --output_folder abra_output/
```

It takes around 30min on one BAM from ~1500 positions in 10,000X.  
The process can be parallelized with the option __--threads__ (_e.g. --threads 6_).

## STEP 2: coverage quality control with QC3

[QC3](https://github.com/slzhao/QC3) is a quality control tools for 3 type of data, sequencing, alignment and variant calling.  
The tool has been optimized by the IARC bioinformatic platform to decrease computation time when used on target-sequencing data which generates highly covered positions.  

Command line example:

```
qc3.pl -m b -i list_bam.txt -o output_folder -r myBedFile.bed -nod -cm 2 -no_batch -d -d_cumul 1000,5000,10000
```

This example computes (amongst other things) the median depth, for each postion in the target bed file for each sample in the _list_bam.txt_ file. Thanks to option __-nod__, untarget depth is not computed at all.  
We recommand to reserve a lot of memory and to use option `-t` to run the tool on multiple threads (_e.g. -t 24_).  
To contruct the _list_bam.txt_ file, a possibility is to run:

```
find /whole_path_to_bam_files/*bam > list_bam.txt
```

Samples are sequenced in technical duplicates, and those with a median coverage less than a particular threshold in at least one of the two libraries need to be remove to avoid technical artifacts.  
[This R script](https://github.com/tdelhomme/target-seq/blob/master/bin/QC3_analysis.r) extracts those samples.
Once the bad samples are identified, the next step is to remove them in the folder used in the variant calling with needlestack (a good practice would be to create a new bad folder and then make symlinks before removing bad samples).  


## STEP 3: variant calling with needlestack

## STEP 4: variant annotation with annovar

## STEP 5: post-filtering on bad samples/positions
