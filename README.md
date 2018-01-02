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
The process can be parallelized with the option __--threads__.

## STEP 2: coverage quality control with QC3

## STEP 3: variant calling with needlestack

## STEP 4: variant annotation with annovar

## STEP 5: post-filtering on bad samples/positions
