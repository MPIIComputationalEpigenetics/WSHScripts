#!/usr/bin/python

#####################################################################################################
#' calculate_heterogeneity_scores.py
#' This script computes Intra-Sample Heterogeneity Scores from bisulfite sequencing data. This is
#' the script used for the calculation in the Kiel Cohort data set, for which we had sequencing
#' reads already mapped to the humane reference genome hg38 with bsmap. The script has two parts:
#' in the first part convert bsmaps output to bismark like and in the second step we compute the 
#' ISH scores FDRP, qFDRP, PDR, Epipolymorphism, Entropy and MHL.
#####################################################################################################

import os
import pypiper
import argparse
import re
import pandas
import sys
import string

#' Part of the looper setup. We add two additional arguments to the parser, one is the sample id of the currently processed sample and the second is the path to the bam file containing the mapped reads (preferably with bsmap). These two arguments are passed through  config/pipeline_interface.yaml to map column names in the sample anntotation sheet to the name of the argument here.

parser = argparse.ArgumentParser(description="Pipeline")
parser.add_argument("--sample_id","-o",help="id of sample to be analyzed")
parser.add_argument("--bam_name",help="path to bam file of sample to be analyzed")
parser = pypiper.add_pypiper_args(parser,groups=["pypiper","looper"])
args = parser.parse_args()

manager = pypiper.PipelineManager(name="HETEROGENEITY",
	outfolder=args.output_parent,
	args=args)

#####################################################################################################
#' PART I: Preprocessing
#####################################################################################################
if not os.path.exists(args.output_parent + "/" + args.sample_id):
	os.makedirs(args.output_parent + "/" + args.sample_id)

#' Run samtools to create fastq from bam file
#fastq = args.output_parent + "/" + args.sample_id  + "/" + args.sample_id + '.fastq'
#cmd = manager.config.tools.samtools + ' fastq ' + ' ' + args.bam_name + ' > ' + fastq
#manager.run(cmd,fastq)

# Map the simulated fastq file to the reference with bimark
#mapped = args.output_parent + "/" + args.sample_id + "/mapped/"
#manager.run("mkdir "+mapped,mapped)
#cmd = manager.config.tools.bismark + " -o " + mapped + " " + manager.config.resources.bisulfite_folder + " --" + manager.config.parameters.bismark.bowtie + " --path_to_bowtie " + manager.config.parameters.bismark.bowtie_path + " --samtools_path " + manager.config.parameters.bismark.samtools + " " + fastq + " --multicore " + str(manager.config.parameters.bismark.cores)
#manager.run(cmd,lock_name=args.sample_id+"locker")

# Map the simulated fastq file to the reference with bsmap
#mapped = args.output_parent + "/" + args.sample_id + "/mapped/"
#manager.run("mkdir "+mapped,mapped)
#cmd = manager.config.tools.bismark + " -o " + mapped + "output.bam -d " + manager.config.resources.bisulfite_folder + " -a " + fastq + " -D C-CGG -w 100 -v 0.08 -r 1 -p 4 -n 0 -s 12 -S 0 -f 5 -q 0 -u -V 2"
#manager.run(cmd,lock_name=args.sample_id+"locker")

# Create the cov files for investigation of methylation states
#covs = args.output_parent + "/" + args.sample_id + '/covs/'
#cmd = 'mkdir ' + covs + '; ' + manager.config.tools.bismark_extractor + ' -o ' + covs + ' --samtools_path ' + manager.config.parameters.bismark.samtools + ' --single-end ' + ' --bedGraph ' + mapped + "*.bam --multicore " + str(manager.config.parameters.bismark.cores)
#manager.run(cmd,covs)

# Sort bam file
#bam = args.output_parent + "/" + args.sample_id + "/bam/"
#sorted_bam = bam + "sorted.bam"
#cmd = "mkdir " + bam + '; ' + manager.config.tools.samtools + " sort " + mapped + "*.bam > " + sorted_bam + "; rm -r " + mapped + " ; rm " + fastq
#manager.run(cmd,lock_name=args.sample_id+'locker')

# Index bam file
#cmd = manager.config.tools.samtools + " index " + sorted_bam
#manager.run(cmd,lock_name=args.sample_id+'locker')
sample_folder = args.output_parent + "/" + args.sample_id + "/"
#' Use script to convert bsmap style bam file to bismark style
bismark_bam = sample_folder + "bismark_bam.bam"
cmd = " ".join([manager.config.tools.rscript, manager.config.parameters.rscript.convert_bam, args.bam_name, sample_folder, "bismark_bam.bam", str(manager.config.parameters.rscript.cores)])
manager.run(cmd,lock_name=args.sample_id+"locker1")

#####################################################################################################
#' PART II: ISH Score Calculation
#' Those were the preprocessing steps needed to be able to comput the heterogeneity scores. Now
#' the calculation actually takes place for FDRP, qFDRP, PDR, Epipolymorphism, Entropy and MHL
#####################################################################################################

#' Calculate FDRP for the output bam file with the corresponding script. Arguments are specified in calculate_heterogeneity_scores.yaml.
fdrp = sample_folder + "FDRP/"
cmd = manager.config.tools.rscript + ' ' + manager.config.parameters.rscript.fdrp_script + ' FDRP ' + fdrp + " " + bismark_bam + ' ' + manager.config.resources.rnb_set + " " + str(manager.config.parameters.rscript.cores)
manager.run(cmd,lock_name=args.sample_id+'locker')

#' Calculate qFDRP for the output bam file with the corresponding script. Arguments are specified in calculate_heterogeneity_scores.yaml.
qfdrp = sample_folder + "qFDRP/"
cmd = manager.config.tools.rscript + ' ' + manager.config.parameters.rscript.qfdrp_script + ' qFDRP ' + qfdrp + " " + bismark_bam + ' '+ manager.config.resources.rnb_set + " " + str(manager.config.parameters.rscript.cores)
manager.run(cmd,lock_name=args.sample_id+'locker')

#' Calculate PDR for the output bam file with the corresponding script. Arguments are specified in calculate_heterogeneity_scores.yaml.
pdr = sample_folder + "PDR/"
cmd = manager.config.tools.rscript + ' ' + manager.config.parameters.rscript.pdr_script + ' PDR ' + pdr + " " + bismark_bam + ' ' + manager.config.resources.rnb_set + " " + str(manager.config.parameters.rscript.cores)
manager.run(cmd,lock_name=args.sample_id+'locker')

#' Calculate haplotpye information from bismark output with the corresponding script. Arguments are specified in calculate_heterogeneity_scores.yaml. We here use the (updated) scripts from the MHL publication (doi:10.1038/ng.3805).
hapinfo_output = sample_folder + "/hapinfo.txt"
os.environ['PATH'] = manager.config.tools.samtools + ':' + os.environ['PATH']
cmd = " ".join([manager.config.tools.perl, manager.config.parameters.mhl.to_hapinfo, manager.config.resources.roi, bismark_bam, "bismark", manager.config.resources.roi, '>', hapinfo_output])
manager.run(cmd,hapinfo_output)

#' Calculate MHL from haplotype information
mhl_output = sample_folder + "/mhl.txt"
cmd = 'echo ' + hapinfo_output + " > " + sample_folder + "info_list.txt; " + " ".join([manager.config.tools.perl, manager.config.parameters.mhl.hapinfo_to_mhl,sample_folder + "info_list.txt", '>', mhl_output])
manager.run(cmd,mhl_output)

#' Run methclone (10.1186/s13059-014-0472-5) software with arguments specified in calculate_heterogeneity_scores.yaml
methclone_output = sample_folder + "methclone_tmp.txt.gz"
cmd = " ".join([manager.config.tools.methclone, bismark_bam, bismark_bam, methclone_output, "methclone", str(manager.config.parameters.methclone.meth_diff), str(manager.config.parameters.methclone.distance_cutoff), str(manager.config.parameters.methclone.coverage)])
manager.run(cmd,methclone_output)

#' Calculate epipolymorphism from methclone's output with custom R scripts.
epipoly_output = sample_folder + "epipoly.csv"
cmd = " ".join([manager.config.tools.rscript, manager.config.parameters.rscript.conversion_epipoly, methclone_output, 'epipoly', sample_folder])
manager.run(cmd,epipoly_output)

#' Calculate methylation entropy from methclone's output with custom R scripts.
entropy_output = args.output_parent + "/" + args.sample_id + 'entropy.csv'
cmd = " ".join([manager.config.tools.rscript, manager.config.parameters.rscript.conversion_entropy, methclone_output, 'entropy', sample_folder]) + ' ; rm ' + methclone_output + "; rm " + bismark_bam
manager.run(cmd,entropy_output)

#' Cleanup
os.system("rm " + methclone_output)
os.system("rm " + bismark_bam + "*")
os.system("rm " + hapinfo_output)
os.system("rm " + sample_folder + "info_list.txt")

#" Stop pipeline
manager.stop_pipeline()
