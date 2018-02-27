#####################################################################################################
#' calculate_heterogeneity_scores.py
#' This script computes Intra-Sample Heterogeneity Scores from bisulfite sequencing data. This is
#' the script used for the calculation in the Ewing Sarcoma data set, for which we had raw 
#' sequencing reads in the form of fastq files. The script has two parts: in the first part
#' we process the reads and align them to a reference genome and in the second step we compute the 
#' ISH scores FDRP, qFDRP, PDR, Epipolymorphism, Entropy and MHL.
#####################################################################################################

#!/usr/bin/python

import os
import pypiper
import argparse
import re
import pandas
import sys
import string

#' Part of the looper setup. We add two additional arguments to the parser, one is the sample id of the currently processed sample and the second is the path to the fastq file containing the raw reads for the corresponding sample. These two arguments are passed through  config/pipeline_interface.yaml to map column names in the sample anntotation sheet to the name of the argument here.
parser = argparse.ArgumentParser(description="Pipeline")
parser.add_argument("--sample_id","-o",help="id of sample to be analyzed")
parser.add_argument("--fastq_name",help="path to bam file of sample to be analyzed")
parser = pypiper.add_pypiper_args(parser,groups=["pypiper","looper"])
args = parser.parse_args()

manager = pypiper.PipelineManager(name="HETEROGENEITY",
	outfolder=args.output_parent,
	args=args)

#' We create a folder for each sample individually to avoid problems with lock files
if not os.path.exists(args.output_parent + "/" + args.sample_id):
	os.makedirs(args.output_parent + "/" + args.sample_id)

sample_folder = args.output_parent + "/" + args.sample_id + "/"

#####################################################################################################
#' PART I: Preprocessing
#####################################################################################################

#' Trim reads with Trim Galore. The arguments for TrimGalore! are specified in calculate_heterogeneoty_scores.yaml
trimmed = sample_folder + "trimmed/"
os.environ['PATH'] = manager.config.tools.cutadapt + ":" + os.environ['PATH']
os.system("mkdir "+trimmed)
cmd = " ".join([manager.config.tools.trim_galore, args.fastq_name, "-a", manager.config.parameters.trim_galore.adapater, "-q", str(manager.config.parameters.trim_galore.quality), "-stringency",str(manager.config.parameters.trim_galore.stringency), "-e", str(manager.config.parameters.trim_galore.error_rate), "--length", str(manager.config.parameters.trim_galore.length), "--output_dir", trimmed,"--rrbs"])
manager.run(cmd, lock_name=args.sample_id+"locker")

#' Map reads to the reference genome with bsmap. Arguments are specified in calculate_heterogeneity_scores.yaml
os.environ['PATH'] = manager.config.parameters.bsmap.samtools_path + ":" + os.environ['PATH']
mapped = sample_folder + "mapped/"
mapped_file = mapped+"mapped.bam"
to_analyze = trimmed + "*.fq.gz"
cmd = " ".join(["mkdir", mapped, ";",
manager.config.tools.bsmap,
"-a", to_analyze,
"-d",manager.config.resources.reference_genome,
manager.config.parameters.bsmap.rrbs,
"-w", str(manager.config.parameters.bsmap.best_hits),
"-v", str(manager.config.parameters.bsmap.mismatch),
"-r", str(manager.config.parameters.bsmap.repeat),
"-p", args.cores,
"-n", str(manager.config.parameters.bsmap.map_strand),
"-s", str(manager.config.parameters.bsmap.seed),
"-S", str(manager.config.parameters.bsmap.random_seed),
"-f", str(manager.config.parameters.bsmap.filter),
"-q", str(manager.config.parameters.bsmap.quality_threshold),
"-u",
"-o", mapped_file,
";rm -r", trimmed])
manager.run(cmd,lock_name=args.sample_id+'locker')

#' Sort bam file with samtools. Now further arguments need to be specified.
bam = sample_folder + "bam/"
sorted_bam = bam + "sorted.bam"
cmd = "mkdir " + bam + '; ' + manager.config.tools.samtools + " sort " + mapped + "*.bam > " + sorted_bam + "; rm -r " + mapped
manager.run(cmd,lock_name=args.sample_id+'locker')

#' Index bam file with samtools. Now further arguments need to be specified.
cmd = manager.config.tools.samtools + " index " + sorted_bam
manager.run(cmd,lock_name=args.sample_id+'locker')

#' Since methclone requires the bam file to contain a methylation call string of the form Z..z...h.. as introduced by bismark, we add this string to the mapped bam file with a custom R script. 
bismark_bam = sample_folder + "bismark_bam.bam"
cmd = " ".join([manager.config.tools.rscript, manager.config.parameters.rscript.convert_bam, sorted_bam, sample_folder, "bismark_bam.bam", args.cores, ";rm",sorted_bam])
manager.run(cmd,lock_name=args.sample_id+"locker1")

#####################################################################################################
#' PART II: ISH Score Calculation
#' Those were the preprocessing steps needed to be able to comput the heterogeneity scores. Now
#' the calculation actually takes place for FDRP, qFDRP, PDR, Epipolymorphism, Entropy and MHL
#####################################################################################################

#' Calculate FDRP for the output bam file with the corresponding script. Arguments are specified in calculate_heterogeneity_scores.yaml.
fdrp = sample_folder + "FDRP/"
cmd = " ".join([manager.config.tools.rscript,manager.config.parameters.rscript.fdrp_script,"FDRP",fdrp ,bismark_bam, manager.config.resources.rnb_set,args.cores])
manager.run(cmd,lock_name=args.sample_id+'locker')

#' Calculate qFDRP for the output bam file with the corresponding script. Arguments are specified in calculate_heterogeneity_scores.yaml.
qfdrp = sample_folder + "qFDRP/"
cmd = " ".join([manager.config.tools.rscript,manager.config.parameters.rscript.qfdrp_script,"qFDRP",qfdrp ,bismark_bam, manager.config.resources.rnb_set,args.cores])
manager.run(cmd,lock_name=args.sample_id+'locker')

#' Calculate PDR for the output bam file with the corresponding script. Arguments are specified in calculate_heterogeneity_scores.yaml.
pdr = sample_folder + "PDR/"
cmd = " ".join([manager.config.tools.rscript,manager.config.parameters.rscript.pdr_script,"PDR",pdr ,bismark_bam, manager.config.resources.rnb_set,args.cores])
manager.run(cmd,lock_name=args.sample_id+'locker')

#' Calculate haplotpye information from bismark output with the corresponding script. Arguments are specified in calculate_heterogeneity_scores.yaml. We here use the (updated) scripts from the MHL publication (doi:10.1038/ng.3805).
mhl_folder = sample_folder + "MHL/"
hapinfo_output = mhl_folder + "hapinfo.txt"
roi = sample_folder + "roi.bed"
os.system("mkdir "+mhl_folder)
cmd = " ".join([manager.config.tools.perl, manager.config.parameters.mhl.to_hapinfo, manager.config.resources.roi, bismark_bam, "bismark", manager.config.resources.roi, ">", hapinfo_output])
manager.run(cmd,lock_name=args.sample_id+'locker')

#' Calculate MHL from haplotype information also with the scripts from the publication of the MHL score.
mhl_output = mhl_folder + "mhl.txt"
cmd = "echo " + hapinfo_output + " > " + mhl_folder + "info_list.txt; " + " ".join([manager.config.tools.perl, manager.config.parameters.mhl.hapinfo_to_mhl, mhl_folder + "info_list.txt", ">", mhl_output])
manager.run(cmd,mhl_output)

#' Run methclone (10.1186/s13059-014-0472-5) software with arguments specified in calculate_heterogeneity_scores.yaml
epiallele_output = sample_folder
methclone_output = epiallele_output + "methclone_tmp.txt.gz"
cmd = " ".join([manager.config.tools.methclone, bismark_bam, bismark_bam, methclone_output, "methclone", str(manager.config.parameters.methclone.meth_diff), str(manager.config.parameters.methclone.distance_cutoff), str(manager.config.parameters.methclone.coverage)])
manager.run(cmd,methclone_output)

#' Calculate epipolymorphism from methclone's output with custom R scripts.
epipoly_output = epiallele_output + "epipoly.csv"
cmd = " ".join([manager.config.tools.rscript, manager.config.parameters.rscript.conversion_epipoly, methclone_output, "epipoly", epiallele_output])
manager.run(cmd,epipoly_output)

#' Calculate methylation entropy from methclone's output with custom R scripts.
entropy_output = epiallele_output + "entropy.csv"
cmd = " ".join([manager.config.tools.rscript, manager.config.parameters.rscript.conversion_entropy, methclone_output, "entropy", epiallele_output])
manager.run(cmd,entropy_output)

#' Cleanup
os.system("rm " + methclone_output)
os.system("rm -r " + bam)
os.system("rm " + bismark_bam + "*")

#" Stop pipeline
manager.stop_pipeline()
