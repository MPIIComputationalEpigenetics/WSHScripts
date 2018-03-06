#!/TL/deep/projects/work/mage/tools/Python-2.7.13/python

#' This is a python script responsible for calculating MHL and Epipolymorphism scores for reads 
#' processed by any read alignment or methylation calling tool. It uses the power of the pipeline
#' management system pypiper to provide several neccessary functions of a pipeline.

import os
import pypiper
import argparse

parser = argparse.ArgumentParser(description="Pipeline")
parser.add_argument('--bam_folder','-b',help='folder containing input bam files')
parser.add_argument('--output_folder','-o',help='folder to store the output files')
parser.add_argument('--epipoly_folder',help='folder to store the output files')
parser.add_argument('--mhl_folder',help='folder to store the output files')
parser.add_argument('--sample_name','-s',help='name of the sample as in bam_folder')
parser.add_argument('--output_name','-n',help='name of the output file')
parser = pypiper.add_pypiper_args(parser,groups=["pypiper","looper"])
args = parser.parse_args()

manager = pypiper.PipelineManager(name="HETEROGENEITY",
	outfolder=args.output_folder,
	args=args)

#' Run samtools to create fastq from bam file
fastq = args.output_folder + args.output_name + '.fastq'
cmd = manager.config.tools.samtools + ' fastq ' + ' ' +args.bam_folder + args.sample_name + '.bam' + ' > ' + fastq
manager.run(cmd,fastq)

#' Run bismarck to realign the reads
bismark_bam = args.output_folder + args.output_name + '_bismark.bam'
cmd = manager.config.tools.bismark + ' -o ' + args.output_folder + ' ' +manager.config.resources.reference + ' --' + manager.config.parameters.bismark.bowtie + ' --path_to_bowtie ' + manager.config.parameters.bismark.bowtie_path + ' --samtools_path ' + manager.config.parameters.bismark.samtools + ' --multicore ' + args.cores + ' ' + fastq
manager.run(cmd,bismark_bam)

#' Sort bam file
sorted_bam = args.output_folder + args.output_name + '_bismark_sorted' + '.bam'
cmd = manager.config.tools.samtools + ' sort ' + bismark_bam + ' > ' + sorted_bam + ';  ' + manager.config.tools.samtools + ' index ' + sorted_bam
manager.run(cmd,sorted_bam)

#' Run methclone software
methclone_output = args.output_folder + args.output_name + '_tmp.txt.gz'
cmd = " ".join([manager.config.tools.methclone, sorted_bam, sorted_bam, methclone_output, args.output_name, str(manager.config.parameters.methclone.meth_diff), str(manager.config.parameters.methclone.distance_cutoff), str(manager.config.parameters.methclone.coverage)])
manager.run(cmd,methclone_output)

#' Calculate epipolymorphism from methclone's output and remove this
epipoly_output = args.epipoly_folder + args.output_name + '.csv'
cmd = " ".join([manager.config.tools.rscript, manager.config.parameters.rscript.conversion, methclone_output, args.output_name, args.epipoly_folder]) + ' ; rm ' + methclone_output
manager.run(cmd,epipoly_output)

#' Calculate haplotpye information from bismark output
hapinfo_output = args.mhl_folder + args.output_name + '_hapinfo.txt'
os.environ['PATH'] = manager.config.parameters.bismark.samtools + ':' + os.environ['PATH']
cmd = " ".join([manager.config.tools.perl, manager.config.parameters.mhl.to_hapinfo, manager.config.resources.roi, sorted_bam, manager.config.tools.bismark, manager.config.resources.chrom_sizes, manager.config.resources.cpg_annotation, '>', hapinfo_output])
manager.run(cmd,hapinfo_output)

#' Calculate MHL from haplotype information
mhl_output = args.mhl_folder + args.output_name + '.txt'
cmd = 'echo ' + hapinfo_output + " > " + args.output_folder + "info_list.txt; " + " ".join([manager.config.tools.perl, manager.config.parameters.mhl.hapinfo_to_mhl, args.output_folder + "info_list.txt", '>', mhl_output])
manager.run(cmd,mhl_output)

#' Cleanup
cmd = 'rm ' + hapinfo_output + '; rm ' + sorted_bam + '; rm ' + sorted_bam + '.bai' + '; rm ' + bismark_bam + '; rm ' + fastq + '; rm ' + methclone_output + '; rm ' + args.output_folder + 'info_list.txt' + '; rm ' + args.output_folder + args.output_name + '_bismark_SE_report.txt'
manager.run(cmd,lock_name=args.output_name)

#' Stop pipeline
manager.stop_pipeline()
