#!/usr/bin/python

#####################################################################################################
#' simulate_erosion.py
#' This script runs the pipeline to create regions showing DNA methylation erosion. It is based on
#' two parameters: alpha determines the probability of creating a methylated CpG at any position and
#' gamma replicates those reads gamma-times to simulate stochasticity of read selection.
#####################################################################################################

import os
import pypiper
import argparse
import re
import pandas
import sys
import random

parser = argparse.ArgumentParser(description="Pipeline")
parser.add_argument("--sample_id","-o",help="id of the analyzed region")
parser.add_argument("--chr",help="chromosome on which the specified region lies")
parser.add_argument("--start",help="start position of specified region")
parser.add_argument("--end",help="end position of specified region")
parser = pypiper.add_pypiper_args(parser,groups=["pypiper","looper"])
args = parser.parse_args()

manager = pypiper.PipelineManager(name="SIMULATION",
	outfolder=args.output_parent,
	args=args)

if not os.path.exists(args.output_parent + "/" + args.sample_id):
	os.makedirs(args.output_parent + "/" + args.sample_id)

pipe_folder = os.path.dirname(sys.argv[0])  + "/"

# Select the region we want to analyze
files = os.listdir(manager.config.resources.genome_folder)
args.chr = str(args.chr)
sel_file = [s for s in files if (args.chr+".fa") == s]
if(len(sel_file) > 0):
	sel_file = sel_file[0]
else:
	sys.exit("Specified chromosome " + args.chr + " not in the reference genome")

start_line = int(round(int(args.start)/50,0))
end_line = int(round(int(args.end)/50,0))
count = 0
result_list = list()
for line in open(manager.config.resources.genome_folder+sel_file):
	count += 1
	if count >= start_line and count <= end_line:
		result_list.append(line)
		continue
	if count > end_line:
		break
	if line[1] == "N" and count < end_line:
		start_line += 1
		end_line += 1
		result_list = list()

genome_location = args.output_parent + "/" + args.sample_id + "/simulated.fa"
out_file = open(genome_location,'wb')
out_file.write(">" + args.chr + "\n" + ''.join(result_list))
out_file.close()

erosion_start = random.randint(2,len(result_list))
erosion_end = random.randint(erosion_start,len(result_list)-1)
while erosion_end - erosion_start < 3:
	erosion_start = random.randint(1,len(result_list))
	erosion_end = random.randint(erosion_start,len(result_list))

first_region = result_list[1:erosion_start]
erosion = result_list[erosion_start+1:erosion_end-1]
second_region = result_list[erosion_end:count]
first_region = ">" + args.chr + "\n" + ''.join(first_region)
erosion = ">" + args.chr + "\n" + ''.join(erosion)
second_region = ">" + args.chr + "\n" + ''.join(second_region)
genome_location = args.output_parent + "/" + args.sample_id + '/'
if not os.path.exists(genome_location+'first/'):
	os.makedirs(genome_location+'first/')

out_file = open(genome_location+'first/first_region.fa','wb')
out_file.write(first_region)
out_file.close()
if not os.path.exists(genome_location+'erosion/'):
	os.makedirs(genome_location+'erosion/')

out_file = open(genome_location+'/erosion/erosion.fa','wb')
out_file.write(erosion)
out_file.close()
if not os.path.exists(genome_location+'second/'):
	os.makedirs(genome_location+'second/')

out_file = open(genome_location+'/second/second_region.fa','wb')
out_file.write(second_region)
out_file.close()
cell_meth = random.randint(0,100)

cell_meth = random.randint(1,100)
erosion_score = 100-abs(50-cell_meth)
is_negative = random.randint(0,100)>50
num_repeat = random.randint(2,10)
out_file = open(args.output_parent + "/" + args.sample_id +'/erosion_score.txt','wb')
out_string = str(erosion_score) + "\n" + str(num_repeat) + "\n" + str(erosion_start*50) + "\n" + str(erosion_end*50)
if is_negative:
	out_string = "Negative control\n" + out_string
else:
	out_string = "Positive control\n" + out_string

out_file.write(out_string)
out_file.close()

background_meth = 0

if is_negative:
	cell_meth = background_meth

num_reads = manager.config.parameters.sherman.number/num_repeat

frac_first = float(erosion_start)/float(len(result_list))
frac_second = float(len(result_list)-erosion_end)/float(len(result_list))
frac_dmr = float(erosion_end-erosion_start)/float(len(result_list))

# Run Sherman to simulate reads
simulated = args.output_parent + "/" + args.sample_id + "/simulated/"
manager.run("mkdir " + simulated,simulated)
cmd = manager.config.tools.sherman + " --length " + str(manager.config.parameters.sherman.length) + " --number_of_seqs " + str(int(round(manager.config.parameters.sherman.number*frac_second))) + " --CG_conversion " + str(background_meth) + " --CH_conversion " + str(manager.config.parameters.sherman.CH) + " --genome_folder " + args.output_parent + "/" + args.sample_id + "/second/ --outfolder " + simulated + " --outname second " + " --error_rate " + str(manager.config.parameters.sherman.error_rate) + " --quality " + str(manager.config.parameters.sherman.quality) + " --rrbs"
manager.run(cmd,lock_name=args.sample_id+'locker')

# Run Sherman to simulate reads
simulated = args.output_parent + "/" + args.sample_id + "/simulated/"
cmd = manager.config.tools.sherman + " --length " + str(manager.config.parameters.sherman.length) + " --number_of_seqs " + str(int(round(manager.config.parameters.sherman.number*frac_first))) + " --CG_conversion " + str(background_meth) + " --CH_conversion " + str(manager.config.parameters.sherman.CH) + " --genome_folder " + args.output_parent + "/" + args.sample_id + "/first/ --outfolder " + simulated + " --outname first " + " --error_rate " + str(manager.config.parameters.sherman.error_rate) + " --quality " + str(manager.config.parameters.sherman.quality) + " --rrbs"
manager.run(cmd,lock_name=args.sample_id+'locker')

simulated = args.output_parent + "/" + args.sample_id + "/simulated/"
cmd = manager.config.tools.sherman + " --length " + str(manager.config.parameters.sherman.length) + " --number_of_seqs " + str(int(round(num_reads*frac_dmr))) + " --CG_conversion " + str(cell_meth) + " --CH_conversion " + str(manager.config.parameters.sherman.CH) + " --genome_folder " + args.output_parent + "/" + args.sample_id + "/erosion/ --outfolder " + simulated + " --outname erosion " + " --error_rate " + str(manager.config.parameters.sherman.error_rate) + " --quality " + str(manager.config.parameters.sherman.quality) + " --rrbs"
manager.run(cmd,lock_name=args.sample_id+'locker')

for run in range(1,num_repeat+1):
	# Merge the fastq files
	cmd = "cat " + simulated + "erosion.fastq >> " + simulated + "all_erosion.fastq"
	manager.run(cmd,lock_name=args.sample_id+'locker')

os.remove(simulated + "erosion.fastq")

# Merge the fastq files
cmd = "cat " + simulated + "*.fastq > " + simulated + "/merged.fq"
manager.run(cmd,lock_name=args.sample_id+'locker')

# Map the simulated fastq file to the reference with bimark
mapped = args.output_parent + "/" + args.sample_id + "/mapped/"
cmd = "mkdir " + mapped + "; cd " + mapped + "; " + manager.config.tools.bismark + " -o " + mapped + " " + manager.config.resources.bisulfite_folder + " --" + manager.config.parameters.bismark.bowtie + " --path_to_bowtie " + manager.config.parameters.bismark.bowtie_path + " --samtools_path " + manager.config.parameters.bismark.samtools + " " + simulated + "*.fq"
manager.run(cmd,mapped)

# Create the cov files for investigation of methylation states
covs = args.output_parent + "/" + args.sample_id + '/covs/'
cmd = 'mkdir ' + covs + '; ' + manager.config.tools.bismark_extractor + ' -o ' + covs + ' --samtools_path ' + manager.config.parameters.bismark.samtools + ' --single-end ' + ' --bedGraph ' + mapped + '*.sam'
manager.run(cmd,covs)

# Create and store rnbSet
cmd = 'rm ' + covs + '*.txt; rm ' + covs + '*.png; rm ' + covs + '*.bedGraph; '+ 'echo "sample_id,filename\nTest,merged.fq_bismark_bt2.bismark.cov" > ' + args.output_parent + "/" + args.sample_id + '/sample_annotation.csv; ' +  manager.config.tools.rscript + ' ' + pipe_folder + manager.config.parameters.rscript.create_rnb + ' ' + args.output_parent + "/" + args.sample_id + "; rm " + args.output_parent + "/" + args.sample_id + "/sample_annotation.csv"
manager.run(cmd,lock_name=args.sample_id+'locker')

# Write sam as bam file
bam = args.output_parent + "/" + args.sample_id + "/bam/"
cmd = "mkdir " + bam + '; ' +manager.config.tools.samtools + " view -b " + mapped + "merged.fq_bismark_bt2.sam > " + bam + "simulated.bam; rm -r " + mapped
manager.run(cmd,bam)

# Sort bam file
sorted_bam = bam + "sorted.bam"
cmd = manager.config.tools.samtools + " sort " + bam + "* > " + sorted_bam + "; rm " + bam + "simulated.bam; rm -Rf " + simulated + "; rm -Rf " + args.output_parent + "/" + args.sample_id + "/first; rm -Rf " + args.output_parent + "/" + args.sample_id + "/second"
manager.run(cmd,lock_name=args.sample_id+'locker')

# Index bam file
cmd = manager.config.tools.samtools + " index " + sorted_bam
manager.run(cmd,lock_name=args.sample_id+'locker')

# Calculate FDRP for the output bam file
fdrp = args.output_parent + "/" + args.sample_id + '/FDRP/'
cmd = manager.config.tools.rscript + ' ' + pipe_folder + manager.config.parameters.rscript.fdrp_script + ' FDRP ' + fdrp + " " + sorted_bam + ' ' + args.output_parent + "/" + args.sample_id + '/rnbSet.zip ' + str(manager.config.parameters.rscript.cores)
manager.run(cmd,lock_name=args.sample_id+'locker')

# Calculate qFDRP for the output bam file
qfdrp = args.output_parent + "/" + args.sample_id + '/qFDRP/'
cmd = manager.config.tools.rscript + ' ' + pipe_folder + manager.config.parameters.rscript.qfdrp_script + ' qFDRP ' + qfdrp + " " + sorted_bam + ' '+ args.output_parent + "/" + args.sample_id + '/rnbSet.zip ' + str(manager.config.parameters.rscript.cores)
manager.run(cmd,lock_name=args.sample_id+'locker')

# Calculate PDR for the output bam file
pdr = args.output_parent + "/" + args.sample_id + '/PDR/'
cmd = manager.config.tools.rscript + ' ' + pipe_folder + manager.config.parameters.rscript.pdr_script + ' PDR ' + pdr + " " + sorted_bam + ' ' + args.output_parent + "/" + args.sample_id + '/rnbSet.zip ' + str(manager.config.parameters.rscript.cores)
manager.run(cmd,lock_name=args.sample_id+'locker')

#' Create MHL region of interest
roi = args.output_parent + "/" + args.sample_id
roi += "/roi.bed"
hapinfo_output = args.output_parent + "/" + args.sample_id + '/hapinfo.txt'
os.environ['PATH'] = manager.config.parameters.bismark.samtools + ':' + os.environ['PATH']
cmd = " ".join([manager.config.tools.perl, pipe_folder + manager.config.parameters.mhl.to_hapinfo, roi, sorted_bam, "bismark", roi, '>', hapinfo_output])
manager.run(cmd,hapinfo_output)

#' Calculate MHL from haplotype information
mhl_output = args.output_parent + "/" + args.sample_id + '/mhl.txt'
cmd = 'echo ' + hapinfo_output + " > " + args.output_parent + "/" + args.sample_id + "/info_list.txt; " + " ".join([manager.config.tools.perl, pipe_folder + manager.config.parameters.mhl.hapinfo_to_mhl, args.output_parent + "/" + args.sample_id + "/info_list.txt", '>', mhl_output])
manager.run(cmd,mhl_output)

# Run methclone software
methclone_output = args.output_parent + "/" + args.sample_id + '/methclone_tmp.txt.gz'
cmd = " ".join([manager.config.tools.methclone, sorted_bam, sorted_bam, methclone_output, "methclone", str(manager.config.parameters.methclone.meth_diff), str(manager.config.parameters.methclone.distance_cutoff), str(manager.config.parameters.methclone.coverage)])
manager.run(cmd,methclone_output)

#' Calculate epipolymorphism from methclone's output and remove this
epipoly_output = args.output_parent + "/" + args.sample_id + 'epipoly.csv'
cmd = " ".join([manager.config.tools.rscript, pipe_folder + manager.config.parameters.rscript.conversion_epipoly, methclone_output, 'epipoly', args.output_parent + "/" + args.sample_id])
manager.run(cmd,epipoly_output)

#' Calculate methylation entropy from methclone's output and remove this
entropy_output = args.output_parent + "/" + args.sample_id + 'entropy.csv'
cmd = " ".join([manager.config.tools.rscript, pipe_folder + manager.config.parameters.rscript.conversion_entropy, methclone_output, 'entropy', args.output_parent + "/" + args.sample_id + "/"])
manager.run(cmd,entropy_output)

#" Stop pipeline
manager.stop_pipeline()
