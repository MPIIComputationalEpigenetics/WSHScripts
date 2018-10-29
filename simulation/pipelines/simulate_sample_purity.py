#!/usr/bin/python

#####################################################################################################
#' simulate_sample_purity.py
#' This script runs the pipeline to create a sample with artificially created cellular contamination.
#' The final sample purity is drawn as a random number in the vector (0.5,0.6,0.7,0.8,0.9,1.0).
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

conversion_rate = random.randint(0,100)
if conversion_rate > 50:
	# Unmethylated region
	conversion_rate = 100
else:
	# Methylated region
	conversion_rate = 0

sample_purity = random.randint(5,10)
is_negative = random.randint(0,100)>50
out_file = args.output_parent + "/" + args.sample_id + '/purity.txt'
out_file = open(out_file,'wb')
out_string = str(sample_purity*10)+"%"
if is_negative:
	out_string = "Negative control\n" + out_string

out_file.write(out_string)
out_file.close()

genome_location = args.output_parent + "/" + args.sample_id + "/simulated.fa"
out_file = open(genome_location,'wb')
out_file.write(">" + args.chr + "\n" + ''.join(result_list))
out_file.close()

num_reads = round(manager.config.parameters.sherman.number/10)
op_start = random.randint(2,len(result_list))
op_end = random.randint(op_start,len(result_list)-1)
first_region = result_list[1:op_start]
dmr_region = result_list[op_start:op_end]
second_region = result_list[op_end:count]
first_region = ">" + args.chr + "\n" + ''.join(first_region)
second_region = ">" + args.chr + "\n" + ''.join(second_region)
dmr_region = ">" + args.chr + "\n" + ''.join(dmr_region)
contamination_location = args.output_parent + "/" + args.sample_id + '/contamination/' 
if not os.path.exists(contamination_location):
	os.makedirs(contamination_location)

if not os.path.exists(contamination_location+'first/'):
	os.makedirs(contamination_location+'first/')

out_file = open(contamination_location+'first/first_region.fa','wb')
out_file.write(first_region)
out_file.close()
if not os.path.exists(contamination_location+'second/'):
	os.makedirs(contamination_location+'second/')

out_file = open(contamination_location+'/second/second_region.fa','wb')
out_file.write(second_region)
out_file.close()
if not os.path.exists(contamination_location+'dmr/'):
	os.makedirs(contamination_location+'dmr/')

out_file = open(contamination_location+'/dmr/dmr_region.fa','wb')
out_file.write(dmr_region)
out_file.close()
out_file = open(contamination_location+'/dmr/dmr_location.txt','wb')
out_file.write(str(op_start*50) + "\n" + str(op_end*50))
out_file.close()
cell_meth = random.randint(0,100)

frac_first = float(op_start)/float(len(result_list))
frac_second = float(len(result_list)-op_end)/float(len(result_list))
frac_dmr = float(op_end-op_start)/float(len(result_list))

for ct in range(1,11):

	if sample_purity < ct:
		if is_negative:
			contamination_rate = conversion_rate
		else:
			if conversion_rate == 100:
				contamination_rate = 0
			else:
				contamination_rate=100
		
			# Run Sherman to simulate reads
		simulated = contamination_location + "/simulated/"
		manager.run("mkdir " + simulated,simulated)
		cmd = manager.config.tools.sherman + " --length " + str(manager.config.parameters.sherman.length) + " --number_of_seqs " + str(int(round(num_reads*frac_first))) + " --CG_conversion " + str(conversion_rate) + " --CH_conversion " + str(manager.config.parameters.sherman.CH) + " --genome_folder " + contamination_location + "/first/ --outfolder " + simulated + " --outname first " + " --error_rate " + str(manager.config.parameters.sherman.error_rate) + " --quality " + str(manager.config.parameters.sherman.quality) + " --rrbs"
		manager.run(cmd,lock_name=args.sample_id+'locker')

		# Run Sherman to simulate reads
		cmd = manager.config.tools.sherman + " --length " + str(manager.config.parameters.sherman.length) + " --number_of_seqs " + str(int(round(num_reads*frac_second))) + " --CG_conversion " + str(conversion_rate) + " --CH_conversion " + str(manager.config.parameters.sherman.CH) + " --genome_folder " + contamination_location + "/second/ --outfolder " + simulated + " --outname second " + " --error_rate " + str(manager.config.parameters.sherman.error_rate) + " --quality " + str(manager.config.parameters.sherman.quality) + " --rrbs"
		manager.run(cmd,lock_name=args.sample_id+'locker')

		# Run Sherman to simulate reads
		cmd = manager.config.tools.sherman + " --length " + str(manager.config.parameters.sherman.length) + " --number_of_seqs " + str(int(round(num_reads*frac_dmr))) + " --CG_conversion " + str(contamination_rate) + " --CH_conversion " + str(manager.config.parameters.sherman.CH) + " --genome_folder " + contamination_location + "/dmr/ --outfolder " + simulated + " --outname dmr " + " --error_rate " + str(manager.config.parameters.sherman.error_rate) + " --quality " + str(manager.config.parameters.sherman.quality) + " --rrbs"
		manager.run(cmd,lock_name=args.sample_id+'locker')

		# Merge the fastq files
		cmd = "cat " + simulated + "*.fastq > " + args.output_parent + "/" + args.sample_id + "/merged_" + str(ct) + ".fq"
		manager.run(cmd,lock_name=args.sample_id+'locker')
		
	else:  
	# Run Sherman to simulate reads
		simulated = args.output_parent + "/" + args.sample_id + "/"
		cmd = manager.config.tools.sherman + " --length " + str(manager.config.parameters.sherman.length) + " --number_of_seqs " + str(int(num_reads)) + " --CG_conversion " + str(conversion_rate) + " --CH_conversion " + str(manager.config.parameters.sherman.CH) + " --genome_folder " + args.output_parent + "/" + args.sample_id + " --outfolder " + simulated + " --outname " + str(ct) + " --error_rate " + str(manager.config.parameters.sherman.error_rate) + " --quality " + str(manager.config.parameters.sherman.quality) + " --rrbs"
		manager.run(cmd,lock_name=args.sample_id+'locker')

		cmd = "mv " + simulated + str(ct) + ".fastq " + simulated + str(ct) + ".fq "
		manager.run(cmd,lock_name=args.sample_id+"locko")

# Merge all cell types to a single file
#cmd = "cat " + args.output_parent + "/" + args.sample_id + "/*.fastq > "  + args.output_parent + "/" + args.sample_id + "/all_merged.fq; rm " + args.output_parent + "/" + args.sample_id + "/simulated.fa"
#manager.run(cmd,lock_name=args.sample_id+'lockard')

# Map the simulated fastq file to the reference with bimark
mapped = args.output_parent + "/" + args.sample_id + "/mapped/"
cmd = "mkdir " + mapped 
manager.run(cmd,lock_name=args.sample_id+'lockard')
cmd = manager.config.tools.bismark + " -o " + mapped + " " + manager.config.resources.bisulfite_folder + " --" + manager.config.parameters.bismark.bowtie + " --path_to_bowtie " + manager.config.parameters.bismark.bowtie_path + " --samtools_path " + manager.config.parameters.bismark.samtools + " " + args.output_parent + "/" + args.sample_id + "/*.fq"
manager.run(cmd,lock_name=args.sample_id+'lockard2')

cmd = manager.config.tools.samtools + " merge " + mapped + "all_merged.sam " + mapped + "*.sam"
manager.run(cmd,lock_name=args.sample_id+"lockard3")

# Create the cov files for investigation of methylation states
covs = args.output_parent + "/" + args.sample_id + '/covs/'
cmd = 'mkdir ' + covs + '; ' + manager.config.tools.bismark_extractor + ' -o ' + covs + ' --samtools_path ' + manager.config.parameters.bismark.samtools + ' --single-end ' + ' --bedGraph ' + mapped + '*.sam'
manager.run(cmd,covs)

# Create and store rnbSet
os.environ['LD_LIBRARY_PATH'] = manager.config.parameters.misc.ld
cmd = 'rm ' + covs + '*.txt; rm ' + covs + '*.png; rm ' + covs + '*.bedGraph; '+ 'echo "sample_id,filename\nTest,all_merged.bismark.cov" > ' + args.output_parent + "/" + args.sample_id + '/sample_annotation.csv; ' +  manager.config.tools.rscript + ' ' + pipe_folder + manager.config.parameters.rscript.create_rnb + ' ' + args.output_parent + "/" + args.sample_id + "; rm " + args.output_parent + "/" + args.sample_id + "/sample_annotation.csv"
manager.run(cmd,lock_name=args.sample_id+'locker')

# Write sam as bam file
bam = args.output_parent + "/" + args.sample_id + "/bam/"
cmd = "mkdir " + bam + '; ' +manager.config.tools.samtools + " view -b " + mapped + "all_merged.sam > " + bam + "simulated.bam"#; rm -r " + mapped
manager.run(cmd,bam)

# Sort bam file
sorted_bam = bam + "sorted.bam"
cmd = manager.config.tools.samtools + " sort " + bam + "* > " + sorted_bam + "; rm " + bam + "simulated.bam"
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
cmd = " ".join([manager.config.tools.perl,  pipe_folder + manager.config.parameters.mhl.to_hapinfo, roi, sorted_bam, "bismark", roi, '>', hapinfo_output])
manager.run(cmd,hapinfo_output)

#' Calculate MHL from haplotype information
mhl_output = args.output_parent + "/" + args.sample_id + '/mhl.txt'
cmd = 'echo ' + hapinfo_output + " > " + args.output_parent + "/" + args.sample_id + "/info_list.txt; " + " ".join([manager.config.tools.perl,  pipe_folder + manager.config.parameters.mhl.hapinfo_to_mhl, args.output_parent + "/" + args.sample_id + "/info_list.txt", '>', mhl_output])
manager.run(cmd,mhl_output)

# Run methclone software
methclone_output = args.output_parent + "/" + args.sample_id + '/methclone_tmp.txt.gz'
cmd = " ".join([manager.config.tools.methclone, sorted_bam, sorted_bam, methclone_output, "methclone", str(manager.config.parameters.methclone.meth_diff), str(manager.config.parameters.methclone.distance_cutoff), str(manager.config.parameters.methclone.coverage)])
manager.run(cmd,methclone_output)

#' Calculate epipolymorphism from methclone's output and remove this
epipoly_output = args.output_parent + "/" + args.sample_id + 'epipoly.csv'
cmd = " ".join([manager.config.tools.rscript,  pipe_folder + manager.config.parameters.rscript.conversion_epipoly, methclone_output, 'epipoly', args.output_parent + "/" + args.sample_id]) + ' ; rm ' + methclone_output
manager.run(cmd,epipoly_output)

#' Calculate methylation entropy from methclone's output and remove this
entropy_output = args.output_parent + "/" + args.sample_id + 'entropy.csv'
cmd = " ".join([manager.config.tools.rscript,  pipe_folder + manager.config.parameters.rscript.conversion_entropy, methclone_output, 'entropy', args.output_parent + "/" + args.sample_id]) + ' ; rm ' + methclone_output
manager.run(cmd,entropy_output)

#' cleanup
cmd = 'rm ' + methclone_output + '; rm -Rf ' + bam + '; rm -Rf ' + simulated + '; rm -Rf ' + args.output_parent + "/" + args.sample_id + '/first' + args.output_parent + "/" + args.sample_id + '/second'
manager.run(cmd,lock_name="locki")

#" Stop pipeline
manager.stop_pipeline()
