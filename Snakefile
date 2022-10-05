

## Configuration file
import os
if len(config) == 0:
  if os.path.isfile("./config.yaml"):
    configfile: "./config.yaml"
  else:
    sys.exit("Make sure there is a config.yaml file in " + os.getcwd() + " or specify one with the --configfile commandline parameter.")

## Make sure that all expected variables from the config file are in the config dictionary
configvars = ['annotation', 'organism', 'build', 'release', 'txome', 'genome', 'gtf', 'salmonindex', 'salmonk', 'STARindex', 'HISAT2index', 'readlength', 'fldMean', 'fldSD', 'metatxt', 'design', 'contrast', 'genesets', 'ncores', 'FASTQ', 'fqext1', 'fqext2', 'fqsuffix', 'output', 'useCondaR', 'Rbin', 'run_trimming', 'run_STAR', 'run_HISAT2', 'run_DRIMSeq', 'run_camera']
for k in configvars:
	if k not in config:
		config[k] = None

## If any of the file paths is missing, replace it with ""
def sanitizefile(str):
	if str is None:
		str = ''
	return str
	
config['txome'] = sanitizefile(config['txome'])
config['gtf'] = sanitizefile(config['gtf'])
config['genome'] = sanitizefile(config['genome'])
config['STARindex'] = sanitizefile(config['STARindex'])
config['HISAT2index'] = sanitizefile(config['HISAT2index'])
config['salmonindex'] = sanitizefile(config['salmonindex'])
config['metatxt'] = sanitizefile(config['metatxt'])
config['metacsv'] = os.path.splitext(os.path.basename(config['metatxt']))[0]+'.csv'

# ## Read metadata
# if not os.path.isfile(config["metatxt"]):
#   sys.exit("Metadata file " + config["metatxt"] + " does not exist.")

# import pandas as pd
# samples = pd.read_csv(config["metatxt"], sep='\t')
from pathlib import Path
from pandas_path import path
import pandas as pd

samples = pd.DataFrame((str(s) for s in Path('../data/FASTQ').rglob('*.fastq.gz')), columns = ["filepath"])
samples['lane'] = samples['filepath'].path.name.str.replace("_R\d.*","", regex = True)
samples['name'] = samples['lane'].str.replace("_S1_L00.*", "", regex = True)
# samples['name'] = samples['lane'].str.replace("_[0-9].*", "", regex = True)

samples['type'] = config['read_type']

# if not set(['names','type']).issubset(samples.columns):
#   sys.exit("Make sure 'names' and 'type' are columns in " + config["metatxt"])


## Sanitize provided input and output directories
import re
def getpath(str):
	if str in ['', '.', './']:
		return ''
	if str.startswith('./'):
		regex = re.compile('^\./?')
		str = regex.sub('', str)
	if not str.endswith('/'):
		str += '/'
	return str

proj_dir = os.path.abspath(os.path.normpath(getpath(config["proj_dir"])))
outputdir = os.path.abspath(getpath(config["output"])) + "/"
FASTQdir = getpath(config["FASTQ"])
datadir = os.path.abspath(os.path.normpath(getpath(config["proj_dir"]))) + "/data/"

clones_filename = "bulk_clones_final.png"

## Define the conda environment for all rules using R
if config["useCondaR"] == True:
	Renv = "envs/environment_R.yaml"
else:
	Renv = "envs/environment.yaml"

## Define the R binary
Rbin = config["Rbin"]

# The config.yaml files determines which steps should be performed
	
def dbtss_output(wildcards):
	input = []
	input.extend(expand(outputdir + "dbtss_coverage/{sample}_dbtss_coverage_over_10.txt", sample = samples.name[samples.type == 'PE'].values.tolist()))
	return input
	
def bigwigoutput(wildcards):
  input = []
  input.extend(expand(outputdir + "STAR/{sample}_Aligned.sortedByCoord.out.bw", sample = samples.name[samples.type == 'PE'].values.tolist()))
  return input

def dexseqoutput(wildcards):
  input = []
  input.extend(expand(outputdir + "dexseq/{sample}.txt", sample = samples.name[samples.type == 'PE'].values.tolist()))
  return input
  
def jbrowse_output(wildcards):
  input = []
  input.extend(expand("/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/trackList.json"))
  # input.extend(expand("/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/output"))
  input.extend(expand("/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/samples/{sample}.bw", sample = samples.name[samples.type == 'PE'].values.tolist()))
  input.extend(expand("/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/samples/{sample}.bam", sample = samples.name[samples.type == 'PE'].values.tolist()))
  input.extend(expand("/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/samples/{sample}.bam.bai", sample = samples.name[samples.type == 'PE'].values.tolist()))
  
  input.append("/var/www/html/jbrowse/" + os.path.basename(proj_dir) + "/trackList.json")
  return input

# def kb_output(wildcards):
# 	input = []
# 	input.extend(expand(outputdir + "kb/unfiltered/adata.h5ad", sample = samples.name[samples.type == 'PE'].values.tolist()))
# 	return input

def cellranger_output(wildcards):
  input = []
  input.extend(expand(outputdir + "cellranger/{sample}/outs/possorted_genome_bam.bam", sample = samples.name[samples.type == 'PE'].values.tolist()))
  return input

def numbat_output(wildcards):
  input = []
  input.extend(expand(outputdir + "numbat/{sample}/done.txt", sample = samples.name[samples.type == 'PE'].values.tolist()))
  # input.extend(expand(outputdir + "numbat/{sample}/{sample}_report.html", sample = samples.name[samples.type == 'PE'].values.tolist()))
  return input
 
def infercnv_output(wildcards):
  input = []
  input.extend(expand(outputdir + "infercnv/{sample}/infercnv.png", sample = samples.name[samples.type == 'PE'].values.tolist()))
  return input
 
  
	
## ------------------------------------------------------------------------------------ ##
## Target definitions
## ------------------------------------------------------------------------------------ ##
## Run all analyses
rule all:
	input:
		# outputdir + "MultiQC/multiqc_report.html",
		cellranger_output,
		numbat_output,
		infercnv_output
		# bigwigoutput,
		# outputdir + "seurat/unfiltered_seu.rds",
		# outputdir + "kallisto/adata.h5ad",
		# stringtie_output,
		# outputdir + "seurat/legacy_unfiltered_seu.rds",
		# dexseqoutput,
		# dbtss_output,
		# jbrowse_output,
		# loom_file = outputdir + "velocyto/" + os.path.basename(proj_dir) + ".loom",
		# velocyto_seu = outputdir + "velocyto/" + "unfiltered_seu.rds",
		# loom = outputdir + "scenic/unfiltered.loom",
		# final_loom = outputdir + "scenic/unfiltered-final.loom"

rule setup:
	input:
		outputdir + "Rout/pkginstall_state.txt",
		outputdir + "Rout/softwareversions.done"

## Install R packages
rule pkginstall:
	input:
		script = "scripts/install_pkgs.R"
	output:
	  outputdir + "Rout/pkginstall_state.txt"
	params:
		flag = config["annotation"],
		ncores = config["ncores"],
		organism = config["organism"]
	priority:
		50
	conda:
		Renv
	log:
		outputdir + "Rout/install_pkgs.Rout"
	benchmark:
	  outputdir + "benchmarks/install_pkgs.txt"
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args outtxt='{output}' ncores='{params.ncores}' annotation='{params.flag}' organism='{params.organism}'" {input.script} {log}'''

## FastQC on original (untrimmed) files
rule runfastqc:
	input:
		expand(outputdir + "FastQC/{sample}_" + str(config["fqext1"]) + "_fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "FastQC/{sample}_" + str(config["fqext2"]) + "_fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "FastQC/{sample}_fastqc.zip", sample = samples.name[samples.type == 'SE'].values.tolist())

## Trimming and FastQC on trimmed files
rule runtrimming:
	input:
		expand(outputdir + "FastQC/{sample}_" + str(config["fqext1"]) + "_val_1_fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "FastQC/{sample}_" + str(config["fqext2"]) + "_val_2_fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "FastQC/{sample}_trimmed_fastqc.zip", sample = samples.name[samples.type == 'SE'].values.tolist())

## Salmon quantification
rule runsalmonquant:
	input:
		expand(outputdir + "salmon/{sample}/quant.sf", sample = samples.name.values.tolist())

## STAR alignment
rule runstar:
	input:
		expand(outputdir + "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai", sample = samples.name.values.tolist()),
		expand(outputdir + "STAR/{sample}_Aligned.sortedByCoord.out.bw", sample = samples.name.values.tolist())
		
## HISAT2 alignment
rule runhisat2:
	input:
		expand(outputdir + "HISAT2/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai", sample = samples.name.values.tolist()),
		expand(outputdir + "HISAT2bigwig/{sample}_Aligned.sortedByCoord.out.all.bw", sample = samples.name.values.tolist())
		
## DBTSS coverage
rule rundbtss:
	input:
		expand(outputdir + "dbtss_coverage/{sample}/{sample}_dbtss_coverage.txt", sample = samples.name.values.tolist())

## List all the packages that were used by the R analyses
rule listpackages:
	log:
		outputdir + "Rout/list_packages.Rout"
	params:
		Routdir = outputdir + "Rout",
		outtxt = outputdir + "R_package_versions.txt",
		script = "scripts/list_packages.R"
	conda:
		Renv
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args Routdir='{params.Routdir}' outtxt='{params.outtxt}'" {params.script} {log}'''

## Print the versions of all software packages
rule softwareversions:
	output:
		touch(outputdir + "Rout/softwareversions.done")
	conda:
		"envs/environment.yaml"
	shell:
		"echo -n 'ARMOR version ' && cat version; "
		"salmon --version; kallisto version; trim_galore --version; "
		"echo -n 'cutadapt ' && cutadapt --version; "
		"fastqc --version; STAR --version; hisat2 --version; samtools --version; multiqc --version; "
		"bedtools --version"

## ------------------------------------------------------------------------------------ ##
## Reference preparation
## ------------------------------------------------------------------------------------ ##

## Generate Salmon index from merged cDNA and ncRNA files
rule salmonindex:
	input:
		txome = config["txome"]
	output:
		config["salmonindex"] + "/hash.bin"
	log:
		outputdir + "logs/salmon_index.log"
	benchmark:
		outputdir + "benchmarks/salmon_index.txt"
	params:
		salmonk = config["salmonk"],
		salmonoutdir = config["salmonindex"],
		anno = config["annotation"]
	conda:
		"envs/environment.yaml"
	shell:
	  """
	  if [ {params.anno} == "Gencode" ]; then
      echo 'Salmon version:\n' > {log}; salmon --version >> {log};
  	  salmon index -t {input.txome} -k {params.salmonk} -i {params.salmonoutdir} --gencode --type quasi

    else
  	  echo 'Salmon version:\n' > {log}; salmon --version >> {log};
      salmon index -t {input.txome} -k {params.salmonk} -i {params.salmonoutdir} --type quasi
    fi
    """
    
# kb ref -i transcriptome.idx -g transcripts_to_genes.txt -f1 cdna.fa dna.primary_assembly.fa.gz gtf.gz

## Generate linkedtxome mapping
rule linkedtxome:
	input:
		txome = config["txome"],
		gtf = config["gtf"],
		salmonidx = config["salmonindex"] + "/hash.bin",
		script = "scripts/generate_linkedtxome.R",
		install = outputdir + "Rout/pkginstall_state.txt"
	log:
		outputdir + "Rout/generate_linkedtxome.Rout"
	benchmark:
		outputdir + "benchmarks/generate_linkedtxome.txt"
	output:
		config["salmonindex"] + ".json"
	params:
		flag = config["annotation"],
		organism = config["organism"],
		release = str(config["release"]),
		build = config["build"]
	conda:
		Renv
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args transcriptfasta='{input.txome}' salmonidx='{input.salmonidx}' gtf='{input.gtf}' annotation='{params.flag}' organism='{params.organism}' release='{params.release}' build='{params.build}' output='{output}'" {input.script} {log}'''

## Generate STAR index
# rule starindex:
# 	input:
# 		genome = config["genome"],
# 		gtf = config["gtf"]
# 	output:
# 		config["STARindex"] + "/SA",
# 		config["STARindex"] + "/chrNameLength.txt"
# 	log:
# 		outputdir + "logs/STAR_index.log"
# 	benchmark:
# 		outputdir + "benchmarks/STAR_index.txt"
# 	params:
# 		STARindex = config["STARindex"],
# 		readlength = config["readlength"]
# 	conda:
# 		"envs/environment.yaml"
# 	threads:
# 		config["ncores"]
# 	shell:
# 		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
# 		"STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.STARindex} "
# 		"--genomeFastaFiles {input.genome} --sjdbGTFfile {input.gtf} --sjdbOverhang {params.readlength}"

## ------------------------------------------------------------------------------------ ##
## Quality control
## ------------------------------------------------------------------------------------ ##

		
## Rseqc gene body coverage plot
rule genebodycoverage:
	input:
		bigwig = outputdir + "HISAT2bigwig/{sample}_Aligned.sortedByCoord.out.all.bw"
	output:
		coverage_txt = outputdir + "rseqc/{sample}.geneBodyCoverage.txt"
	params:
		sample = outputdir + "rseqc/{sample}",
		bed = config["bed"]
	log:
		outputdir + "logs/genebodycoverage_{sample}.log"
	benchmark:
		outputdir + "benchmarks/genebodycoverage_{sample}.txt"
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"echo 'geneBody_coverage.py version:\n' > {log}; geneBody_coverage.py --version >> {log}; "
		"geneBody_coverage2.py -r {params.bed} -i {input.bigwig}  -o {params.sample}"

# The config.yaml files determines which steps should be performed
def multiqc_input(wildcards):
	input = []
	input.extend(expand(outputdir + "FastQC/{sample}_fastqc.zip", sample = samples.name[samples.type == 'SE'].values.tolist()))
	input.extend(expand(outputdir + "FastQC/{sample}_" + str(config["fqext1"]) + "_fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()))
	input.extend(expand(outputdir + "FastQC/{sample}_" + str(config["fqext2"]) + "_fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()))
	if config["run_genebodycoverage"]:
	  input.extend(expand(outputdir + "rseqc/{sample}." + "geneBodyCoverage.txt", sample = samples.name[samples.type == 'PE'].values.tolist()))
	if config["run_SALMON"]:
		input.extend(expand(outputdir + "salmon/{sample}/quant.sf", sample = samples.name.values.tolist()))
	if config["run_trimming"]:
		# input.extend(expand(outputdir + "FASTQtrimmed/{sample}_trimmed.fq.gz", sample = samples.name[samples.type == 'SE'].values.tolist()))
		input.extend(expand(outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext1"]) + "_val_1.fq.gz", sample = samples.name[samples.type == 'PE'].values.tolist()))
		input.extend(expand(outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext2"]) + "_val_2.fq.gz", sample = samples.name[samples.type == 'PE'].values.tolist()))
		input.extend(expand(outputdir + "FastQC/{sample}_trimmed_fastqc.zip", sample = samples.name[samples.type == 'SE'].values.tolist()))
		input.extend(expand(outputdir + "FastQC/{sample}_" + str(config["fqext1"]) + "_val_1_fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()))
		input.extend(expand(outputdir + "FastQC/{sample}_" + str(config["fqext2"]) + "_val_2_fastqc.zip", sample = samples.name[samples.type == 'PE'].values.tolist()))
	if config["run_STAR"]:
		input.extend(expand(outputdir + "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai", sample = samples.name.values.tolist()))
	if config["run_HISAT2"]:
		input.extend(expand(outputdir + "HISAT2/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai", sample = samples.name.values.tolist()))
	return input

## Determine the input directories for MultiQC depending on the config file
def multiqc_params(wildcards):
	param = [outputdir + "FastQC"]
	if config["run_SALMON"]:
		param.append(outputdir + "salmon")
	if config["run_trimming"]:
		param.append(outputdir + "FASTQtrimmed")
	if config["run_STAR"]:
		param.append(outputdir + "STAR")
	if config["run_HISAT2"]:
		param.append(outputdir + "HISAT2")
	return param

## MultiQC
rule multiqc:
	input:
		multiqc_input
	output:
		outputdir + "MultiQC/multiqc_report.html"
	params:
		inputdirs = multiqc_params,
		MultiQCdir = outputdir + "MultiQC"
	log:
		outputdir + "logs/multiqc.log"
	benchmark:
		outputdir + "benchmarks/multiqc.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'MultiQC version:\n' > {log}; multiqc --version >> {log}; "
		"multiqc {params.inputdirs} -f -o {params.MultiQCdir}"


## ------------------------------------------------------------------------------------ ##
## Adapter trimming
## ------------------------------------------------------------------------------------ ##
# TrimGalore!
# rule trimgaloreSE:
# 	input:
# 		fastq = FASTQdir + "{sample}." + str(config["fqsuffix"]) + ".gz"
# 	output:
# 		outputdir + "FASTQtrimmed/{sample}_trimmed.fq.gz"
# 	params:
# 		FASTQtrimmeddir = outputdir + "FASTQtrimmed"
# 	log:
# 		outputdir + "logs/trimgalore_{sample}.log"
# 	benchmark:
# 		outputdir + "benchmarks/trimgalore_{sample}.txt"
# 	conda:
# 		"envs/environment.yaml"
# 	shell:
# 		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
# 		"trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt {input.fastq}"

rule trimgalorePE:
	input:
		fastq1 = FASTQdir + "{sample}_" + str(config["fqext1"]) + "." + str(config["fqsuffix"]) + ".gz",
		fastq2 = FASTQdir + "{sample}_" + str(config["fqext2"]) + "." + str(config["fqsuffix"]) + ".gz"
	output:
		intermediate_fastq1 = outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext1"]) + "_val_1.fq.gz",
		intermediate_fastq2 = outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext2"]) + "_val_2.fq.gz"
	params:
		FASTQtrimmeddir = outputdir + "FASTQtrimmed"
	log:
		outputdir + "logs/trimgalore_{sample}.log"
	benchmark:
		outputdir + "benchmarks/trimgalore_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
		"trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt "
		"--paired {input.fastq1} {input.fastq2}; "

## ------------------------------------------------------------------------------------ ##
## HISAT2 mapping
## ------------------------------------------------------------------------------------ ##
## Genome mapping with HISAT2
rule HISAT2PE:
	input:
		fastq1 = outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext1"]) + "_val_1.fq.gz" if config["run_trimming"] else FASTQdir + "{sample}_" + str(config["fqext1"]) + "." + str(config["fqsuffix"]) + ".gz",
		fastq2 = outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext2"]) + "_val_2.fq.gz" if config["run_trimming"] else FASTQdir + "{sample}_" + str(config["fqext2"]) + "." + str(config["fqsuffix"]) + ".gz"
	output:
		bam = temp(outputdir + "HISAT2/{sample}/{sample}_Aligned.out.bam")
	threads:
		config["ncores"]
	log:
		version = outputdir + "logs/HISAT2_{sample}.log",
		stats = outputdir + "HISAT2" + "/HISAT2_{sample}_stats.txt"
	benchmark:
		outputdir + "benchmarks/HISAT2_{sample}.txt"
	params:
		HISAT2index = config["HISAT2index"],
		HISAT2dir = outputdir + "HISAT2"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'hisat2 --version:\n' > {log.version}; hisat2 --version >> {log.version}; "
		"hisat2 --new-summary --pen-noncansplice 20 --threads {threads} --mp 1,0 --sp 3,1 -x {params.HISAT2index} -1 {input.fastq1} -2 {input.fastq2} 2> {log.stats} | samtools view -Sbo {output.bam}"

# convert and sort sam files
rule bamsort:
	input:
		bam = outputdir + "HISAT2/{sample}/{sample}_Aligned.out.bam"
	output:
		sorted_bam = outputdir + "HISAT2/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	log:
		outputdir + "logs/samtools_sort_{sample}.log"
	benchmark:
		outputdir + "benchmarks/samtools_sort_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'samtools version:\n' > {log}; samtools --version >> {log}; "
		"samtools sort -O bam -o {output.sorted_bam} {input.bam}"

## Index bam files
rule bamindexhisat2:
	input:
		bam = outputdir + "HISAT2/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		outputdir + "HISAT2/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
	log:
		outputdir + "logs/samtools_index_{sample}.log"
	benchmark:
		outputdir + "benchmarks/samtools_index_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'samtools version:\n' > {log}; samtools --version >> {log}; "
		"samtools index {input.bam}"

## Convert gdna BAM files to bigWig
rule bigwighisat2:
	input:
		bam = outputdir + "HISAT2/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		outputdir + "HISAT2bigwig/{sample}_Aligned.sortedByCoord.out.all.bw"
	params:
	  prefix = outputdir + "HISAT2bigwig/{sample}_Aligned.sortedByCoord.out"
		# HISAT2bigwigdir = outputdir + "HISAT2bigwig"
	log:
		outputdir + "logs/bigwig_{sample}.log"
	benchmark:
		outputdir + "benchmarks/bigwig_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
	  "megadepth {input.bam} --threads {threads} --bigwig --prefix {params.prefix}"

## ------------------------------------------------------------------------------------ ##
## Stringtie
## ------------------------------------------------------------------------------------ ##
# Transcript assembly using StringTie
rule stringtie:
	input:
		bam = outputdir + "HISAT2/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		gtf = outputdir + "stringtie/{sample}/{sample}.gtf"
	log:
		outputdir + "logs/stringtie_{sample}.log"
	benchmark:
		outputdir + "benchmarks/stringtie_{sample}.txt"
	threads:
		config["ncores"]
	params:
		stringtiegtf = config["gtf"],
		stringtiedir = outputdir + "stringtie"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'stringtie version:\n' > {log}; stringtie --version >> {log}; "
		"stringtie {input.bam} -G {params.stringtiegtf} -x MT -eB -o {output.gtf}"
>>>>>>> 270a8379b4439f52fab5e013fc019bda77a25874

## ------------------------------------------------------------------------------------ ##
## DEXSeq
## ------------------------------------------------------------------------------------ ##

# Transcript assembly using dexseq
rule dexseq:
	input:
		bam = outputdir + "HISAT2/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		txt = outputdir + "dexseq/{sample}.txt"
	log:
		outputdir + "logs/dexseq_{sample}.log"
	benchmark:
		outputdir + "benchmarks/dexseq_{sample}.txt"
	threads:
		config["ncores"]
	params:
		dexseqgtf = config["exon_collapsed_gff"],
		dexseqdir = outputdir + "dexseq"
	conda:
		"envs/environment.yaml"
	shell:
		# "echo 'dexseq version:\n' > {log}; dexseq --version >> {log}; "
		"python scripts/dexseq_count.py -r pos -p yes -s no -f bam {params.dexseqgtf} {input.bam} {output.txt}"

## ------------------------------------------------------------------------------------ ##
## kb abundance estimation
## ------------------------------------------------------------------------------------ ##

rule kallisto_batch:
  input:
    metatxt = config["metatxt"]
  output:
    outputdir + "kallisto/smart-seq_batch.txt"
	log:
		outputdir + "logs/kallisto.log"
	benchmark:
		outputdir + "benchmarks/kallisto.txt"
	threads:
		config["ncores"]
	run:
		# make batch
		# samples = pd.read_csv(config["metatxt"], sep='\t')
		batch_table = samples[['names', 'type']]
		batch_table = batch_table.assign(R1="../data/FASTQ/"+batch_table['names']+"_R1.fastq.gz", R2="../data/FASTQ/"+batch_table['names']+"_R2.fastq.gz")
		batch_table = batch_table.drop(['type'], axis=1)
		batch_table.to_csv(output[0], sep = "\t", index = False, header = False)

# Estimate abundances with kb
rule kallistoPE:
	input:
	  kallistobatch = outputdir + "kallisto/smart-seq_batch.txt"
	output:
		outputdir + "kallisto/matrix.abundance.mtx"
	log:
		outputdir + "logs/kallisto.log"
	benchmark:
		outputdir + "benchmarks/kallisto.txt"
	threads:
		config["ncores"]
	params:
		kallistoindex = config["kallistoindex"],
		kallistodir = outputdir + "kallisto",
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'kallisto version:\n' > {log}; kallisto version >> {log}; "
		"kallisto pseudo --quant -i {params.kallistoindex} -o {params.kallistodir} -b {input.kallistobatch} -t {threads}"
		
# convert kallisto output to anndata
rule kallisto2adata:
	input:
	  matrix = outputdir + "kallisto/matrix.abundance.mtx",
	  cells = outputdir + "kallisto/matrix.cells",
	  transcripts = outputdir + "kallisto/transcripts.txt",
	  script = "scripts/make_adata.py"
	output:
		outputdir + "kallisto/adata.h5ad"
	log:
		outputdir + "logs/kallisto.log"
	benchmark:
		outputdir + "benchmarks/kallisto.txt"
	threads:
		config["ncores"]
	params:
		kallistodir = outputdir + "kallisto/",
	conda:
		"envs/environment.yaml"
	shell:
		"python {input.script} -m {input.matrix} -c {input.cells} -t {input.transcripts} --outdir {params.kallistodir}"

## ------------------------------------------------------------------------------------ ##
## Salmon abundance estimation
## ------------------------------------------------------------------------------------ ##
# Estimate abundances with Salmon
rule salmonSE:
	input:
		index = config["salmonindex"] + "/hash.bin",
		fastq = outputdir + "FASTQtrimmed/{sample}_trimmed.fq.gz" if config["run_trimming"] else FASTQdir + "{sample}." + str(config["fqsuffix"]) + ".gz"
	output:
		outputdir + "salmon/{sample}/quant.sf"
	log:
		outputdir + "logs/salmon_{sample}.log"
	benchmark:
		outputdir + "benchmarks/salmon_{sample}.txt"
	threads:
		config["ncores"]
	params:
		salmonindex = config["salmonindex"],
		fldMean = config["fldMean"],
		fldSD = config["fldSD"],
		salmondir = outputdir + "salmon"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'Salmon version:\n' > {log}; salmon --version >> {log}; "
		"salmon quant -i {params.salmonindex} -l A -r {input.fastq} "
		"-o {params.salmondir}/{wildcards.sample} --seqBias --gcBias "
		"--fldMean {params.fldMean} --fldSD {params.fldSD} -p {threads}"

rule salmonPE:
	input:
		index = config["salmonindex"] + "/hash.bin",
		fastq1 = outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext1"]) + "_val_1.fq.gz" if config["run_trimming"] else FASTQdir + "{sample}_" + str(config["fqext1"]) + "." + str(config["fqsuffix"]) + ".gz",
		fastq2 = outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext2"]) + "_val_2.fq.gz" if config["run_trimming"] else FASTQdir + "{sample}_" + str(config["fqext2"]) + "." + str(config["fqsuffix"]) + ".gz"
	output:
		outputdir + "salmon/{sample}/quant.sf"
	log:
		outputdir + "logs/salmon_{sample}.log"
	benchmark:
		outputdir + "benchmarks/salmon_{sample}.txt"
	threads:
		config["ncores"]
	params:
		salmonindex = config["salmonindex"],
		fldMean = config["fldMean"],
		fldSD = config["fldSD"],
		salmondir = outputdir + "salmon"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'Salmon version:\n' > {log}; salmon --version >> {log}; "
		"salmon quant -i {params.salmonindex} -l A -1 {input.fastq1} -2 {input.fastq2} "
		"-o {params.salmondir}/{wildcards.sample} --seqBias --gcBias "
		"--fldMean {params.fldMean} --fldSD {params.fldSD} -p {threads}"

## ------------------------------------------------------------------------------------ ##
## STAR mapping
## ------------------------------------------------------------------------------------ ##
## Genome mapping with STAR
# rule starSE:
# 	input:
# 		index = config["STARindex"] + "/SA",
# 		fastq = outputdir + "FASTQtrimmed/{sample}_trimmed.fq.gz" if config["run_trimming"] else FASTQdir + "{sample}." + str(config["fqsuffix"]) + ".gz"
# 	output:
# 		outputdir + "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
# 	threads:
# 		config["ncores"]
# 	log:
# 		outputdir + "logs/STAR_{sample}.log"
# 	benchmark:
# 		outputdir + "benchmarks/STAR_{sample}.txt"
# 	params:
# 		STARindex = config["STARindex"],
# 		STARdir = outputdir + "STAR"
# 	conda:
# 		"envs/environment.yaml"
# 	shell:
# 		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
# 		"STAR --genomeDir {params.STARindex} --readFilesIn {input.fastq} "
# 		"--runThreadN {threads} --outFileNamePrefix {params.STARdir}/{wildcards.sample}/{wildcards.sample}_ "
# 		"--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c"
# 
# rule unzip_star_input:
# 	input:
# 		fastq1 = outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext1"]) + "_val_1.fq.gz" if config["run_trimming"] else FASTQdir + "{sample}_" + str(config["fqext1"]) + "." + str(config["fqsuffix"]) + ".gz",
# 		fastq2 = outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext2"]) + "_val_2.fq.gz" if config["run_trimming"] else FASTQdir + "{sample}_" + str(config["fqext2"]) + "." + str(config["fqsuffix"]) + ".gz"
# 	output:
# 		temp(outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext1"]) + "_val_1.fq" if config["run_trimming"] else FASTQdir + "{sample}_" + str(config["fqext1"]) + "." + str(config["fqsuffix"])),
# 		temp(outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext2"]) + "_val_2.fq" if config["run_trimming"] else FASTQdir + "{sample}_" + str(config["fqext2"]) + "." + str(config["fqsuffix"]))
# 	params:
# 		processes = "-p"+str(config["ncores"])
# 	shell:
# 		"unpigz --keep {params.processes} {input.fastq1}; "
# 		"unpigz --keep {params.processes} {input.fastq2}"

rule starPE:
	input:
		index = config["STARindex"] + "/SA",
		fastq1 = outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext1"]) + "_val_1.fq" if config["run_trimming"] else FASTQdir + "{sample}_" + str(config["fqext1"]) + "." + str(config["fqsuffix"]),
		fastq2 = outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext2"]) + "_val_2.fq" if config["run_trimming"] else FASTQdir + "{sample}_" + str(config["fqext2"]) + "." + str(config["fqsuffix"])
	output:
		outputdir + "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	threads:
		config["ncores"]
	log:
		outputdir + "logs/STAR_{sample}.log"
	benchmark:
		outputdir + "benchmarks/STAR_{sample}.txt"
	params:
		STARindex = config["STARindex"],
		STARdir = outputdir + "STAR"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --genomeDir {params.STARindex} --readFilesIn {input.fastq1} {input.fastq2} "
		"--runThreadN {threads} --outFileNamePrefix {params.STARdir}/{wildcards.sample}/{wildcards.sample}_ "
		"--outSAMtype BAM SortedByCoordinate --outSAMattrIHstart 0"
		
rule prep_cellranger_input:
	input:
		fastq1 = outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext1"]) + "_val_1.fq.gz" if config["run_trimming"] else FASTQdir + "{sample}_" + str(config["fqext1"]) + "." + str(config["fqsuffix"]),
		fastq2 = outputdir + "FASTQtrimmed/{sample}_" + str(config["fqext2"]) + "_val_2.fq.gz" if config["run_trimming"] else FASTQdir + "{sample}_" + str(config["fqext2"]) + "." + str(config["fqsuffix"])
	output:
		fastq1 = outputdir + "FASTQtrimmed/{sample}/{sample}" + "_S1_L001_R1_001.fastq.gz",
		fastq2 = outputdir + "FASTQtrimmed/{sample}/{sample}" + "_S1_L001_R2_001.fastq.gz"
	threads:
		config["ncores"]
	conda:
		"envs/environment.yaml"
	shell:
		"cp {input.fastq1} {output.fastq1}; "
		"cp {input.fastq2} {output.fastq2}"
		
rule cellranger:
	input:
		fastqdir = FASTQdir + "{sample}"
	output:
		bam = outputdir + "cellranger/{sample}/outs/possorted_genome_bam.bam",
		barcodes = outputdir + "cellranger/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
		matrix_dir = directory(outputdir + "cellranger/{sample}/outs/filtered_feature_bc_matrix")
	threads:
		config["ncores"]
	log:
		outputdir + "logs/cellranger_{sample}.log"
	benchmark:
		outputdir + "benchmarks/cellranger_{sample}.txt"
	params:
		cellranger_transcriptome = config["cellranger_transcriptome"],
		cellranger_dir = outputdir + "cellranger"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'cellranger version:\n' > {log}; cellranger --version >> {log}; "
		"rm -rf {params.cellranger_dir}/{wildcards.sample}; "
		"cellranger count --id={wildcards.sample} --fastqs={input.fastqdir} "
		" --sample={wildcards.sample} --transcriptome={params.cellranger_transcriptome} --localmem 40; "
		"mv {wildcards.sample} {params.cellranger_dir}/{wildcards.sample}"

## Index bam files
rule bamindex:
	input:
		bam = outputdir + "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		outputdir + "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
	log:
		outputdir + "logs/samtools_index_{sample}.log"
	benchmark:
		outputdir + "benchmarks/samtools_index_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'samtools version:\n' > {log}; samtools --version >> {log}; "
		"samtools index {input.bam}"

## Convert BAM files to bigWig
rule bigwig:
	input:
		bam = outputdir + "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
		chrl = config["STARindex"] + "/chrNameLength.txt"
	output:
		outputdir + "STAR/{sample}_Aligned.sortedByCoord.out.bw"
	params:
		STARbigwigdir = outputdir + "STAR"
	log:
		outputdir + "logs/bigwig_{sample}.log"
	benchmark:
		outputdir + "benchmarks/bigwig_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'bedtools version:\n' > {log}; bedtools --version >> {log}; "
		"bedtools genomecov -split -ibam {input.bam} -bg | LC_COLLATE=C sort -k1,1 -k2,2n > "
		"{params.STARbigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph; "
		"bedGraphToBigWig {params.STARbigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph "
		"{input.chrl} {output}; rm -f {params.STARbigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph"

## ------------------------------------------------------------------------------------ ##
## Transcript quantification
## ------------------------------------------------------------------------------------ ##

## tximeta
rule tximeta:
	input:
	  outputdir + "Rout/pkginstall_state.txt",
		expand(outputdir + "salmon/{sample}/quant.sf", sample = samples.name.values.tolist()),
		metatxt = config["metatxt"],
		salmonidx = config["salmonindex"] + "/hash.bin",
		json = config["salmonindex"] + ".json",
		script = "scripts/run_tximeta.R"
	output:
		outputdir + "outputR/tximeta.rds"
	log:
		outputdir + "Rout/tximeta.Rout"
	benchmark:
		outputdir + "benchmarks/tximeta.txt"
	params:
		salmondir = outputdir + "salmon",
		flag = config["annotation"],
		organism = config["organism"]
	conda:
		Renv
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args salmondir='{params.salmondir}' json='{input.json}' metafile='{input.metatxt}' outrds='{output}' annotation='{params.flag}' organism='{params.organism}'" {input.script} {log}'''

## rna velocity
rule velocyto:
	input:
	  outputdir + "Rout/pkginstall_state.txt",
		bam_files = expand(outputdir + "HISAT2/{sample}/{sample}_Aligned.sortedByCoord.out.bam", sample = samples.name.values.tolist()),
	output:
		loom_file = outputdir + "velocyto/" + os.path.basename(proj_dir) + ".loom"
	log:
		outputdir + "Rout/velocyto.Rout"
	benchmark:
		outputdir + "benchmarks/velocyto.txt"
	params:
	  repeat_mask = config["repeat_mask"],
	  proj_name = os.path.basename(proj_dir),
	  gtf = config["gtf"],
	  loom_dir = outputdir + "velocyto/"
	conda:
		Renv
	shell:
		"velocyto run-smartseq2 -o {params.loom_dir} -m {params.repeat_mask} -e {params.proj_name} {input.bam_files} {params.gtf}"
		

## rna velocity on a seurat object
rule velocyto_seurat:
  input:
  	loom_file = outputdir + "velocyto/" + os.path.basename(proj_dir) + ".loom",
  	seu_file = outputdir + "seurat/unfiltered_seu.rds",
  	script = "scripts/compute_velocity.R"
  output:
    outputdir + "velocyto/unfiltered_seu.rds"
	log:
		outputdir + "Rout/velocyto.Rout"
	benchmark:
		outputdir + "benchmarks/veloctyo_seurat.txt"
	params:
		organism = config["organism"]
	conda:
		Renv
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args loom_path='{input.loom_file}' proj_dir='{proj_dir}' outrds='{output}' organism='{params.organism}'" {input.script} {log}'''


## ------------------------------------------------------------------------------------ ##
## Transposable element quantification
## ------------------------------------------------------------------------------------ ##

## scTE
rule scTE:
	input:
		bam_files = expand(outputdir + "HISAT2/{sample}/{sample}_Aligned.sortedByCoord.out.bam", sample = samples.name.values.tolist()),
	output:
		te_csv = outputdir + "scTE/" + os.path.basename(proj_dir) + "_tes.csv"
	log:
		outputdir + "scTE/scTE.out"
	benchmark:
		outputdir + "benchmarks/scTE.txt"
	params:
	  te_build = config["te_build"],
	  te_index = config["te_index"]
	conda:
		"envs/environment.yaml"
	shell:
		"scTE -i {input.bam_files} -o {output} -g {params.te_build} -x {params.te_index} -CB False -UMI False"

## ------------------------------------------------------------------------------------ ##
## numbat
## ------------------------------------------------------------------------------------ ##

# include: "rules/numbat.smk"
rule pileup_and_phasing:
	input:
		bam = outputdir + "cellranger/{sample}/outs/possorted_genome_bam.bam",
		barcodes = outputdir + "cellranger/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
	output:
		allele_table = outputdir + "numbat/{sample}_allele_counts.tsv.gz"
	log:
		outputdir + "logs/pileup_and_phase_{sample}.log"
	benchmark:
		outputdir + "benchmarks/pileup_and_phase_{sample}.txt"
	threads:
		config["ncores"]
	params:
		gmap = config["gmap"],
		snpvcf = config["snpvcf"],
		paneldir = config["paneldir"],
		numbatdir = outputdir + "numbat",
		outdir = outputdir + "numbat/{sample}",
		script = "scripts/pileup_and_phase.R"
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args eagle='eagle' smartseq='FALSE' cellTAG='CB' UMItag='Auto' label='{wildcards.sample}' samples='{wildcards.sample}' bams='{input.bam}' barcodes='{input.barcodes}' gmap='{params.gmap}' snpvcf='{params.snpvcf}' paneldir='{params.paneldir}' numbatdir='{params.numbatdir}' outdir='{params.outdir}' ncores='{threads}'" {params.script} {log}'''

rule numbat:
	input:
		allele_table = outputdir + "numbat/{sample}_allele_counts.tsv.gz",
		matrix_dir = outputdir + "cellranger/{sample}/outs/filtered_feature_bc_matrix"
		# matrix_dir = datadir + "{sample}"
	output:
		outputdir + "numbat/{sample}/done.txt"
	log:
		outputdir + "logs/numbat_{sample}.log"
	benchmark:
		outputdir + "benchmarks/numbat_{sample}.txt"
	threads:
		config["ncores"]
	params:
		numbatdir = outputdir + "numbat/",
		max_entropy = config["max_entropy"],
		min_LLR = config["min_LLR"],
		max_iter = config["max_iter"],
		read_prop = config["read_prop"],
		tau = config["tau"],
		numbat_t = config["numbat_t"],
		prof= outputdir + "numbat/{sample}/log.prof",
		normal_reference_mat = config["normal_reference_mat"],
		script = "scripts/run_numbat.R"
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args normal_reference_mat='{params.normal_reference_mat}' tau='{params.tau}' read_prop='{params.read_prop}' max_iter='{params.max_iter}' min_LLR='{params.min_LLR}' t='{params.numbat_t}' max_entropy='{params.max_entropy}' allele_df='{input.allele_table}' matrix_dir='{input.matrix_dir}' out_dir='{params.numbatdir}/{wildcards.sample}' ncores='{threads}' rprof_out='{params.prof}'" {params.script} {log}'''

rule pagoda_prep:
	input:
		matrix_dir = datadir + "{sample}",
		script = "scripts/prep_pagoda.R"
	output:
		prepped_pagoda = outputdir + "pagoda/{sample}.rds"
		# aneuploid_pagoda = outputdir + "pagoda/{sample}_aneuploid.rds",
		# webapp = outputdir + "numbat/{sample}_app.rds",
		# aneuploid_webapp = outputdir + "numbat/{sample}_aneuploid_app.rds"
	log:
		outputdir + "logs/pagoda_{sample}.log"
	benchmark:
		outputdir + "benchmarks/pagoda_{sample}.txt"
	threads:
		config["ncores"]
	params:
		numbat_outdir = outputdir + "numbat/{sample}/"
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args ncores='{threads}' numbat_outdir='{params.numbat_outdir}' matrix_dir='{input.matrix_dir}' outrds='{output.prepped_pagoda}'" {input.script} {log}'''

rule numbat_report:
	input:
		clone_file = outputdir + "numbat/{sample}/" + clones_filename,
		prepped_pagoda = outputdir + "pagoda/{sample}.rds",
		script = "scripts/render_rmarkdown.R"
	output:
		report = outputdir + "numbat/{sample}/{sample}_report.html"
	log:
		outputdir + "logs/report_{sample}.log"
	benchmark:
		outputdir + "benchmarks/report_{sample}.txt"
	threads:
		config["ncores"]
	params:
		numbat_dir = outputdir + "numbat/{sample}/",
		rmarkdown_file = os.getcwd() + "/scripts/plot_numbat_output.Rmd"
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args rmarkdown_file='{params.rmarkdown_file}' clone_file='{input.clone_file}' report='{output.report}' prepped_pagoda='{input.prepped_pagoda}' numbat_dir='{params.numbat_dir}'" {input.script} {log}'''

rule run_infercnv:
	input:
		script = "scripts/run_infercnv.R",
		matrix_dir = outputdir + "cellranger/{sample}/outs/filtered_feature_bc_matrix"
		# matrix_dir = datadir + "{sample}"
	output:
		annotations_path = outputdir + "infercnv/{sample}/annotations.tsv",
		infercnv_plot = outputdir + "infercnv/{sample}/infercnv.png",
	log:
		outputdir + "logs/infercnv_{sample}.log"
	benchmark:
		outputdir + "benchmarks/infercnv_{sample}.txt"
	threads:
		config["ncores"]
	params:
		normal_reference_mat = config["normal_reference_mat"],
		out_dir = outputdir + "infercnv/{sample}/"
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args matrix_dir='{input.matrix_dir}' threads='{threads}' normal_reference_mat='{params.normal_reference_mat}' annotations_path='{output.annotations_path}' out_dir='{params.out_dir}' " {input.script} {log}'''


## ------------------------------------------------------------------------------------ ##
## Input variable check
## ------------------------------------------------------------------------------------ ##
def geneset_param(wildcards):
	if config["run_camera"]:
                gs = config["genesets"].replace(" ", "") if config["genesets"] is not None else "NOTDEFINED"
		return "genesets='" + gs + "'"
	else:
		return ""


## check design matrix and contrasts
rule checkinputs:
    input:
        "config.yaml",
        script = "scripts/check_input.R"
    output:
        outputdir + "Rout/check_input.txt"
    log:
        outputdir + "Rout/check_input.Rout"
    benchmark:
    	outputdir + "benchmarks/check_input.txt"
    params:
        gtf = config["gtf"],
        genome = config["genome"],
        txome = config["txome"],
        fastqdir = config["FASTQ"],
        metatxt = config["metatxt"],
        design = config["design"].replace(" ", "") if config["design"] is not None else "NOTDEFINED",
        contrast = config["contrast"].replace(" ", "") if config["contrast"] is not None else "NOTDEFINED",
        annotation = config["annotation"].replace(" ", "") if config["annotation"] is not None else "NOTDEFINED",
        genesets = geneset_param,
        fqsuffix = str(config["fqsuffix"]),
        fqext1 = str(config["fqext1"]),
        fqext2 = str(config["fqext2"]),
        run_camera = str(config["run_camera"]),
        organism = config["organism"]    
    conda:
	    Renv
    shell:
        '''{Rbin} CMD BATCH --no-restore --no-save "--args metafile='{params.metatxt}' design='{params.design}' contrast='{params.contrast}' outFile='{output}' gtf='{params.gtf}' genome='{params.genome}' fastqdir='{params.fastqdir}' fqsuffix='{params.fqsuffix}' fqext1='{params.fqext1}' fqext2='{params.fqext2}' txome='{params.txome}' run_camera='{params.run_camera}' organism='{params.organism}' {params.genesets} annotation='{params.annotation}'" {input.script} {log};
        cat {output}
        '''
       

## ------------------------------------------------------------------------------------ ##
## Differential expression
## ------------------------------------------------------------------------------------ ##
rule edgeR:
	input:
		outputdir + "Rout/pkginstall_state.txt",
		rds = outputdir + "seurat/unfiltered_seu.rds",
		script = "scripts/run_render.R",
		template = "scripts/edgeR_dge.Rmd"
	output:
		html = outputdir + "outputR/edgeR_dge.html",
		rds = outputdir + "outputR/edgeR_dge.rds"
	params:
		directory = outputdir + "outputR",
		organism = config["organism"],        
                design = config["design"].replace(" ", "") if config["design"] is not None else "",
                contrast = config["contrast"].replace(" ", "") if config["contrast"] is not None else "",
		genesets = geneset_param
	log:
		outputdir + "Rout/run_dge_edgeR.Rout"
	benchmark:
		outputdir + "benchmarks/run_dge_edgeR.txt"
	conda:
		Renv
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args se='{input.rds}' organism='{params.organism}' design='{params.design}' contrast='{params.contrast}' {params.genesets} rmdtemplate='{input.template}' outputdir='{params.directory}' outputfile='edgeR_dge.html'" {input.script} {log}'''

## ------------------------------------------------------------------------------------ ##
## Differential transcript usage
## ------------------------------------------------------------------------------------ ##
## DRIMSeq
rule DRIMSeq:
	input:
	  outputdir + "Rout/pkginstall_state.txt",
		rds = outputdir + "seurat/unfiltered_seu.rds",
		script = "scripts/run_render.R",
		template = "scripts/DRIMSeq_dtu.Rmd"
	output:
		html = outputdir + "outputR/DRIMSeq_dtu.html",
		rds = outputdir + "outputR/DRIMSeq_dtu.rds"
	params:
		directory = outputdir + "outputR",
		organism = config["organism"],
		ncores = config["ncores"], design = config["design"].replace(" ", "") if config["design"] is not None else "", contrast = config["contrast"].replace(" ", "") if config["contrast"] is not None else ""
	log:
		outputdir + "Rout/run_dtu_drimseq.Rout"
	benchmark:
		outputdir + "benchmarks/run_dtu_drimseq.txt"
	conda:
		Renv
	threads:
		config["ncores"]
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args se='{input.rds}' design='{params.design}' contrast='{params.contrast}' ncores='{params.ncores}' rmdtemplate='{input.template}' outputdir='{params.directory}' outputfile='DRIMSeq_dtu.html'" {input.script} {log}'''

## ------------------------------------------------------------------------------------ ##
## shiny app
## ------------------------------------------------------------------------------------ ##
def shiny_input(wildcards):
	input = [outputdir + "Rout/pkginstall_state.txt"]
	if config["run_STAR"]:
		input.extend(expand(outputdir + "STAR/{sample}_Aligned.sortedByCoord.out.bw", sample = samples.name.values.tolist()))
	if config["run_HISAT2"]:
		input.extend(expand(outputdir + "HISAT2bigwig/{sample}_Aligned.sortedByCoord.out.all.bw", sample = samples.name.values.tolist()))
	return input

def shiny_params(wildcards):
	param = ["outputdir='" + outputdir + "outputR'"]
	if config["run_STAR"]:
		param.append("bigwigdir='" + outputdir + "STAR'")
	if config["run_HISAT2"]:
		param.append("bigwigdir='" + outputdir + "HISAT2bigwig'")
	return param

## shiny
rule shiny:
	input:
		shiny_input,
		rds = outputdir + "outputR/DRIMSeq_dtu.rds" if config["run_DRIMSeq"]
			else outputdir + "outputR/edgeR_dge.rds",
		script = "scripts/run_render.R",
		gtf = config["gtf"],
		template = "scripts/prepare_shiny.Rmd"
	output:
		html = outputdir + "outputR/prepare_shiny.html",
		rds = outputdir + "outputR/shiny_sce.rds"
	params:
		p = shiny_params
	log:
		outputdir + "Rout/prepare_shiny.Rout"
	benchmark:
		outputdir + "benchmarks/prepare_shiny.txt"
	conda:
		Renv
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args se='{input.rds}' gtffile='{input.gtf}' rmdtemplate='{input.template}' outputfile='prepare_shiny.html' {params.p}" {input.script} {log}'''
		
## ------------------------------------------------------------------------------------ ##
## dbtss coverage mapping
## ------------------------------------------------------------------------------------ ##
## compute coverage from dbtss
rule dbtss:
	input:
		sorted_bam = outputdir + "HISAT2/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		coverage_txt = outputdir + "dbtss_coverage/{sample}_dbtss_coverage_over_10.txt"
	threads:
		config["ncores"]
	log:
		outputdir + "logs/dbtss_{sample}.log"
	benchmark:
		outputdir + "benchmarks/dbtss_{sample}.txt"
	params:
		dbtss_bed = config["dbtss_bed"],
		genome_file = config["chrom_sizes"]
	conda:
		"envs/environment.yaml"
	shell:
	  "bedtools coverage -sorted -s -g {params.genome_file} -a {params.dbtss_bed} -b {input.sorted_bam} | awk '$4 > 10 {{print}}' > {output.coverage_txt}"

## ------------------------------------------------------------------------------------ ##
## configure jbrowse
## ------------------------------------------------------------------------------------ ##

include: "rules/jbrowse.smk"

## ------------------------------------------------------------------------------------ ##
## Success and failure messages
## ------------------------------------------------------------------------------------ ##
onsuccess:
	print("Success! The Snakemake workflow is completed.")

onerror:
	print("Error! The Snakemake workflow aborted.")

