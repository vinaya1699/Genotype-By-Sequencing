# Example Usage : 
# python Python_GBS.py -org Phaseolus_vulgaris -t 20 -d /path/to/working_directory 
print(r"""
   ____   ____   ____
  / ___| | __ ) | ___|
 | |  _  |  _ \ |___ \
 | |_| | | |_) | ___) |
  \____| |____/ |____/

        Genotype By Sequencing (GBS) Automation Pipeline
        : Vinaya Kadam
""")

print(r"""
      # Example Usage : 
# python Python_GBS.py -org Phaseolus_vulgaris -t 20 -d /path/to/working_directory 
""")

import pandas as pd
import numpy as np
import os
import sys
import glob
import subprocess
import argparse
import datetime
import time
import shutil
import logging


script_start = time.time()
print(f"Script started at: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# Argument parser for command-line options
parser = argparse.ArgumentParser(description="GBS Automation Script")
parser.add_argument('-d', '--Working_Directory', type=str,required=True, help='Input working directory containing raw FASTQ files')
parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads')
parser.add_argument('-org', '--organism', type=str, required=True, help='Organism name (reference fasta prefix)')
parser.add_argument('--annotation-file', type=str, default=None, help='Explicit annotation file ( GTF/GFF/GFF3 )')
parser.add_argument('--annotation-format', type=str, choices=['gtf', 'gff', 'gff3'], default=None, help='Annotation format for snpEff build; autodetect if not set')
# parser.add_argument('--qual', type=float, default=30, help='Minimum QUAL value')
# parser.add_argument('--min_dp', type=int, default=5, help='Minimum INFO/DP value')
# parser.add_argument('--f_missing', type=float, default=0.8, help='Maximum F_MISSING value')
# parser.add_argument('--min_af', type=float, default=0.05, help='Minimum allele frequency')

args = parser.parse_args()
threads = args.threads
organism = args.organism

# Change to input directory
os.chdir(args.Working_Directory)
print(f"Changed working directory to: {args.Working_Directory}")


ref_dir = "0_Reference_Genome"
ref_fasta = f"{ref_dir}/{organism}.fasta"
picard = "/Analysis3/Vinaya/picard.jar"
gatk = "/apps/gatk-4.2.6.1/gatk"  
gvcf_dir = "4_Variant_Calling"
snpEff = "/apps/snpEff5.0/snpEff.jar"
snpEff_data="/mnt/tnas/SnpEff5.0_Database/data"
snpEff_config = "/mnt/tnas/SnpEff5.0_Database/"

Input_Raw_files = glob.glob("1_RawData/*_R1.fastq.gz")
print(Input_Raw_files)

# Convert list to pandas Series, then use replace and unique
files_series = pd.Series(Input_Raw_files)
# Remove directory path and suffix to get sample names
samples = files_series.str.replace(r'.*/|_R1\.fastq\.gz$', '', regex=True).unique()

print("Sample names extracted:", samples)

# Create directories if they don't exist
folders = [
    "2_Clean_data",
    "3_Alignment",
    "1_RawData/Fastqc_Output",
    "2_Clean_data/Fastqc_Output",
    "4_Variant_Calling"
]
for folder in folders:
    os.makedirs(folder, exist_ok=True)
    os.chmod(folder, 0o777)

print("Directories created and permissions set.")

print(r"""

""")


print("############################### fastqc ##################################")

print(f"Using {threads} threads for FastQC analysis.")
# Run FastQC for each sample
for sample in samples:
    r1 = f"1_RawData/{sample}_R1.fastq.gz"
    r2 = f"1_RawData/{sample}_R2.fastq.gz"
    log = f"1_RawData/Fastqc_Output/{sample}_fastqc.log"
    cmd = [
        "fastqc",
        "-t", str(threads),
        r1, r2,
        "-o", "1_RawData/Fastqc_Output/"
    ]
    with open(log, "w") as logfile:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            print(line, end='')      # Print to console
            logfile.write(line)      # Write to log file
            print("FastQC analysis completed for all samples.")
        process.wait()


print("############################### fastp ##################################")

print(f"Using {threads} threads for fastp analysis.")
for sample in samples:
    r1 = f"1_RawData/{sample}_R1.fastq.gz"
    r2 = f"1_RawData/{sample}_R2.fastq.gz"
    r1_clean = f"2_Clean_data/{sample}_R1.fq.gz"
    r2_clean = f"2_Clean_data/{sample}_R2.fq.gz"
    html = f"2_Clean_data/{sample}_fastp.html"
    json = f"2_Clean_data/{sample}_fastp.json"
    log = f"2_Clean_data/{sample}_fastp.log"
    cmd = [
        "fastp",
        "-i", r1,
        "-I", r2,
        "-o", r1_clean,
        "-O", r2_clean,
        "-h", html,
        "-j", json,
        "-w", str(threads),
        "-c",
        "--detect_adapter_for_pe"
    ]
    with open(log, "w") as logfile:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            print(line, end='')      # Print to console
            logfile.write(line)      # Write to log file
            print("fastp analysis completed for all samples.")
        process.wait()

print("############################### Clean Data QC ##################################")

print(f"Using {threads} threads for FastQC analysis.")
# Run FastQC for each sample
for sample in samples:
    r1 = f"2_Clean_data/{sample}_R1.fq.gz"
    r2 = f"2_Clean_data/{sample}_R2.fq.gz"
    log = f"2_Clean_data/Fastqc_Output/{sample}_fastqc.log"
    cmd = [
        "fastqc",
        "-t", str(threads),
        r1, r2,
        "-o", "2_Clean_data/Fastqc_Output/"
    ]
    with open(log, "w") as logfile:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            print(line, end='')      # Print to console
            logfile.write(line)      # Write to log file
            print("FastQC analysis completed for all clean data samples.")
        process.wait()

print("############################### Indexing Ref. Genome ##################################")

# List of required index extensions
index_extensions = [".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"]

# Build full file paths for the required index files
missing_files = []
for ext in index_extensions:
    file_path = os.path.join(ref_dir, os.path.splitext(os.path.basename(ref_fasta))[0] + ext)
    if not os.path.exists(file_path):
        missing_files.append(file_path)

if missing_files:
    print(f"Missing index files: {missing_files}")
    print("Indexing reference genome...")
    cmd = ["bwa", "index", ref_fasta]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    for line in process.stdout:
        print(line, end='')
    process.wait()
    print("Reference genome indexing completed.")
else:
    print("All index files found. Reference genome is already indexed.")

print("############################### Alignment & MarkDuplicates ##################################")


alignment_dir = "3_Alignment"

for sample in samples:
    r1 = f"2_Clean_data/{sample}_R1.fq.gz"
    r2 = f"2_Clean_data/{sample}_R2.fq.gz"
    bam = f"{alignment_dir}/{sample}.bam"
    marked_bam = f"{alignment_dir}/{sample}_sorted_rdgrp_dup_marked.bam"
    metrics = f"{alignment_dir}/{sample}_dup_mark-metrics.txt"
    read_group = f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"

    # 1. Mapping and sorting
    bwa_cmd = [
        "bwa", "mem",
        "-R", read_group,
        "-M",
        "-t", str(threads),
        ref_fasta,
        r1, r2
    ]
    samtools_cmd = [
        "samtools", "sort",
        "-@", str(threads),
        "-o", bam,
        "-"
    ]
    print(f"Running mapping and sorting for sample: {sample}")
    bwa = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    samtools = subprocess.Popen(samtools_cmd, stdin=bwa.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    bwa.stdout.close()
    samtools.communicate()
    bwa.wait()
    samtools.wait()

    # 2. Mark duplicates (final BAM)
    print(f"Running MarkDuplicates for sample: {sample}")
    markdup_cmd = [
        "java", "-Xmx16g", "-jar", picard, "MarkDuplicates",
        f"I={bam}",
        f"O={marked_bam}",
        f"M={metrics}",
        "QUIET=TRUE"
    ]
    subprocess.run(markdup_cmd, check=True)
    print(f"Final BAM generated: {marked_bam}")

print("All samples processed. Only one final BAM per sample is generated.")

print("############################### GATK HaplotypeCaller ##################################")


# Ensure output directory exists
os.makedirs(gvcf_dir, exist_ok=True)

# Index reference FASTA with samtools faidx if .fai does not exist
fai_file = f"{ref_fasta}.fai"
if not os.path.exists(fai_file):
    print(f"Indexing reference genome with samtools faidx: {ref_fasta}")
    cmd = ["samtools", "faidx", ref_fasta]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    for line in process.stdout:
        print(line, end='')
    process.wait()
    print("Reference genome indexed with samtools faidx.")
else:
    print(f"Reference genome index (.fai) already exists: {fai_file}")

# reference dictionary if .dict does not exist
dict_file = os.path.splitext(ref_fasta)[0] + ".dict"

if not os.path.exists(dict_file):
    print(f"Generating dictionary file: {dict_file}")
    cmd = [
        "java", "-jar", picard, "CreateSequenceDictionary",
        f"R={ref_fasta}",
        f"O={dict_file}"
    ]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    for line in process.stdout:
        print(line, end='')
    process.wait()
    print("Reference dictionary was created using Picard.")
else:
    print(f"Reference dictionary (.dict) already exists: {dict_file}")

for sample in samples:
    marked_bam = f"3_Alignment/{sample}_sorted_rdgrp_dup_marked.bam"
    # bam_index = f"{marked_bam}.bai"
    bam_index = f"{marked_bam}.csi"
    # Index BAM if not present
    if not os.path.exists(bam_index):
        print(f"Indexing BAM file for sample {sample}: {marked_bam}")
        # cmd = ["samtools", "index", marked_bam]
        cmd = ["samtools", "index",  "-c",marked_bam]
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            print(line, end='')
        process.wait()
        print(f"BAM index created: {bam_index}")
    else:
        print(f"BAM index already exists: {bam_index}")

    gvcf = f"{gvcf_dir}/{sample}.gvcf"
    log = f"{gvcf_dir}/{sample}.gvcf.log"
    cmd = [
        gatk, "HaplotypeCaller",
        "-I", marked_bam,
        "-R", ref_fasta,
        "-O", gvcf,
        "-ERC", "GVCF"
    ]
    print(f"Running GATK HaplotypeCaller for sample: {sample}")
    with open(log, "w") as logfile:
        process = subprocess.Popen(cmd, stdout=logfile, stderr=subprocess.STDOUT, text=True)
        process.wait()
    print(f"GVCF generated: {gvcf}")

print("GATK HaplotypeCaller completed for all samples.")

print("############################### Joint Genotype Calling ##################################")

# 1. Create interval.list (list of all contigs/chromosomes)
interval_list = "interval.list"
if not os.path.exists(interval_list):
    contigs = set()
    with open(ref_fasta, "r") as fasta:
        for line in fasta:
            if line.startswith(">"):
                # Skip contigs with mitochondrion or chloroplast in the header
                if "mitochondrion" in line.lower() or "chloroplast" in line.lower():
                    continue
                contig = line[1:].strip().split()[0]  # Take everything before first space
                contigs.add(contig)
    with open(interval_list, "w") as out:
        for contig in sorted(contigs):
            out.write(contig + "\n")
    print(f"interval.list created: {interval_list}")
else:
    print(f"interval.list already exists: {interval_list}")

# 2. Create samplenames.txt (sample_name <tab> gvcf_path)
sample_map = "samplenames.txt"
with open(sample_map, "w") as f:
    for sample in samples:
        gvcf = f"{gvcf_dir}/{sample}.gvcf"
        f.write(f"{sample}\t{os.path.abspath(gvcf)}\n")
print(f"samplenames.txt created: {sample_map}")

# 3. Run GenomicsDBImport
genomicsdb_dir = "gatk_db"
genomicsdb_cmd = [
    gatk, "GenomicsDBImport",
    "--genomicsdb-workspace-path", genomicsdb_dir,
    "--batch-size", "40",
    "--reader-threads", str(threads),
    "--sample-name-map", sample_map,
    "-L", interval_list,
    "--consolidate", "false"
]
print("Running GenomicsDBImport...")
process = subprocess.Popen(genomicsdb_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
for line in process.stdout:
    print(line, end='')
process.wait()
print("GenomicsDBImport completed.")

# 4. Run GenotypeGVCFs
jointvcf = "joint_genotyped.vcf.gz"
genotypecmd = [
    gatk, "--java-options", "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true",
    "GenotypeGVCFs", "-R", ref_fasta,
    "-V", f"gendb://{genomicsdb_dir}", "-O", jointvcf
]
print("Running GenotypeGVCFs...")

# CHANGE: Set environment variable for TileDB file locking
env = os.environ.copy()
env["TILEDB_DISABLE_FILE_LOCKING"] = "1"

process = subprocess.Popen(genotypecmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)
for line in process.stdout:
    print(line, end="")
process.wait()
print(f"Joint genotyping completed. Output: {jointvcf}")

print(f"Script ended at: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Total script duration: {time.time() - script_start:.2f} seconds")

# After joint genotyping, add:
print("############################### SNP Extraction and Filtering ##################################")

jointvcf = "joint_genotyped.vcf.gz"
snps_vcf = "snps_only.vcf"
filtered_vcf = "Filtered_SNP.vcf.gz"

# 1. Extract SNPs only
snps_cmd = [
    gatk, "SelectVariants", "-V", jointvcf, "-select-type", "SNP", "-O", snps_vcf
]
print("Extracting SNPs only...")
process = subprocess.Popen(snps_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
for line in process.stdout:
    print(line, end="")
process.wait()
print(f"SNPs extracted to {snps_vcf}")

# 2. Filter SNPs with user parameters
# filter_expr = f"QUAL <= {args.qual} && INFO/DP <= {args.min_dp} && F_MISSING=={args.f_missing}"
filter_cmd = [
    gatk, "VariantFiltration", "-V", snps_vcf,
        "-filter", "QD < 2.0", "--filter-name", "QD2",
        "-filter", "QUAL < 30.0", "--filter-name", "QUAL30",
        "-filter", "SOR > 3.0", "--filter-name", "SOR3",
        "-filter", "FS > 60.0", "--filter-name", "FS60",
        "-filter", "MQ < 40.0", "--filter-name", "MQ40",
        "-filter", "MQRankSum < -12.5", "--filter-name", "MQRankSum-12.5",
        "-filter", "ReadPosRankSum < -8.0", "--filter-name", "ReadPosRankSum-8",
        "-O", filtered_vcf
]
print("Filtering SNPs with user-defined parameters...")
process = subprocess.Popen(filter_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
for line in process.stdout:
    print(line, end="")
process.wait()
print(f"Filtered SNPs written to {filtered_vcf}")

print("############################### SNP Annotation ##################################")

organism_snpEff_dir = os.path.join(snpEff_data, organism)

# Create organism directory if it does not exist
if not os.path.exists(organism_snpEff_dir):
    print(f"Creating directory for snpEff organism data at {organism_snpEff_dir}")
    os.makedirs(organism_snpEff_dir, exist_ok=True)

# Copy fasta to snpEff data directory, rename to sequences.fa
ref_fasta_dest = os.path.join(organism_snpEff_dir, "sequences.fa")
print(f"Copying reference fasta {ref_fasta} to {ref_fasta_dest}")
shutil.copy2(ref_fasta, ref_fasta_dest)

# Copy annotation, support GTF/GFF/GFF3
annotation_file = args.annotation_file
annotation_format = args.annotation_format
if annotation_file:
    if not os.path.isabs(annotation_file):
        annotation_file = os.path.join(args.Working_Directory, annotation_file)
    if not os.path.exists(annotation_file):
        raise FileNotFoundError(f"Annotation file not found: {annotation_file}")
    _, ext = os.path.splitext(annotation_file)
    ext = ext.lower()
    if ext == ".gtf":
        annotation_format = "gtf"
    elif ext in [".gff3", ".gff"]:
        annotation_format = "gff3" if ext == ".gff3" else "gff"
    else:
        raise ValueError(f"Unsupported annotation file extension: {ext}")
else:
    if annotation_format:
        annotation_format = annotation_format.lower()
        if annotation_format not in ["gtf", "gff", "gff3"]:
            raise ValueError(f"Unsupported annotation format: {annotation_format}")
        candidate = f"{os.path.splitext(ref_fasta)[0]}.{annotation_format}"
        if annotation_format == "gff3" and not os.path.exists(candidate):
            candidate = f"{os.path.splitext(ref_fasta)[0]}.gff"  # fallback
        if not os.path.exists(candidate):
            raise FileNotFoundError(f"Annotation file not found for --annotation-format {annotation_format}: {candidate}")
        annotation_file = candidate
    else:
        preferred = [".gtf", ".gff3", ".gff"]
        for suf in preferred:
            candidate = f"{os.path.splitext(ref_fasta)[0]}{suf}"
            if os.path.exists(candidate):
                annotation_file = candidate
                break
        if not annotation_file:
            raise FileNotFoundError("No annotation file found (searched .gtf, .gff3, .gff)")
        _, ext = os.path.splitext(annotation_file)
        if ext.lower() == ".gtf":
            annotation_format = "gtf"
        elif ext.lower() == ".gff3":
            annotation_format = "gff3"
        else:
            annotation_format = "gff"

if annotation_format == "gtf":
    ann_dest = os.path.join(organism_snpEff_dir, "genes.gtf")
    snpEff_build_mode = "-gtf22"
elif annotation_format in ["gff", "gff3"]:
    ann_dest = os.path.join(organism_snpEff_dir, "genes.gff")
    snpEff_build_mode = "-gff3"
else:
    raise ValueError(f"Unsupported annotation format after detection: {annotation_format}")

print(f"Copying annotation {annotation_file} to {ann_dest} (format {annotation_format})")
shutil.copy2(annotation_file, ann_dest)

print("Reference genome and annotation copied and renamed for snpEff.")

# Build snpEff database
# Check if the database already exists by looking for the 'snpEffectPredictor.bin' file
snpEff_db_file = os.path.join(organism_snpEff_dir, "snpEffectPredictor.bin")
if not os.path.exists(snpEff_db_file):
    print(f"Building snpEff database for organism: {organism} with {annotation_format} annotations")
    build_cmd = [
        "java", "-jar", snpEff,
        "build", snpEff_build_mode, "-v", "-c", snpEff_config, organism
    ]
    process = subprocess.Popen(build_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    for line in process.stdout:
        print(line, end="")
    process.wait()
    print("snpEff database build completed.")
else:
    print(f"snpEff database already exists for organism: {organism}")   

# Annotate filtered VCF
filtered_vcf = "Filtered_SNP.vcf.gz"
annotated_vcf = "Annotated_SNP.vcf"
annotate_cmd = [
    "java", "-jar", snpEff,
    "-c", snpEff_config,
    "-v", organism,     
    filtered_vcf
]
print(f"Annotating filtered VCF: {filtered_vcf}")
with open(annotated_vcf, "w") as outfile:
    process = subprocess.Popen(annotate_cmd, stdout=outfile, stderr=subprocess.STDOUT, text=True)
    process.wait()

print(f"Annotation completed. Output written to {annotated_vcf}")      

#####################################################################################################
print(f"Script ended at: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Total script duration: {time.time() - script_start:.2f} seconds")
