## This script prepares target PAH-region genotype data for haplotype and IBD analysis.
## It first restricts the target VCF to high-quality biallelic SNPs without missing genotypes in the PAH locus.
## The same genomic region is extracted from a phased 1000 Genomes chromosome 12 reference panel.
## Only variants shared between the target samples and the reference panel are retained.
## The target haplotypes are then phased with SHAPEIT5 using the matched reference panel and genetic map.
## Finally, the phased VCF and PLINK-formatted genetic map are used as input for hap-IBD to infer shared IBD segments around PAH.

## ============================================================
## Target sample VCF files and output files
## ============================================================

# Raw VCF containing PAH-region variants for the target individuals
RAW="data/PAH_variants.vcf.gz"

# Cleaned target VCF after restricting to region, biallelic SNPs, and no missing genotypes
CLEAN="data/PAH_variants_clean.vcf.gz"

# Target VCF restricted to variants overlapping the 1000 Genomes reference panel
CLEANo="data/PAH_variants_clean_overlap.vcf.gz"

# Output from SHAPEIT5 phasing in BCF format
PHASED="data/PAH_variants_clean_overlap_phased.bcf"

# Same phased target data converted to compressed VCF format
PHASEDv="data/PAH_variants_clean_overlap_phased.vcf.gz"

## ============================================================
## Reference panel and genetic map
## ============================================================

# 1000 Genomes chr12 reference panel, already phased
REF=/projects/popgen/data/genomes/1000Genotypes/vcf/vcfGRCh38/CCDG_14151_B01_GRM_WGS_2020-08-05_chr12.filtered.shapeit2-duohmm-phased.vcf.gz

# Reference panel restricted to the PAH region and biallelic SNPs
REFREG=data/1000g.chr12.region.snps.vcf.gz

# Reference panel restricted to variants overlapping the target data
REFo=data/1000g.chr12.region_overlab.snps.vcf.gz

# SHAPEIT5 genetic map in gmap format
MAP=/maps/projects/aa-AUDIT/people/bcn627/PAH/prog/shapeit5/resources/maps/b38/chr12.b38.gmap.gz

# Same genetic map converted to PLINK-style map format for downstream IBD tools
MAPplink=/maps/projects/aa-AUDIT/people/bcn627/PAH/data/chr12.b38.plink.map

# Convert the genetic map to PLINK map format:
# columns: chromosome, marker ID, genetic position, physical position
# The first line of the gmap file is skipped as header
awk 'NR==1 {next} {print "chr12", "chr12:"$1, $3, $1}' <(zcat $MAP) > $MAPplink

## ============================================================
## IBD output prefix
## ============================================================

# Prefix for hap-IBD output files
IBDout="results/PAH_ibd"

## ============================================================
## Region of interest: PAH locus plus surrounding region
## ============================================================

CHR=chr12
START=98376431
END=107419587

# Region string used by bcftools and SHAPEIT5
REGION=${CHR}:${START}-${END}

## ============================================================
## Load required software modules
## ============================================================

module load gsl/2.5 perl bcftools htslib

## ============================================================
## Inspect raw target VCF
## ============================================================

# Index the raw VCF using a tabix index
bcftools index -t $RAW

# Print summary statistics for the raw VCF
# SN lines contain general numbers such as number of samples and records
bcftools stats $RAW | grep ^SN 


## ============================================================
## Clean target VCF
## ============================================================

# Restrict target data to:
# - the PAH region
# - biallelic sites only: -m2 -M2
# - SNPs only: -v snps
# N_MISSING<=0 means all samples must have called genotypes at the site.
bcftools view -i 'N_MISSING<=0' -r $REGION -m2 -M2 -v snps $RAW -Oz -o $CLEAN

# Index the cleaned target VCF
bcftools index -f -t $CLEAN

# Check summary statistics after filtering
bcftools stats $CLEAN | grep ^SN 

## ============================================================
## Prepare 1000 Genomes reference panel in the same region
## ============================================================

# Restrict reference panel to:
# - the PAH region
# - biallelic sites only
# - SNPs only
bcftools view -r $REGION -m2 -M2 -v snps $REF -Oz -o $REFREG

# Index the regional reference VCF
bcftools index -f -t $REFREG

# Check summary statistics for the regional reference panel
bcftools stats $REFREG | grep ^SN 


## ============================================================
## Keep only variants shared between target and reference
## ============================================================

# Find variants present in both target and reference files.
#
# -p tmp       writes output files to tmp/
# -n=2        keeps variants present in both input files
# -w1,2       writes records from both files:
#             tmp/0000.vcf = records from $CLEAN
#             tmp/0001.vcf = records from $REFREG
bcftools isec -p tmp -n=2 -w1,2 $CLEAN $REFREG

# Count overlapping variants written from the target file
bcftools view -H tmp/0000.vcf | wc -l

# Save the overlapping reference-panel variants
bcftools view -Oz -o $REFo tmp/0001.vcf

# Count overlapping variants written from the reference file
bcftools view -H tmp/0001.vcf | wc -l

# Index the overlapping reference panel
bcftools index -f -t $REFo

# Save the overlapping target variants
bcftools view -Oz -o $CLEANo tmp/0000.vcf

# Index the overlapping target VCF
bcftools index -f -t $CLEANo

# Check summary statistics for the final target VCF used for phasing
bcftools stats $CLEANo | grep ^SN 


## ============================================================
## Phase target haplotypes using SHAPEIT5 and 1000 Genomes reference
## ============================================================

# Phase the target samples using SHAPEIT5 phase_common.
#
# --input      target VCF restricted to variants overlapping the reference
# --reference 1000 Genomes reference panel restricted to the same variants
# --region    PAH region
# --map       genetic map
# --output    phased target output in BCF format
# --thread    number of CPU threads
./prog/shapeit5/phase_common/bin/phase_common \
    --input $CLEANo \
    --reference $REFo \
    --region $REGION \
    --map $MAP \
    --output $PHASED \
    --thread 8

# Index the phased BCF output
bcftools index -f "$PHASED"

# Check that the phased file header can be read
bcftools view -h "$PHASED" >/dev/null

# Inspect the first few phased variant records
bcftools view -H "$PHASED" | head


## ============================================================
## Convert phased BCF to phased VCF
## ============================================================

# Convert phased BCF to compressed VCF format
bcftools view -Oz -o $PHASEDv $PHASED

# Index the compressed phased VCF
bcftools index -f -t $PHASEDv


## ============================================================
## Target sample VCF files and output files
## ============================================================

# Raw VCF containing PAH-region variants for the target individuals
RAW="data/PAH_variants.vcf.gz"

# Cleaned target VCF after restricting to region, biallelic SNPs, and no missing genotypes
CLEAN="data/PAH_variants_clean.vcf.gz"

# Target VCF restricted to variants overlapping the 1000 Genomes reference panel
CLEANo="data/PAH_variants_clean_overlap.vcf.gz"

# Output from SHAPEIT5 phasing in BCF format
PHASED="data/PAH_variants_clean_overlap_phased.bcf"

# Same phased target data converted to compressed VCF format
PHASEDv="data/PAH_variants_clean_overlap_phased.vcf.gz"

## ============================================================
## Reference panel and genetic map
## ============================================================

# 1000 Genomes chr12 reference panel, already phased
REF=/projects/popgen/data/genomes/1000Genotypes/vcf/vcfGRCh38/CCDG_14151_B01_GRM_WGS_2020-08-05_chr12.filtered.shapeit2-duohmm-phased.vcf.gz

# Reference panel restricted to the PAH region and biallelic SNPs
REFREG=data/1000g.chr12.region.snps.vcf.gz

# Reference panel restricted to variants overlapping the target data
REFo=data/1000g.chr12.region_overlab.snps.vcf.gz

# SHAPEIT5 genetic map in gmap format
MAP=/maps/projects/aa-AUDIT/people/bcn627/PAH/prog/shapeit5/resources/maps/b38/chr12.b38.gmap.gz

# Same genetic map converted to PLINK-style map format for downstream IBD tools
MAPplink=/maps/projects/aa-AUDIT/people/bcn627/PAH/data/chr12.b38.plink.map

# Convert the genetic map to PLINK map format:
# columns: chromosome, marker ID, genetic position, physical position
# The first line of the gmap file is skipped as header
awk 'NR==1 {next} {print "chr12", "chr12:"$1, $3, $1}' <(zcat $MAP) > $MAPplink

## ============================================================
## IBD output prefix
## ============================================================

# Prefix for hap-IBD output files
IBDout="results/PAH_ibd"

## ============================================================
## Region of interest: PAH locus plus surrounding region
## ============================================================

CHR=chr12
START=98376431
END=107419587

# Region string used by bcftools and SHAPEIT5
REGION=${CHR}:${START}-${END}

## ============================================================
## Load required software modules
## ============================================================

module load gsl/2.5 perl bcftools htslib

## ============================================================
## Inspect raw target VCF
## ============================================================

# Index the raw VCF using a tabix index
bcftools index -t $RAW

# Print summary statistics for the raw VCF
# SN lines contain general numbers such as number of samples and records
bcftools stats $RAW | grep ^SN 


## ============================================================
## Clean target VCF
## ============================================================

# Restrict target data to:
# - the PAH region
# - biallelic sites only: -m2 -M2
# - SNPs only: -v snps
# N_MISSING<=0 means all samples must have called genotypes at the site.
bcftools view -i 'N_MISSING<=0' -r $REGION -m2 -M2 -v snps $RAW -Oz -o $CLEAN

# Index the cleaned target VCF
bcftools index -f -t $CLEAN

# Check summary statistics after filtering
bcftools stats $CLEAN | grep ^SN 

## ============================================================
## Prepare 1000 Genomes reference panel in the same region
## ============================================================

# Restrict reference panel to:
# - the PAH region
# - biallelic sites only
# - SNPs only
bcftools view -r $REGION -m2 -M2 -v snps $REF -Oz -o $REFREG

# Index the regional reference VCF
bcftools index -f -t $REFREG

# Check summary statistics for the regional reference panel
bcftools stats $REFREG | grep ^SN 


## ============================================================
## Keep only variants shared between target and reference
## ============================================================

# Find variants present in both target and reference files.
#
# -p tmp       writes output files to tmp/
# -n=2        keeps variants present in both input files
# -w1,2       writes records from both files:
#             tmp/0000.vcf = records from $CLEAN
#             tmp/0001.vcf = records from $REFREG
bcftools isec -p tmp -n=2 -w1,2 $CLEAN $REFREG

# Count overlapping variants written from the target file
bcftools view -H tmp/0000.vcf | wc -l

# Save the overlapping reference-panel variants
bcftools view -Oz -o $REFo tmp/0001.vcf

# Count overlapping variants written from the reference file
bcftools view -H tmp/0001.vcf | wc -l

# Index the overlapping reference panel
bcftools index -f -t $REFo

# Save the overlapping target variants
bcftools view -Oz -o $CLEANo tmp/0000.vcf

# Index the overlapping target VCF
bcftools index -f -t $CLEANo

# Check summary statistics for the final target VCF used for phasing
bcftools stats $CLEANo | grep ^SN 


## ============================================================
## Phase target haplotypes using SHAPEIT5 and 1000 Genomes reference
## ============================================================

# Phase the target samples using SHAPEIT5 phase_common.
#
# --input      target VCF restricted to variants overlapping the reference
# --reference 1000 Genomes reference panel restricted to the same variants
# --region    PAH region
# --map       genetic map
# --output    phased target output in BCF format
# --thread    number of CPU threads
./prog/shapeit5/phase_common/bin/phase_common \
    --input $CLEANo \
    --reference $REFo \
    --region $REGION \
    --map $MAP \
    --output $PHASED \
    --thread 8

# Index the phased BCF output
bcftools index -f "$PHASED"

# Check that the phased file header can be read
bcftools view -h "$PHASED" >/dev/null

# Inspect the first few phased variant records
bcftools view -H "$PHASED" | head


## ============================================================
## Convert phased BCF to phased VCF
## ============================================================

# Convert phased BCF to compressed VCF format
bcftools view -Oz -o $PHASEDv $PHASED

# Index the compressed phased VCF
bcftools index -f -t $PHASEDv

## ============================================================
## Infer IBD using hap-ibd
## ============================================================

java -Xmx16g -jar prog/hap-ibd.jar gt=$PHASEDv  map=$MAPplink out=$IBDout min-seed=1.0 min-output=1.0 min-markers=50 nthreads=8

