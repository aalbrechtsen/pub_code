# PAH phasing and IBD analysis 

This script if for preparing PAH-region genotype data for haplotype phasing and IBD analysis.
The script filters target variants, extracts the same region from a phased 1000 Genomes chromosome 12 reference panel, keeps shared variants, phases target haplotypes with SHAPEIT5, and infers IBD segments with hap-IBD.

```bash

# ============================================================
# Input and output files
# ============================================================

# Raw VCF containing PAH-region variants for the target individuals.
RAW="data/PAH_variants.vcf.gz"

# Cleaned target VCF after restricting to the PAH region, biallelic SNPs,
# and variants without missing genotypes.
CLEAN="data/PAH_variants_clean.vcf.gz"

# Target VCF restricted to variants overlapping the 1000 Genomes reference panel.
CLEANo="data/PAH_variants_clean_overlap.vcf.gz"

# Output from SHAPEIT5 phasing in BCF format.
PHASED="data/PAH_variants_clean_overlap_phased.bcf"

# Same phased target data converted to compressed VCF format.
PHASEDv="data/PAH_variants_clean_overlap_phased.vcf.gz"


# ============================================================
# Reference panel and genetic maps
# ============================================================

# Phased 1000 Genomes chromosome 12 reference panel.
REF="/projects/popgen/data/genomes/1000Genotypes/vcf/vcfGRCh38/CCDG_14151_B01_GRM_WGS_2020-08-05_chr12.filtered.shapeit2-duohmm-phased.vcf.gz"

# Reference panel restricted to the PAH region and biallelic SNPs.
REFREG="data/1000g.chr12.region.snps.vcf.gz"

# Reference panel restricted to variants overlapping the target data.
REFo="data/1000g.chr12.region_overlap.snps.vcf.gz"

# SHAPEIT5 genetic map in gmap format.
MAP="/maps/projects/aa-AUDIT/people/bcn627/PAH/prog/shapeit5/resources/maps/b38/chr12.b38.gmap.gz"

# PLINK-style genetic map for hap-IBD.
MAPplink="/maps/projects/aa-AUDIT/people/bcn627/PAH/data/chr12.b38.plink.map"


# ============================================================
# Region of interest and output prefix
# ============================================================

CHR="chr12"
START=98376431
END=107419587
REGION="${CHR}:${START}-${END}"

# Prefix for hap-IBD output files.
IBDout="results/PAH_ibd"


# ============================================================
# Setup
# ============================================================

mkdir -p data results tmp
module load gsl/2.5 perl bcftools htslib


# ============================================================
# Convert genetic map to PLINK format
# ============================================================

# hap-IBD requires columns:
# chromosome, marker ID, genetic position, physical position.
awk 'NR==1 {next} {print "chr12", "chr12:"$1, $3, $1}' <(zcat "$MAP") > "$MAPplink"


# ============================================================
# Inspect and index raw target VCF
# ============================================================

bcftools index -f -t "$RAW"
bcftools stats "$RAW" | grep '^SN'


# ============================================================
# Clean target VCF
# ============================================================

# Keep biallelic SNPs in the PAH region with no missing genotypes.
bcftools view \
    -i 'N_MISSING<=0' \
    -r "$REGION" \
    -m2 -M2 \
    -v snps \
    "$RAW" \
    -Oz \
    -o "$CLEAN"

bcftools index -f -t "$CLEAN"
bcftools stats "$CLEAN" | grep '^SN'


# ============================================================
# Extract the same region from the reference panel
# ============================================================

bcftools view \
    -r "$REGION" \
    -m2 -M2 \
    -v snps \
    "$REF" \
    -Oz \
    -o "$REFREG"

bcftools index -f -t "$REFREG"
bcftools stats "$REFREG" | grep '^SN'


# ============================================================
# Keep only variants shared between target and reference
# ============================================================

rm -rf tmp/isec_pah
mkdir -p tmp/isec_pah

# tmp/isec_pah/0000.vcf: overlapping variants from the target file.
# tmp/isec_pah/0001.vcf: overlapping variants from the reference file.
bcftools isec \
    -p tmp/isec_pah \
    -n=2 \
    -w1,2 \
    "$CLEAN" \
    "$REFREG"

# Count overlapping target and reference variants.
bcftools view -H tmp/isec_pah/0000.vcf | wc -l
bcftools view -H tmp/isec_pah/0001.vcf | wc -l

# Save and index the overlapping reference-panel variants.
bcftools view -Oz -o "$REFo" tmp/isec_pah/0001.vcf
bcftools index -f -t "$REFo"

# Save and index the overlapping target variants.
bcftools view -Oz -o "$CLEANo" tmp/isec_pah/0000.vcf
bcftools index -f -t "$CLEANo"
bcftools stats "$CLEANo" | grep '^SN'


# ============================================================
# Phase target haplotypes with SHAPEIT5
# ============================================================

./prog/shapeit5/phase_common/bin/phase_common \
    --input "$CLEANo" \
    --reference "$REFo" \
    --region "$REGION" \
    --map "$MAP" \
    --output "$PHASED" \
    --thread 8

bcftools index -f "$PHASED"
bcftools view -h "$PHASED" >/dev/null
bcftools view -H "$PHASED" | head


# ============================================================
# Convert phased BCF to compressed VCF
# ============================================================

bcftools view -Oz -o "$PHASEDv" "$PHASED"
bcftools index -f -t "$PHASEDv"


# ============================================================
# Infer IBD using hap-IBD
# ============================================================

java -Xmx16g \
    -jar prog/hap-ibd.jar \
    gt="$PHASEDv" \
    map="$MAPplink" \
    out="$IBDout" \
    min-seed=1.0 \
    min-output=1.0 \
    min-markers=50 \
    nthreads=8
```
