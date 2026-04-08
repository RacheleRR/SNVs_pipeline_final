#!/bin/bash
# =============================================================================
# VARIANT FILTERING PIPELINE

# =============================================================================
# Pipeline overview:
#   1. MPC v2 annotation 
#   2. Allele frequency filter (gnomAD AF < 1%)
#   3. Split into PTV and Missense branches
#   4. Recompress + index
#   5. ClinVar annotation
#   6. Pathogenicity filtering
#   7. Export to TSV
#   8. QC variant counts
# =============================================================================

# --- PATHS & VARIABLES -------------------------------------------------------

BASE="/media/rachele/DATA/SNVs"

# Input file (VEP-annotated, gnomAD MPC-filtered)
INPUT_RAW="${BASE}/Rubiu.v2-GRCh38.deepvariant.splitted.norm.vep.merged.vcf.gz"

# Log file — captures all terminal output
LOG="${BASE}/pipeline_$(date +%Y%m%d_%H%M%S).log"

# Redirect all stdout and stderr to both terminal and log file
exec > >(tee -a "$LOG") 2>&1
 
echo "Pipeline started: $(date)"
echo "Log file: $LOG"

# Output directories
DIR_ANNOTATED="${BASE}/01_MPC_annotated"
DIR_FILTERED="${BASE}/02_filtered"
DIR_CLINVAR="${BASE}/03_clinvar_annotated"
DIR_PATHOGENIC="${BASE}/04_pathogenic"
DIR_TSV="${BASE}/05_tsv"
DIR_COUNT="${BASE}/06_tsv"

mkdir -p "$DIR_ANNOTATED" "$DIR_FILTERED" "$DIR_CLINVAR" "$DIR_PATHOGENIC" "$DIR_TSV" "$DIR_COUNT"

# ClinVar reference file (update date as needed)
CLINVAR="${BASE}/clinvar_20250729.chr.vcf.gz"

# MPC refrence file 
MPC="${BASE}/MPC_v2.GRCh38.bcf"

# TSV conversion scripts
VCF_TO_TSV_PASTEUR="/home/rachele/SNVs/scripts_with_pasteur_data/vcf_to_tsv_PASTEUR/"
VCF_TO_TSV_CLINVAR="/home/rachele/SNVs/scripts_with_pasteur_data/vcf_to_tsv_ClinVAR"

# ----- BASIC CHECKS-----------------

# Index input file if not already indexed
if [ ! -f "${INPUT_RAW}.tbi" ]; then
  echo "Checking sort order..."
  
  SORTED_INPUT="${BASE}/Rubiu.v2-GRCh38.deepvariant.splitted.norm.vep.merged.sorted.vcf.gz"
  
  if [ ! -f "$SORTED_INPUT" ]; then
    echo "Sorting input file (this may take a while)..."
    bcftools sort \
      -Oz -o "$SORTED_INPUT" \
      --temp-dir "${BASE}/tmp" \
      "$INPUT_RAW"
  fi

  echo "Indexing sorted file..."
  tabix -p vcf "$SORTED_INPUT"
  
  # Update INPUT_RAW to point to sorted file
  INPUT_RAW="$SORTED_INPUT"
  echo "Input file updated to sorted version: $INPUT_RAW"
else
  echo "Input index already exists, skipping"
fi

# --- STEP 1 ADD BCFT TOOLS ANNOTATION -----------------------------------------
# Add the MPC v2 annotations 
echo ">>> STEP 1: MPC annotation"
MPC_ANNOTATED="${DIR_ANNOTATED}/MPC_annotated.vcf.gz"

bcftools annotate \
  -a "$MPC" \
  -c INFO \
  -Oz -o "$MPC_ANNOTATED" \
  "$INPUT_RAW"

tabix -p vcf "$MPC_ANNOTATED"


# --- STEP 2: ALLELE FREQUENCY FILTER -----------------------------------------
# Exclude variants with gnomAD genome AF > 1% (common variants)

echo ">>> STEP 2: Allele frequency filter"
AF_FILTERED="${DIR_FILTERED}/AF_filtered.vcf.gz"

bcftools view \
  --threads 22 \
  -e 'INFO/gnomad_AF > 0.01' \
  -Oz -o "$AF_FILTERED" \
  "$MPC_ANNOTATED"

tabix -p vcf "$AF_FILTERED"

echo "AF-filtered output: $AF_FILTERED"


# --- STEP 3A: MISSENSE BRANCH ------------------------------------------------
# Filter for high-impact missense variants:
#   - MPC >= 2.6 (missense pathogenicity constraint)
#   - AlphaMissense pathogenicity >= 0.565
#   - Canonical transcript only

echo ""
echo ">>> STEP 3A: Missense filter (MPC + AlphaMissense + Canonical)"

MISSENSE_MPC="${DIR_FILTERED}/missense_MPC.vcf.gz"
MISSENSE_FINAL="${DIR_FILTERED}/missense_MPC_AM_canonical.vcf.gz"

# First pass: MPC score threshold
bcftools view \
  --threads 22 \
  -i 'MPC >= 2.6' \
  -Oz -o "$MISSENSE_MPC" \
  "$AF_FILTERED"

# Second pass: AlphaMissense pathogenicity + canonical transcript
filter_vep \
  -i "$MISSENSE_MPC" \
  -o "$MISSENSE_FINAL" \
  --only_matched \
  --filter "am_pathogenicity >= 0.565 and CANONICAL is YES"

echo "Missense output: $MISSENSE_FINAL"


# --- STEP 3B: PTV BRANCH -----------------------------------------------------
# Filter for high-confidence protein-truncating variants (PTVs)
# on canonical transcripts only (LOFTEE HC)

echo ""
echo ">>> STEP 3B: PTV filter (LoF HC + Canonical)"

PTV_FINAL="${DIR_FILTERED}/PTV_HC_canonical.vcf.gz"

filter_vep \
  -i "$AF_FILTERED" \
  -o "$PTV_FINAL" \
  --only_matched \
  --filter "LoF is HC and CANONICAL is YES"

echo "PTV output: $PTV_FINAL"


# --- STEP 4: RECOMPRESS + INDEX ----------------------------------------------
# Ensure BGZF compression and tabix indexing for downstream bcftools steps

echo ""
echo ">>> STEP 4: Recompress and index filtered files"

for input in "$MISSENSE_FINAL" "$PTV_FINAL"; do
  echo "Processing: $input"

  if file "$input" | grep -q "BGZF"; then
    echo "  Already BGZF-compressed, skipping recompression"
  else
    echo "  Recompressing with bgzip..."
    mv "$input" "${input}.tmp"
    bgzip -c "${input}.tmp" > "$input"
    rm "${input}.tmp"
  fi

  echo "  Indexing with tabix..."
  tabix -p vcf "$input"
  echo "  Done: $input"
done


# --- STEP 5: CLINVAR ANNOTATION ----------------------------------------------
# Add ClinVar fields: CLNSIG, CLNDN, GENEINFO, CLNVC etc.

echo ""
echo ">>> STEP 5: ClinVar annotation"

CLINVAR_FIELDS="INFO/CLNSIG,INFO/CLNSIGCONF,INFO/CLNDN,INFO/CLNDNINCL,INFO/CLNDISDB,INFO/GENEINFO,INFO/CLNVC"

PTV_CLINVAR="${DIR_CLINVAR}/PTV_HC_canonical_ClinVar.vcf.gz"
MISSENSE_CLINVAR="${DIR_CLINVAR}/missense_MPC_AM_canonical_ClinVar.vcf.gz"

bcftools annotate \
  --threads 22 \
  -a "$CLINVAR" \
  -c "$CLINVAR_FIELDS" \
  -Oz -o "$PTV_CLINVAR" \
  "$PTV_FINAL"

bcftools annotate \
  --threads 22 \
  -a "$CLINVAR" \
  -c "$CLINVAR_FIELDS" \
  -Oz -o "$MISSENSE_CLINVAR" \
  "$MISSENSE_FINAL"

# Index annotated files
tabix -p vcf "$PTV_CLINVAR"
tabix -p vcf "$MISSENSE_CLINVAR"

echo "ClinVar-annotated PTV:      $PTV_CLINVAR"
echo "ClinVar-annotated Missense: $MISSENSE_CLINVAR"


# --- STEP 6: CLINVAR PATHOGENICITY FILTER ------------------------------------
# Lenient filter: catches all P/LP combinations including mixed annotations

echo ""
echo ">>> STEP 6: ClinVar pathogenicity filter (lenient)"

lenient_filter='INFO/CLNSIG~"^Pathogenic$" ||
  INFO/CLNSIG~"^Likely_pathogenic$" ||
  INFO/CLNSIG~"Pathogenic/Likely_pathogenic" ||
  INFO/CLNSIG~"Likely_pathogenic/Pathogenic" ||
  INFO/CLNSIG~"Pathogenic\\|Likely_pathogenic" ||
  INFO/CLNSIG~"Pathogenic,_low_penetrance" ||
  INFO/CLNSIG~"Likely_pathogenic,_low_penetrance" ||
  INFO/CLNSIG~"Likely_pathogenic|" ||
  INFO/CLNSIG~"|Likely_pathogenic" ||
  INFO/CLNSIG~"Pathogenic|" ||
  INFO/CLNSIG~"|Pathogenic" ||
  INFO/CLNSIG~"Affects" ||
  INFO/CLNSIG~"Pathogenic|Affects"'

for category in "PTV_HC_canonical" "missense_MPC_AM_canonical"; do
  input="${DIR_CLINVAR}/${category}_ClinVar.vcf.gz"
  output="${DIR_PATHOGENIC}/${category}_ClinVar_pathogenic.vcf.gz"

  echo "Filtering: $input"
  bcftools view --threads 22 -i "$lenient_filter" -Oz -o "$output" "$input"
  tabix -p vcf "$output"
  echo "Output: $output"
  echo "Variant count: $(bcftools view -H "$output" | wc -l)"
  echo ""
done


# --- STEP 7: EXPORT TO TSV ---------------------------------------------------

echo ""
echo ">>> STEP 7: Export to TSV"

cd "$VCF_TO_TSV_PASTEUR"
# Non-ClinVar filtered files
bash ./transform_VCF.sh -i "$PTV_FINAL"
bash ./transform_VCF.sh -i "$MISSENSE_FINAL"

cd "$VCF_TO_TSV_CLINVAR"
# ClinVar pathogenic files
bash ./transform_VCF.sh -i "${DIR_PATHOGENIC}/PTV_HC_canonical_ClinVar_pathogenic.vcf.gz"
bash ./transform_VCF.sh -i "${DIR_PATHOGENIC}/missense_MPC_AM_canonical_ClinVar_pathogenic.vcf.gz"


#not sure if we should as welltrasnform the oncenot filtered byclinvar but with clinvar annotation


# --- STEP 8: QC — VARIANT COUNTS ---------------------------------------------
# Summary of variant counts across all output files

echo ""
echo ">>> STEP 8: QC — Variant counts"
echo ""

QC_COUNTS="${DIR_COUNT}/QC_variant_counts.tsv"
echo -e "Step\tFile\tVariant_Count" > "$QC_COUNTS"
 
declare -A step_labels=(
  ["$AF_FILTERED"]="1_AF_filter"
  ["$MISSENSE_MPC"]="2a_Missense_MPC"
  ["$MISSENSE_FINAL"]="2a_Missense_MPC_AM_canonical"
  ["$PTV_FINAL"]="2b_PTV_HC_canonical"
  ["${DIR_PATHOGENIC}/PTV_HC_canonical_ClinVar_pathogenic.vcf.gz"]="5_PTV_pathogenic"
  ["${DIR_PATHOGENIC}/missense_MPC_AM_canonical_ClinVar_pathogenic.vcf.gz"]="5_Missense_pathogenic"
)
 
echo ""
printf "%-35s %-60s %s\n" "Step" "File" "Variant_Count"
printf "%-35s %-60s %s\n" "----" "----" "-------------"
 
for file in \
  "$AF_FILTERED" \
  "$MISSENSE_MPC" \
  "$MISSENSE_FINAL" \
  "$PTV_FINAL" \
  "${DIR_PATHOGENIC}/PTV_HC_canonical_ClinVar_pathogenic.vcf.gz" \
  "${DIR_PATHOGENIC}/missense_MPC_AM_canonical_ClinVar_pathogenic.vcf.gz"; do
    count=$(bcftools view -H "$file" | wc -l)
    step="${step_labels[$file]}"
    fname="$(basename "$file")"
    printf "%-35s %-60s %s\n" "$step" "$fname" "$count"
    echo -e "${step}\t${fname}\t${count}" >> "$QC_COUNTS"
done
 
echo ""
echo "Variant counts saved to: $QC_COUNTS"

# Also dump CLNSIG annotation breakdown per ClinVar-annotated file
CLNSIG_REPORT="${DIR_COUNT}/QC_ClinVar_annotation_counts.tsv"
echo -e "File\tCLNSIG_Annotation\tCount" > "$CLNSIG_REPORT"

for file in "$PTV_CLINVAR" "$MISSENSE_CLINVAR"; do
  category=$(basename "$file" .vcf.gz)
  bcftools query -f '%INFO/CLNSIG\n' "$file" |
    sort | uniq -c |
    awk -v f="$category" 'BEGIN{OFS="\t"} {print f, $2, $1}' >> "$CLNSIG_REPORT"
done

echo ""
echo "ClinVar annotation breakdown saved to: $CLNSIG_REPORT"
echo ""
echo ">>> Pipeline complete: $(date)"
echo ">>> Log saved to: $LOG"