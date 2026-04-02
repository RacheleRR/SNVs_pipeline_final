# CLINVAR renmae chr to match my files 

echo -e "1\tchr1\n2\tchr2\n3\tchr3\n4\tchr4\n5\tchr5\n6\tchr6\n7\tchr7\n8\tchr8\n9\tchr9\n10\tchr10\n11\tchr11\n12\tchr12\n13\tchr13\n14\tchr14\n15\tchr15\n16\tchr16\n17\tchr17\n18\tchr18\n19\tchr19\n20\tchr20\n21\tchr21\n22\tchr22\nX\tchrX\nY\tchrY\nMT\tchrM" > rename_chr.txt
bcftools annotate --rename-chrs rename_chr.txt clinvar_20250729.vcf.gz | bgzip > clinvar_20250729.chr.vcf.gz
tabix clinvar_20250729.chr.vcf.gz


#GENE SETS
# Load my data from CSV files
SCHEMA <- read_csv("/home/rachele/Documents/old/geneset_confrontation/SCHEMA.csv")
GWAS_120 <- read.delim("~/GWAS_120.csv")       
BipEx_Bipolar <- read_csv("/home/rachele/Documents/old/geneset_confrontation/Bipolar_Disorder.csv")


# Clean GENESETS
# Convert gene names to a standard format for easier comparison
# BipEx_Bipolar
convert_col <- gconvert(query = BipEx_Bipolar$Gene, organism = "hsapiens", target = "ENSG", mthreshold = Inf, filter_na = TRUE)
BipEx_Bipolar <- merge(BipEx_Bipolar, convert_col[, c("target", "name")], by.x = "Gene", by.y = "target", all.x = TRUE)
BipEx_Bipolar <- BipEx_Bipolar %>% select(Gene, name, everything())

# SCHEMA
convert_col <- gconvert(query = SCHEMA$Gene, organism = "hsapiens", target = "ENSG", mthreshold = Inf, filter_na = TRUE)
SCHEMA <- merge(SCHEMA, convert_col[, c("target", "name")], by.x = "Gene", by.y = "target", all.x = TRUE)
SCHEMA <- SCHEMA %>% select(Gene, name, everything())

# General clean-up of data
# Remove unnecessary columns, rename columns, and clean up gene names
BipEx_Bipolar <- BipEx_Bipolar %>% subset(select = -Gene) %>% rename('Gene' = name) %>% mutate(Gene = gsub(' ', '', Gene)) %>% setNames(make.names(colnames(.))) 
BipEx_Bipolar_p_val_PTV <- BipEx_Bipolar[!is.na(BipEx_Bipolar$PTV.Fisher.p.val) & BipEx_Bipolar$PTV.Fisher.p.val <= 0.05,]
BipEx_Bipolar_p_val_Missense <- BipEx_Bipolar[!is.na(BipEx_Bipolar$Damaging.Missense.Fisher.p.val) & BipEx_Bipolar$Damaging.Missense.Fisher.p.val <= 0.05,]

BipEx_Bipolar_p_val_Missense <- BipEx_Bipolar_p_val_Missense %>%filter(!is.na(Damaging.Missense.Fisher.odds.ratio), Damaging.Missense.Fisher.odds.ratio > 1)
BipEx_Bipolar_p_val_PTV <- BipEx_Bipolar_p_val_PTV %>%filter(!is.na(PTV.Fisher.odds.ratio), PTV.Fisher.odds.ratio > 1)
SCHEMA <- SCHEMA %>% subset(select = -Gene) %>% rename('Gene' = name) %>% mutate(Gene = gsub(' ', '', Gene)) %>% setNames(make.names(colnames(.)))  
SCHEMA_pVal <- SCHEMA[SCHEMA$P.meta <= 0.05, ]
SCHEMA_qVal <- SCHEMA[SCHEMA$Q.meta <= 0.05, ]
SCHEMA_pVal <- SCHEMA_pVal[ (!is.na(SCHEMA_pVal$OR..Class.I.) & SCHEMA_pVal$OR..Class.I. > 1) | (!is.na(SCHEMA_pVal$OR..Class.II.) & SCHEMA_pVal$OR..Class.II. > 1),]
SCHEMA_qVal <- SCHEMA_qVal[ (!is.na(SCHEMA_qVal$OR..Class.I.) & SCHEMA_qVal$OR..Class.I. > 1) | (!is.na(SCHEMA_qVal$OR..Class.II.) & SCHEMA_qVal$OR..Class.II. >1),]
GWAS_120 <- GWAS_120 %>% rename('Gene' = GENE_name)
BipEx_Bipolar_combined <- bind_rows(
  BipEx_Bipolar_p_val_PTV,
  BipEx_Bipolar_p_val_Missense
) %>%
  distinct()

brain_gene_consensus_filtered_consensus_no_pitular <- brain_gene_consensus_filtered_consensus_no_pitular %>%rename('GENE_name' = Gene)
brain_gene_consensus_filtered_consensus_no_pitular <- brain_gene_consensus_filtered_consensus_no_pitular %>%rename('Gene' = Gene.name)

brain_gene_consensus_ntm_consensus_no_pitular <- brain_gene_consensus_ntm_consensus_no_pitular %>%rename('GENE_name' = Gene)
brain_gene_consensus_ntm_consensus_no_pitular <- brain_gene_consensus_ntm_consensus_no_pitular %>%rename('Gene' = Gene.name)

write_csv(SCHEMA_pVal, "SCHEMA_pval.csv")
write_csv(SCHEMA_qVal, "SCHEMA_qval.csv")
write_csv(GWAS_120, "GWAS_120.csv")
# Combined dataset
write_csv(BipEx_Bipolar_combined, "BipEx_Bipolar.csv")
write_csv(brain_gene_consensus_ntm_consensus_no_pitular, "brain_gene_consensus_ntm_consensus_no_pitular.csv")
write_csv(brain_gene_consensus_filtered_consensus_no_pitular,"brain_gene_consensus_filtered_consensus_no_pitular")