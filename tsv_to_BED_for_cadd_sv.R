# Use to convert annotation containing TSV file to bed files as input for CADD.
# submit as Rscript tsv_to_BED_for_caddsv.R <tsv_file> <samplename_NOCAPS>
library(dplyr)
library(yaml)
library(stringi)
warning("This conversion does not take into account the start 0 for BED and start 1 for VCF files")

args <- commandArgs(trailingOnly = TRUE)
file_tsv <- args[1]
samplename <- args[2]
# file_tsv <- "/Users/sreenivasan/Desktop/CADD-SV/R22-016.GRCh38.pbsv.vcf_annotated.xlsx" ##example

columns_oI <- c("chr", "POS.start", "POS-start", "POS.END", "POS-END", "SVTYPE")
variants <- read.table(file_tsv, sep = "\t", header = TRUE) %>%
    select(any_of(columns_oI))
if (dim(variants)[2] != 4) {stop("Check the column names of the TSV file")}

#convert inversions (INV) and translocations (BND) to deletions for CADD-SV
variants <- variants %>%
    filter(SVTYPE != "BND") %>%
    mutate(SVTYPE = case_when(SVTYPE == "INV" ~ "DEL",
                                TRUE ~ SVTYPE))
if(!all(unique(variants$SVTYPE) %in% c("DEL", "DUP", "INS"))) stop("Check for more SVTYPES")

# remove random chromosomes
variants  <- variants %>%
    filter(!grepl(chr, pattern="_random")) %>%
    filter(!grepl(chr, pattern="_alt")) %>%
    filter(!grepl(chr, pattern="chrUn")) %>%
    filter(!grepl(chr, pattern="chrM"))

# remove duplicates
variants <- variants %>%
    distinct(across(.cols = everything())) %>%
    arrange(across(1:3))

# Write bed files with no more than 10k variants.
outdir <- paste0(dirname(file_tsv), "/caddsv_analysis/", samplename)
dir.create(path = paste0(outdir, "/input"), recursive = TRUE)
filename <- paste0(outdir, "/input/id_", samplename, ".bed")

# Write the complete  bedfile
message("Writing the complete file")
write.table(variants, file = filename, col.names = FALSE, row.names = FALSE, , sep = "\t", quote = FALSE)

# For config.yml file
cadd_yaml_list <- list(dataset = c(samplename))

split <- TRUE
if(split) {
    variants_split <- split(variants, f = seq_along(along.with=variants[, 1]) %/% 9999)
    for (chunk in seq_along(variants_split)){
        message("Writing split bed files, Chunk #", chunk)
        filename <- paste0(outdir, "/input/id_", samplename, chunk, ".bed")
        write.table(variants_split[[chunk]], file = filename, col.names = FALSE, row.names = FALSE, , sep = "\t", quote = FALSE)
        cadd_yaml_list[["dataset"]] <- c(cadd_yaml_list[["dataset"]], paste0(samplename, chunk))
    }
}

# Now create a config.yml file for CADD-SV
filename <- paste0(outdir, "/config.yml")
write_yaml(x = cadd_yaml_list, file = filename)

# Now use the web-service or installation in OMICS to run CADD-SV calls