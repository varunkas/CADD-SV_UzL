# CADD-SV 
## Notes for using in UzL OMICS cluster for PacBIO and Nanopore data

### Conversion of VCF.GZ to BED.

Convert the pbsv.vcf.gz to bedfile using the following script. Do not use the annotated vcf because some variants go missing during the annotation.

        Only selects DEL, DUP, and INS
        /Users/sreenivasan/Documents/Works/scripts/CADD-SV_UzL/vcfgz_to_bed.sh <vcfgz_file> <name>
    
The ```<name>``` will be used as the 5th column in the output bedfile. This way multiple bedfiles can be combined for a single CADD-SV run by using:

        cat xx >> id_combined.unsorted.bed

        # Sort using R, because bedtools sort also sorts based on the name column. Annoying!

        R
        bedfile <- read.table("id_combined.unsorted.bed", sep="\t", header=FALSE)
        bedfile <- bedfile %>% arrange(V1, V2, V3)
        write.table(bedfile, "id_combined.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

Submit this for CADD-SV in the omics cluster using the command: 

    ./cadd_sv_main.sh
    # Note the config.yml might needs editing. 

### Splitting of CADD-SV output
The output ```combined_score.bed``` now contains the scores for all samples. This needs to be split using:

    # Get unique names (col #5) from the first 1000 lines
    names=$(tail -n +2 combined_score.bed | head -n 1000 | cut -f5 | sort | uniq)
    
    names=( "R21-007" "R21-015" "R22-005" "R22-016" )

    # Split them
    for name in $names
    do
        echo $name
        awk -v name=$name 'BEGIN {OFS = "\t"} {if($5 == name) print chr$1, $2, $3, $4, $5, $6}' combined_score.bed > ${name}_cadd-sv_score.bed
    done

**Check the pacbio_custom_analysis repo for information on how to combine the CADD-SV scores with VEP annotations!!**

## CADD-SV – a framework to score the effect of structural variants 

Here, we describe CADD-SV, a method to retrieve a wide set of annotations in the range and vicinity of a SV. Our tool computes summary statistics and uses a trained machine learning model to differentiate deleterious from neutral variants. In training, we use human and chimpanzee derived alleles as proxy-neutral and contrast them with matched simulated variants as proxy-pathogenic. This approach has proven powerful in the interpretation of SNVs (CADD, https://cadd.gs.washington.edu). We show that CADD-SV scores correlate with known pathogenic variants in individual genomes and allelic diversity.


## Pre-requirements

### Conda

The pipeline depends on [Snakemake](https://snakemake.readthedocs.io/en/stable/), a workflow management system that wraps up all scripts and runs them highly automated, in various environments (workstations, clusters, grid, or cloud). Further, we use Conda as software/dependency management tool. Conda can install snakemake and all neccessary software with its dependencies automatically. Conda installation guidelines can be found here:

https://conda.io/projects/conda/en/latest/user-guide/install/index.html

### Snakemake

After installing Conda, you install Snakemake using Conda and the `environment.yaml` provided in this repository. For this purpose, please clone or download and uncompress the repository first. Then change into the root folder of the local repository. 

```bash
git clone https://github.com/kircherlab/CADD-SV
cd CADD-SV
```

We will now initiate the Conda environment, which we will need for getting the Snakemake workflow invoked. Using this environment (`run.caddsv`) snakemake will be installed

```bash
conda env create -n run.caddsv --file environment.yaml
```

The second conda environment (`envs/SV.yml`), containing all packages and tools to run CADD-SV, will be installed automatically during the first run. This can take some time.

### Annotations

CADD-SV depends on various annotations to provide the model with its necessary input features. CADD-SV automatically retrieves and transforms these annotations (see Snakefile) and combines them in bed-format at `/desired-sv-set/matrix.bed`

Annotations can be downloaded and expanded individually. However, to run CADD-SV with the pre-trained model and to minimize runtime and memory failures use the annotation sets as stored at https://kircherlab.bihealth.org/download/CADD-SV/

```bash
wget https://kircherlab.bihealth.org/download/CADD-SV/v1.0/dependencies.tar.gz
tar -xf dependencies.tar.gz
```

## Config

Almost ready to go. After you prepared the files above, you may need to adjust the name of your dataset in the `config.yml`. 

## List of required input files

- Models and scripts as cloned from this GIT repository
- Annotations in the `annotations/` folder
- CADD-SV scores SV in a coordinate sorted BED format on the GRCh38 genome build. The type of SV needs to be included for each variant in the 4th column. We recommend to split files containing more than 10,000 SVs into smaller files. An example input file can be found in `input/`. The file needs to have the suffix `id_`. If you plan to process variants from another genome build or SVs in VCF format, see below.

## Running the pipeline

Ready to go! If you run the pipeline on a cluster see the `cluster.json` for an estimate of minimum resource requirements for the individual jobs. Note that this depends on your dataset size so you may have to adjust this.

To start the pipeline:

```bash
conda activate run.caddsv
# dry run to see if everything works
snakemake  --use-conda --configfile config.yml -j 4 -n
# run the pipeline
snakemake  --use-conda --configfile config.yml -j 4
```

## Output files

The pipeline outputs your SV set containing all annotations in BED format in a folder named `output` containing the CADD-SV and two raw scores in rows 6-8.
Further information about individual annotations are kept in a subfolder named after your input dataset.


# Further Information

## Annotations

CADD-SV integrates different annotations, here some links to its annotation sources. A complete list can be found as Suppl. Table 1 of the manuscript/pre-print. 

##### Integrated Scores
CADD (https://krishna.gs.washington.edu/download/CADD/bigWig/) \
LINSIGHT (http://compgen.cshl.edu/LINSIGHT/LINSIGHT.bw) 

##### Species conservation and constraint metrics
PhastCons (http://hgdownload.cse.ucsc.edu/goldenpath/hg38/) \
Syntenic regions (http://webclu.bio.wzw.tum.de/cgi-bin/syntenymapper/get-species-list.py) \
GERP score (http://mendel.stanford.edu/SidowLab/downloads/gerp/) 

##### Population and disease constraint metrics
pLI score (ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/manuscript_data/forweb_cleaned_exac_r03_march16_z_data_pLI.txt.gz) \
Conserved coding regions (CCR) (https://www.nature.com/articles/s41588-018-0294-6?WT.feed_name=subjects_population-genetics) \
DDD Happloinsufficiency (https://decipher.sanger.ac.uk/files/downloads/HI_Predictions_Version3.bed.gz) 

##### Epigenetic and regulatory activity
Encode Features such as Histon Modifications and DNase and RNase-seq (https://www.encodeproject.org/help/batch-download/) \
GC content (http://hgdownload.cse.ucsc.edu/gbdb/hg38/bbi/gc5BaseBw/gc5Base.bw) \
ChromHMM states of ENCODE cell lines (http://compbio.mit.edu/ChromHMM/) 

##### 3D genome organization
CTCF (http://genome.cshlp.org/content/suppl/2012/08/28/22.9.1680.DC1/Table_S2_Location_of_ChIP-seq_binding_positions_in_19_cell_lines.txt) \
Encode-HiC data (https://www.encodeproject.org/search/?type=Experiment&assay_term_name=HiC&replicates.library.biosample.donor.organism.scientific_name=Homo%20sapiens&status=released) \
Enhancer-promoter-links from FOCS (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5930446/) \
Frequently Interacting Regulatory Elements (FIREs) (https://www.sciencedirect.com/science/article/pii/S2211124716314814) \
Directionality index for HiC data from various datasets (https://www.genomegitar.org/processed-data.html) \
DeepC saliencies score (http://userweb.molbiol.ox.ac.uk/public/rschwess/container_for_ucsc/data/deepC/saliency_scores/saliencies_merged_gm12878_5kb.bw) 

##### Gene and element annotations
Ensembl-gff3 genebuild 96 (ftp://ftp.ensembl.org/pub/release-96/gff3/homo_sapiens/Homo_sapiens.GRCh38.96.chr.gff3.gz) \
Fantom5 enhancers (https://zenodo.org/record/556775#.Xkz3G0oo-70) 


## Converting VCF and other genome builds

  If you want to score SVs in a VCF format or your SVs are not in GRCh38 genomebuild coordinates:
  We provide an environment to handle this. It uses the SURVIVOR tools (https://github.com/fritzsedlazeck/SURVIVOR).
  
  ```bash
  conda env create -n prepBED --file envs/prepBED.yml
  ```
  
  To convert your VCF into BED format run:
  ```
  conda activate prepBED
  SURVIVOR vcftobed input.vcf 0 -1 output.bed
  cut -f1,2,6,11 output.bed > beds/set_id.bed
  
  ```
  
  To lift hg19 coordinates to GRCh38 apply the following steps:
  
  ```
  conda activate prepBED
  liftOver beds/setname_hg19_id.bed dependencies/hg19ToHg38.over.chain.gz beds/setname_id.bed beds/setname_unlifted.bed
  ```
   
