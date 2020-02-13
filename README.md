#CADD-SV is a scoring tool for structural variants that uses SVs that survived millions of years of purifying selection to train a linear model using multiple features (see below) to estimate the effect of novel SVs.

##philip.kleinert@bihealth.de

#annotations

CADD 
https://krishna.gs.washington.edu/download/CADD/bigWig/
mean & max of CADD-Phred Score

CTCF 
http://genome.cshlp.org/content/suppl/2012/08/28/22.9.1680.DC1/Table_S2_Location_of_ChIP-seq_binding_positions_in_19_cell_lines.txt
Chip Seq from 19 celllines from wang et al

Encode-HIC data
https://www.encodeproject.org/search/?type=Experiment&assay_term_name=HiC&replicates.library.biosample.donor.organism.scientific_name=Homo%20sapiens&status=released
HIC processed data for nested tads and TAD domains for Caki2 G401 LNCaP NCI-H460 Panc1 RPMI-7951 SJCRH30 SK-MEL-5 SK-N-MC T47D

enhancer-promoter-links
from FOCS https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5930446/

ensembl-gff3 genebuild
ftp://ftp.ensembl.org/pub/release-96/gff3/homo_sapiens/Homo_sapiens.GRCh38.96.chr.gff3.gz

FIREs frequently interacting regulatory elements
https://www.sciencedirect.com/science/article/pii/S2211124716314814
for "gm12878","msc","mes","imr90","h1"-cell lines

GC content in 5bp resolution
http://hgdownload.cse.ucsc.edu/gbdb/hg38/bbi/gc5BaseBw/gc5Base.bw

genomegitar
https://www.genomegitar.org/processed-data.html
Directionality index for HIC data from various datasets

hic in hg18 (remove)
chromosome.sdsc.edu/mouse/hi-c/hESC.domain.tar.gz
human ESC domains

PhastCons
http://hgdownload.cse.ucsc.edu/goldenpath/hg38/

PLI score
ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/manuscript_data/forweb_cleaned_exac_r03_march16_z_data_pLI.txt.gz

Syntenic regions from
http://webclu.bio.wzw.tum.de/cgi-bin/syntenymapper/get-species-list.py

chromhmm states of encode celllines(Philipp) 
chromHMM_GRCh38.bg.gz

gerp score 
http://mendel.stanford.edu/SidowLab/downloads/gerp/

conserved coding regions (CCR)
https://www.nature.com/articles/s41588-018-0294-6?WT.feed_name=subjects_population-genetics



