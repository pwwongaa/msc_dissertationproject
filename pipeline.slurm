#!/bin/bash
#SBATCH --job-name=bioinfo_pipeline
#SBATCH --account=BISC029044
#SBATCH --partition=teach_cpu
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --output=pipeline_output.log
#SBATCH --error=pipeline_error.log
#SBATCH --time=7-00:00:00

# Load necessary modules
module load nanoplot chopper flye quast busco braker muscle trimal iqtree antismash dbcan RepeatModeler RepeatMasker orthofinder
module add lang/python/anaconda/3.8-2021-TW

# Activate necessary conda environments
source activate busco
source activate orthof
source activate phylog

# Set up directories
BASE_DIR="/user/work/be22534/dissertation_project"
RAW_DIR="${BASE_DIR}/data/raw"
NANOPLOT_DIR="${BASE_DIR}/nanoplot"
CHOPPER_DIR="${BASE_DIR}/chopper"
FLYE_DIR="${BASE_DIR}/flye"
QUAST_BUSCO_DIR="${BASE_DIR}/quast_busco"
BRAKER2_DIR="${BASE_DIR}/braker2"
PHYLOGENOMICS_DIR="${BASE_DIR}/phylogenomics"
FUNCTIONAL_ANNOTATIONS_DIR="${BASE_DIR}/functional_annotations"
REPEAT_DIR="${BASE_DIR}/repeat"
DOWNLOAD_DIR="${BASE_DIR}/downloaded_fasta"

mkdir -p $NANOPLOT_DIR $CHOPPER_DIR $FLYE_DIR $QUAST_BUSCO_DIR $BRAKER2_DIR $PHYLOGENOMICS_DIR/muscle5 $PHYLOGENOMICS_DIR/trimal $PHYLOGENOMICS_DIR/iqtree2 $FUNCTIONAL_ANNOTATIONS_DIR $REPEAT_DIR $DOWNLOAD_DIR

# Step 1: Nanoplot for QC
for FILE in T17_02.fastq T17_15.fastq T17_42.fastq; do
    OUTPUT_DIR="${NANOPLOT_DIR}/nanoplot_output_${FILE%.fastq}"
    NanoPlot --fastq ${RAW_DIR}/${FILE} --plots kde --raw --tsv_stats --outdir $OUTPUT_DIR
done

# Step 2: Chopper
# Preprocess raw data with Chopper (assuming correct Chopper command)
chopper -i ${RAW_DIR} -o ${CHOPPER_DIR}

# Step 3: Flye assembly
for FILE in T17_02.fastq T17_15.fastq T17_42.fastq; do
    OUTPUT_DIR="${FLYE_DIR}/flye_output_${FILE%.fastq}"
    flye --nano-raw ${RAW_DIR}/${FILE} --out-dir $OUTPUT_DIR --threads 12
done

# Step 4: RepeatModeler and RepeatMasker
for FILE in T17_02 T17_15 T17_42; do
    FLYE_OUTPUT="${FLYE_DIR}/flye_output_${FILE}/assembly.fasta"
    REPEAT_OUTPUT="${REPEAT_DIR}/repeat_output_${FILE}"
    
    mkdir -p $REPEAT_OUTPUT
    cd $REPEAT_OUTPUT
    BuildDatabase -name ${FILE} $FLYE_OUTPUT
    RepeatModeler -database ${FILE} -pa 8
    RepeatMasker -gff -pa 8 -lib ${REPEAT_OUTPUT}/consensi.fa.classified $FLYE_OUTPUT
done

# Step 5: BRAKER2 annotations
for FILE in T17_02 T17_15 T17_42; do
    BRAKER_OUTPUT="${BRAKER2_DIR}/braker2_output_${FILE}"
    GENOME_FILE="${FLYE_DIR}/flye_output_${FILE}/assembly.fasta.masked"
    braker.pl --genome=$GENOME_FILE --prot_seq=uniprot_ascomycota-filtered-reviewed.fasta --cores 8 --fungus --AUGUSTUS_ab_initio --softmasking --gff3 --workingdir=$BRAKER_OUTPUT
done

# Step 6: BUSCO after annotations
for FILE in T17_02 T17_15 T17_42; do
    ANNOTATED_GENOME="${BRAKER2_DIR}/braker2_output_${FILE}/augustus.hints.gff3"
    BUSCO_OUTPUT="${QUAST_BUSCO_DIR}/busco_output_${FILE}_annotations"
    
    busco -i $ANNOTATED_GENOME -o $BUSCO_OUTPUT -l sordariomycetes_odb10 -m genome --cpu 4 -f
done

# Step 7: OrthoFinder before phylogenomics
ORTHOFINDER_INPUT="${PHYLOGENOMICS_DIR}/orthofinder_input"
ORTHOFINDER_OUTPUT="${PHYLOGENOMICS_DIR}/orthofinder_output"

mkdir -p $ORTHOFINDER_INPUT $ORTHOFINDER_OUTPUT

# Assuming you have multiple proteomes to analyze
# Copy all proteomes to the ORTHOFINDER_INPUT directory
cp /path/to/proteomes/*.faa $ORTHOFINDER_INPUT

orthofinder -f $ORTHOFINDER_INPUT -t 8 -o $ORTHOFINDER_OUTPUT

# Step 8: Phylogenomics analysis
for FILE in ${ORTHOFINDER_OUTPUT}/*.fa; do
    FILE_NO_EXT=$(basename ${FILE} .fa)
    muscle -align $FILE -output ${PHYLOGENOMICS_DIR}/muscle5/${FILE_NO_EXT}_aln.fa
    trimal -in ${PHYLOGENOMICS_DIR}/muscle5/${FILE_NO_EXT}_aln.fa -out ${PHYLOGENOMICS_DIR}/trimal/${FILE_NO_EXT}_trim.fa -automated1
done

cd ${PHYLOGENOMICS_DIR}/iqtree2
iqtree2 -s concat.fa -m MFP -B 1000 -T AUTO

# Step 9: Functional annotations
for FILE in ${BRAKER2_DIR}/*; do
    ANNOTATED_GENOME="${FILE}/augustus.hints.gff3"
    antismash --input-dir $ANNOTATED_GENOME --output-dir ${FUNCTIONAL_ANNOTATIONS_DIR}/$(basename $FILE)
    run_dbcan --input-dir $ANNOTATED_GENOME --output-dir ${FUNCTIONAL_ANNOTATIONS_DIR}/$(basename $FILE)
done

# Separate Processing for Downloaded fasta.gz Files
# Ensure each fasta.gz file is processed through QUAST, BUSCO, and annotations separately
for FILE in ${DOWNLOAD_DIR}/*.fasta.gz; do
    BASENAME=$(basename ${FILE} .fasta.gz)
    UNZIPPED_FILE="${DOWNLOAD_DIR}/${BASENAME}.fasta"
    
    # Unzip the fasta.gz file
    gunzip -c ${FILE} > ${UNZIPPED_FILE}
    
    # Run QUAST
    QUAST_OUTPUT="${QUAST_BUSCO_DIR}/quast_output_${BASENAME}"
    quast.py ${UNZIPPED_FILE} -o ${QUAST_OUTPUT}
    
    # Run RepeatModeler and RepeatMasker
    REPEAT_OUTPUT="${REPEAT_DIR}/repeat_output_${BASENAME}"
    mkdir -p $REPEAT_OUTPUT
    cd $REPEAT_OUTPUT
    BuildDatabase -name ${BASENAME} ${UNZIPPED_FILE}
    RepeatModeler -database ${BASENAME} -pa 8
    RepeatMasker -gff -pa 8 -lib ${REPEAT_OUTPUT}/consensi.fa.classified ${UNZIPPED_FILE}
    
    # Run BRAKER2 annotations
    BRAKER_OUTPUT="${BRAKER2_DIR}/braker2_output_${BASENAME}"
    GENOME_FILE="${REPEAT_OUTPUT}/masked/${BASENAME}.masked"
    braker.pl --genome=$GENOME_FILE --prot_seq=uniprot_ascomycota-filtered-reviewed.fasta --cores 8 --fungus --AUGUSTUS_ab_initio --softmasking --gff3 --workingdir=$BRAKER_OUTPUT
    
    # Run BUSCO after annotations
    ANNOTATED_GENOME="${BRAKER_OUTPUT}/augustus.hints.gff3"
    BUSCO_OUTPUT="${QUAST_BUSCO_DIR}/busco_output_${BASENAME}_annotations"
    busco -i $ANNOTATED_GENOME -o $BUSCO_OUTPUT -l sordariomycetes_odb10 -m genome --cpu 4 -f
done

echo "Pipeline completed successfully."
