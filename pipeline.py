import os
import subprocess

#?
# checked with ori slurm 

# i. set up required functions
def run_command(command):
    """To run a CML shell command with error handling"""
    print(f"Running command: {command}")
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Command '{command}' failed with error: {e}")
        exit(1)

# ii. set up directories for files management
BASE_DIR = "/user/work/dissertation_project"
RAW_DIR = os.path.join(BASE_DIR, "data/raw")
NANOPLOT_DIR = os.path.join(BASE_DIR, "nanoplot")
CHOPPER_DIR = os.path.join(BASE_DIR, "chopper")
FLYE_DIR = os.path.join(BASE_DIR, "flye")
QUAST_BUSCO_DIR = os.path.join(BASE_DIR, "quast_busco")
BRAKER2_DIR = os.path.join(BASE_DIR, "braker2")
PHYLOGENOMICS_DIR = os.path.join(BASE_DIR, "phylogenomics")
FUNCTIONAL_ANNOTATIONS_DIR = os.path.join(BASE_DIR, "functional_annotations")
REPEAT_DIR = os.path.join(BASE_DIR, "repeat")
DOWNLOAD_DIR = os.path.join(BASE_DIR, "downloaded_fasta")

directories = [
    NANOPLOT_DIR, CHOPPER_DIR, FLYE_DIR, QUAST_BUSCO_DIR, BRAKER2_DIR,
    os.path.join(PHYLOGENOMICS_DIR, "muscle5"), os.path.join(PHYLOGENOMICS_DIR, "trimal"),
    os.path.join(PHYLOGENOMICS_DIR, "iqtree2"), FUNCTIONAL_ANNOTATIONS_DIR, REPEAT_DIR, DOWNLOAD_DIR
]

for directory in directories:
    os.makedirs(directory, exist_ok=True)

# iii. activate necessary conda environments for each part if needed (?)
# conda_envs = ["busco", "orthof", "phylog"]
# for env in conda_envs:
#     run_command(f"source activate {env}")

#? Step 1: Nanoplot for QC
for file in ["TA001.fastq", "TA002.fastq", "TA003.fastq"]:
    output_dir = os.path.join(NANOPLOT_DIR, f"nanoplot_output_{os.path.splitext(file)[0]}")
    #NanoPlot - --raw, output the raw data worksheet
    command = f"NanoPlot --fastq {os.path.join(RAW_DIR, file)} --plots kde --raw --tsv_stats --outdir {output_dir}"
    run_command(command)

#? Step 2: Chopper for quality filtering
for file in ["T17_02.fastq", "T17_15.fastq", "T17_42.fastq"]:
    chopper_command = f"chopper -i {os.path.join(RAW_DIR, file)} -o {CHOPPER_DIR}"
    run_command(chopper_command)


#? Step 3: Flye - denovo whole genome assembly
for file in ["TA001.fastq", "TA002.fastq", "TA003.fastq"]:

    # Step 3a: Flye - whole genome assemlby
    output_dir = os.path.join(FLYE_DIR, f"flye_output_{os.path.splitext(file)[0]}")
    command = f"flye --nano-raw {os.path.join(RAW_DIR, file)} --out-dir {output_dir} --threads 4"
    run_command(command)

    # Step 3b: QUAST for assembly quality assessment
    quast_output_dir = os.path.join(QUAST_BUSCO_DIR, f"quast_output_{os.path.splitext(file)[0]}")
    quast_command = f"quast.py {os.path.join(output_dir, 'assembly.fasta')} -o {quast_output_dir}"
    run_command(quast_command)
    
    # Step 3c: BUSCO for assembly completeness
    busco_output_dir = os.path.join(QUAST_BUSCO_DIR, f"busco_output_{os.path.splitext(file)[0]}_assembly")
    busco_command = f"busco -i {os.path.join(output_dir, 'assembly.fasta')} -o {busco_output_dir} -l sordariomycetes_odb10 -m genome --cpu 4"
    run_command(busco_command)

#? Step 4: RepeatModeler and RepeatMasker - Repeat database building and the repeat softmasking
for file in ["TA001", "TA002", "TA003"]:
    flye_output = os.path.join(FLYE_DIR, f"flye_output_{file}/assembly.fasta")
    repeat_output = os.path.join(REPEAT_DIR, f"repeat_output_{file}")

    os.makedirs(repeat_output, exist_ok=True)
    os.chdir(repeat_output)
    
    #RepeatModeler
    build_db_command = f"BuildDatabase -name {file} -engine ncbi {flye_output}"
    run_command(build_db_command)
    repeat_modeler_command = f"RepeatModeler -database {file} -engine ncbi -thread 8"
    run_command(repeat_modeler_command)

    #RepeatMasker
    repeat_masker_command = f"RepeatMasker -pa 12 -s -gff -nolow -xsmall -lib {os.path.join(repeat_output, 'consensi.fa.classified')} {flye_output}"
    run_command(repeat_masker_command)

#? Step 5: BRAKER2 - gene annotations
for file in ["TA001", "TA002", "TA003"]:
    braker_output = os.path.join(BRAKER2_DIR, f"braker2_output_{file}")
    genome_file = os.path.join(FLYE_DIR, f"flye_output_{file}/assembly.fasta.masked")
    
    braker_command = (
        f"braker.pl --genome={genome_file} \n"
        f"--prot_seq=uniprot_ascomycota-filtered-reviewed.fasta \n"
        f"--cores 8 --fungus --AUGUSTUS_ab_initio --softmasking --gff3 --workingdir={braker_output}"
    )
    run_command(braker_command)

#? Step 6: BUSCO - qc for annotated proteomes
for file in ["TA001", "TA002", "TA003"]:
    annotated_genome = os.path.join(BRAKER2_DIR, f"braker2_output_{file}/augustus.hints.gff3")
    busco_output = os.path.join(QUAST_BUSCO_DIR, f"busco_output_{file}_annotations")

    busco_command = f"busco -i {annotated_genome} -o {busco_output} -l sordariomycetes_odb10 -m genome --cpu 4 -f"
    run_command(busco_command)

#? Step 7: OrthoFinder - Searching single-copy orthologs
orthofinder_input = os.path.join(PHYLOGENOMICS_DIR, "orthofinder_input")
orthofinder_output = os.path.join(PHYLOGENOMICS_DIR, "orthofinder_output")
os.makedirs(orthofinder_input, exist_ok=True)
os.makedirs(orthofinder_output, exist_ok=True)

# Assuming you have multiple proteomes to analyze
# Copy all proteomes to the ORTHOFINDER_INPUT directory
proteomes_path = "/path/to/proteomes/*.faa"
run_command(f"cp {proteomes_path} {orthofinder_input}")

orthofinder_command = f"orthofinder -f {orthofinder_input} -t 8 -o {orthofinder_output}"
run_command(orthofinder_command)

#? Step 8: Phylogenomics analysis
for file in os.listdir(orthofinder_output):
    if file.endswith(".fa"):
        file_no_ext = os.path.splitext(file)[0]
        #Muscle
        muscle_command = f"muscle -align {os.path.join(orthofinder_output, file)} -output {os.path.join(PHYLOGENOMICS_DIR, 'muscle5', f'{file_no_ext}_aln.fa')}"
        run_command(muscle_command)
        #trimAl
        trimal_command = f"trimal -in {os.path.join(PHYLOGENOMICS_DIR, 'muscle5', f'{file_no_ext}_aln.fa')} -out {os.path.join(PHYLOGENOMICS_DIR, 'trimal', f'{file_no_ext}_trim.fa')} -automated1"
        run_command(trimal_command)

#IQ-TREE2
iqtree_command = f"iqtree2 -s {os.path.join(PHYLOGENOMICS_DIR, 'iqtree2', 'concat.fa')} -m MFP -B 1000 -T AUTO"
run_command(iqtree_command)

#? Step 9: Functional annotations
for file in os.listdir(BRAKER2_DIR):
    annotated_genome = os.path.join(BRAKER2_DIR, file, "augustus.hints.gff3")
    #antismash
    antismash_command = f"antismash --input-dir {annotated_genome} --output-dir {os.path.join(FUNCTIONAL_ANNOTATIONS_DIR, file)} --taxon fungi"
    run_command(antismash_command)
    #dbcan
    dbcan_command = f"run_dbcan --input-dir {annotated_genome} --output-dir {os.path.join(FUNCTIONAL_ANNOTATIONS_DIR, file)} --db_dir /user/work/dissertation_project/dbcan_db/ --tools diamond --dbcan_threads 8"
    run_command(dbcan_command)

# Separate Processing for Downloaded fasta.gz Files
for file in os.listdir(DOWNLOAD_DIR):
    if file.endswith(".fasta.gz"):
        basename = os.path.splitext(file)[0]
        unzipped_file = os.path.join(DOWNLOAD_DIR, f"{basename}.fasta")
        
        # Unzip the fasta.gz file
        gunzip_command = f"gunzip -c {os.path.join(DOWNLOAD_DIR, file)} > {unzipped_file}"
        run_command(gunzip_command)
        
        # Run QUAST
        quast_output = os.path.join(QUAST_BUSCO_DIR, f"quast_output_{basename}")
        quast_command = f"quast.py {unzipped_file} -o {quast_output}"
        run_command(quast_command)
        
        # Run RepeatModeler and RepeatMasker
        repeat_output = os.path.join(REPEAT_DIR, f"repeat_output_{basename}")
        os.makedirs(repeat_output, exist_ok=True)
        os.chdir(repeat_output)
        
        build_db_command = f"BuildDatabase -name {basename} {unzipped_file}"
        run_command(build_db_command)
        
        repeat_modeler_command = f"RepeatModeler -database {basename} -pa 8"
        run_command(repeat_modeler_command)
        
        repeat_masker_command = f"RepeatMasker -gff -pa 8 -lib {os.path.join(repeat_output, 'consensi.fa.classified')} {unzipped_file}"
        run_command(repeat_masker_command)
