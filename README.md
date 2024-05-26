# msc_dissertationproject
Fungal WGS ONT bioinformatics pipeline

Bioinformatics Pipeline - Fungal WGS ONT Bioinformatics Pipeline

#Overview

This repository contains scripts for running a comprehensive bioinformatics pipeline for fungal whole-genome sequencing using Oxford Nanopore Technology (ONT). The pipeline performs quality control&filtering, whole genome assembly, genome annotation, phylogenetic analysis, and functional annotations.

#Python Script

pipeline.py: The main Python script that sequentially executes the entire bioinformatics pipeline that can be run using Python

pipeline.slurm: The integrated SLURM script that executes all steps of the bioinformatics pipeline that can be run in terminal

- Quality Control: Perform quality control using NanoPlot.
- Data Preprocessing: Preprocess raw data using Chopper.
- Genome Assembly: Assemble genomes using Flye.
- Repeat Identification and Masking: Identify and mask repeats using RepeatModeler and RepeatMasker.
- Genome Annotation: Annotate genomes using BRAKER2.
- Quality Assessment: Assess genome assembly quality using QUAST.
- Genome Completeness: Evaluate genome completeness using BUSCO.
- Orthologous Gene Identification: Identify orthologous genes using OrthoFinder.
- Phylogenetic Analysis: Perform phylogenetic analysis using MUSCLE, trimAl, and IQ-TREE.
- Functional Annotations: Conduct functional annotations using AntiSMASH and dbCAN.
- Processing Publicly Available FASTA Files: Process downloaded FASTA files, including quality assessment, repeat masking, annotation, and completeness evaluation.

#Usage

Ensure you have the required software and dependencies installed. This pipeline is designed to run on a SLURM-managed cluster or using a Python environment with the necessary bioinformatics tools installed.

#Prerequisites
- SLURM
- Python 3.x
- Conda environments for specific bioinformatics tools
- Bioinformatics tools: NanoPlot, Chopper, Flye, QUAST, BUSCO, BRAKER2, RepeatModeler,  RepeatMasker, OrthoFinder, MUSCLE, trimAl, IQ-TREE, AntiSMASH, dbCAN

# Running the Pipeline
- Set Up Directories: Ensure all necessary directories are created. The script will handle this automatically.
- Prepare Input Files: Place your input files in the appropriate directories.
Running the Python Script

To run the entire pipeline using Python:

$python3 pipeline.py
Running the SLURM Script

To run the entire pipeline, submit the integrated SLURM script:

$sbatch scripts/pipeline.slurm
This README.md file provides a clear overview of the pipeline, the steps involved, and instructions on how to run the scripts.





