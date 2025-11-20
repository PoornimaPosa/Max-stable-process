HPC Workflow for Spatial/Spatiotemporal Max-Stable Modelling in R

This repository contains high-performance R code designed to fit spatial and spatiotemporal max-stable process (MSP) models across multiple climate regions and time windows.
The workflow uses:

Parallel computing with future + furrr

Large-scale model fitting using SpatialExtremes

Automated extraction of return levels, TIC, and F-madograms

Batch processing over 12 regions and 20 moving windows

Because the computation is heavy, this code is intended to be run on an HPC (High-Performance Computing) cluster.

HPC Requirements

R environment
Make sure your cluster has R ≥ 4.0 with the following packages installed:

SpatialExtremes
furrr
future
tidyverse
openxlsx
readxl

You may need to load modules such as:

module load R/4.3.0
module load gdal
module load proj


(Modules vary between clusters.)

Running the Code on an HPC Cluster
Step 1. Edit the code paths

At the top of the script, set your correct working directory:

base_path <- "/path/to/project/"

Step 2. Adjust number of HPC workers

In the main script:

plan(multisession, workers = 12)


Change 12 to the number of CPU cores your SLURM allocation provides.

Submitting a SLURM Job

Each cluster has a different SLURM setup, so you must write your own SLURM submission file.
Below is a template you can adapt.

Create a file, e.g., run_msp.slurm:

#!/bin/bash
#SBATCH --job-name=msp_run
#SBATCH --output=msp_run.log
#SBATCH --error=msp_run.err
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00

module load R/4.3.0

Rscript your_script.R

Modify according to your cluster:

partition

cpus-per-task

mem

time limit

module names
Different HPC systems use different configurations.

Submit the job:

sbatch run_msp.slurm

Output Files

The script automatically generates:

Return levels
results/RL_s_REGION.xlsx
results/RL_st_REGION.xlsx

TIC summaries
results/TIC_summary_s.xlsx
results/TIC_summary_st.xlsx

F-madogram plots
fmadograms/fmado_TYPE_REGION_WINDOW.png

Important Notes

The script is designed to run entirely in parallel, accelerating ~240 model fits.

The user must customise the SLURM file for their own HPC environment.

For very large jobs, setting future::plan(cluster) with HPC node lists is also possible.
