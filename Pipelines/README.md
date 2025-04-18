# --- ILLUMINA PIPELINE ---

## Overview
The pipeline processes raw data generated using Illumina technology and produces two output files:

1. **Gene Expression Matrix (`file.csv`)**  
   A matrix containing samples, genes, and their corresponding gene expression levels.

2. **Sample-Diagnosis Mapping (`file_mapping.csv`)**  
   A file mapping each patient (sample) to their diagnosis (e.g., ASD or control).

## Output Location
Both files are saved in the `Dataset` directory under their respective subfolders.

## How to Run the Pipeline

### Steps:
1. Place the raw dataset in the appropriate folder under `Datasets/<technology> <dataset_id>`.  
   For example, for Illumina technology and dataset `GSE42133`, place the raw data in:  
   `Datasets/Illumina/GSE42133`.

2. Navigate to the root of the project directory in your terminal.

3. Execute the pipeline using the following command:
   ```bash
   Rscript Pipelines/<technology> <dataset_id>