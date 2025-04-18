# Datasets Structure

The datasets are organized by technology:

- **Affymetrix**
- **Agilent**
- **Illumina**


Within each technology folder, all datasets related to that technology are stored. For example:

- **Illumina**
  - Dataset1: `GSE111`
  - Dataset2: `GSE112`
  - ...

Each dataset folder contains:
1. **Raw Data**: The unprocessed dataset.
2. **Pipeline Output (`pipelineOut`)**: The folder containing the results of the pipeline.


### Pipeline Output Details

1. **Gene Expression Matrix (`gene_expression_matrix.csv`)**:  
   A matrix containing gene expression levels for all patients.

2. **Patient-Diagnosis Mapping (`patient_diagnosis_mapping.csv`)**:  
   A file mapping each patient (sample) to their diagnosis (e.g., ASD or control).
