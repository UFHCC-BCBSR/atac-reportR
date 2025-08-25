# ATAC-seq Report Generator (HiPerGator Version)

This Shiny app generates an interactive HTML report for differential chromatin accessibility using output from the **nf-core/atacseq** pipeline.

It supports automated normalization, filtering, and contrast-based differential analysis using `edgeR`, and produces summary tables, volcano plots, contrast comparisons, enrichment analyses, and trackhub export tools. 

The sample sheet can be downloaded, edited, and uploaded to customize which samples are included, how metadata are grouped, and what comparisons are made. 

Optional fields allow the report to be customized for your project.
---

## What You Need

To run the app, you'll need:

- Access to **HiPerGator**
- Output from an `nf-core/atacseq` run on HiPerGator copied to `/blue/cancercenter-dept/privapps/data/atac/<seqID>`:
	- To copy files, do `bash /blue/cancercenter-dept/privapps/data/atac/retrieve-atac-results.bash --output <nfcore-output-dir> --destination /blue/cancercenter-dept/privapps/data/atac/<seqID>`
---

## How to Use it on Hipergator (for development etc.)

### 1. Clone this repository anywhere on HiPerGator:

```bash
https://github.com/UFHCC-BCBSR/atac-reportR.git
cd atac-reportR
```

### 2. : Connect to reserver to use Rstudio on HiPerGator (in a dev session or with SLURM sbatch)

```bash
moduler load R
rserver
```

For help with this, see https://docs.rc.ufl.edu/software/apps/r/rstudio_server/

### 3. : Use "Run App" to run the app
---

## 📝 Sample Sheet and Parameters

- If the sample sheet is not provided, the app will automatically use:
  ```
  <nf-core output dir>/pipeline_info/samplesheet.valid.csv
  ```
- You may **download and modify the sample sheet** to group, exclude, or rename samples.
- All input parameters (contrasts, filtering thresholds, organism, etc.) can be customized in the interface.

---

## 📖 Documentation

More details on how to use the app, the meaning of each parameter, and tips for editing the sample sheet are provided **within the app UI**.

---

## 📦 Technical Details

This app runs entirely inside a container built from:
- R + Bioconductor packages
- Shiny + bs4Dash
- The `atac_edgeR_report.Rmd` analysis script

---

## Maintainer

Heather Kates  
📧 hkates@ufl.edu  
UF Health Cancer Center – BCB-SR Bioinformatics Analyst

