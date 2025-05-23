# ATAC-seq Report Generator (HiPerGator Version)

This Shiny app generates an interactive HTML report for differential chromatin accessibility using output from the **nf-core/atacseq** pipeline.

It supports automated normalization, filtering, and contrast-based differential analysis using `edgeR`, and produces summary tables, volcano plots, contrast comparisons, enrichment analyses, and trackhub export tools. The sample sheet can be edited to customize which samples are included, how metadata are grouped, and what comparisons are made. Optional metadata fields allow the report to be customized for your project.

---

## ✅ What You Need

To run the app, you'll need:

- Access to **HiPerGator**
- Output from an `nf-core/atacseq` run on HiPerGator
- Read access to the prebuilt container image (SIF file):  
  ```
  /blue/cancercenter-dept/CONTAINERS/atacseq_app.sif
  ```

---

## 🚀 How to Use It

### 1. Clone this repository anywhere on HiPerGator:
```bash
git clone https://github.com/HeatherKates/atacReportApp.git
cd atacReportApp
```

### 2. Submit the job:
```bash
sbatch run_atac_apptainer.sbatch
```

This will:

- Launch the app on a compute node
- Print SSH tunneling instructions in the *out file like:
  ```
  ssh -N -L 7291:c1234:7291 yourGatorUser@hpg.rc.ufl.edu
  ```
- Allow you to access the app at:  
  http://localhost:7291

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

The container is mounted via `apptainer` and does **not require any installation of R or packages on HiPerGator**.

---

## 🛠 Troubleshooting

If you get a message about missing access to the `.sif` file, ensure your HiPerGator account has permission to read from:

```
/blue/cancercenter-dept/CONTAINERS/atacseq_app.sif
```

Contact hkates@ufl.edu if access needs to be granted.

---

## Maintainer

Heather Kates  
📧 hkates@ufl.edu  
UF Health Cancer Center – BCB-SR Bioinformatics Analyst

