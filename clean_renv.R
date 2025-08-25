# Step 0: Install core bootstrapping packages globally (renv + BiocManager)
install.packages(c("renv", "BiocManager"), repos = "https://cloud.r-project.org")

# Step 1: Clean isolated cache (optional)
Sys.setenv(RENV_PATHS_CACHE = tempfile())

# Step 2: Set Bioconductor 3.22 globally
Sys.setenv(RENV_CONFIG_BIOCONDUCTOR_VERSION = "3.22")
options(repos = BiocManager::repositories(version = "3.22"))

# Step 3: Start a clean renv environment
renv::init(bare = TRUE)

# Step 4: Install packages (CRAN + BioC)
renv::install("bioc::DESeq2")
renv::install("ggplot2")

# Step 5: Snapshot
renv::snapshot(force = TRUE)
