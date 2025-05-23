#!/bin/bash
#SBATCH --job-name=shinyapp
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=atacapp_%j.out
#SBATCH --account=timgarrett

# === Load Apptainer ===
module load apptainer

# === Choose port ===
PORT=7291

# === Find current node name ===
NODE=$(hostname)

# === Print SSH tunnel instructions ===
echo "======================================================"
echo "Shiny app is launching on HiPerGator node: $NODE"
echo ""
echo "To access it from your local browser, run this command"
echo "in a new terminal on your local machine:"
echo ""
echo "ssh -N -L ${PORT}:${NODE}:${PORT} ${USER}@hpg.rc.ufl.edu"
echo ""
echo "Then open your browser and go to:"
echo "http://localhost:${PORT}"
echo "======================================================"

APPDIR=$(pwd)
export APPDIR

# === Run the Shiny app inside the container ===
apptainer exec \
  --bind "$(pwd)":/app \
  --bind /blue/cancercenter-dept/CONTAINERS/shared_R_lib:/usr/local/lib/R/site-library \
  --env HOST_APPDIR="$APPDIR" \
  /blue/cancercenter-dept/CONTAINERS/atacseq_app.sif \
  Rscript -e "shiny::runApp('/app', host = '0.0.0.0', port = ${PORT})" > atacReportApp_Render_${SLURM_JOB_ID}.log 2>&1

