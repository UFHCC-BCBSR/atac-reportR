# === CONFIG ===
PUBLIC_WEB_DIR="/orange/cancercenter-dept/web/public/BCB-SR"
PUBLIC_WEB_URL="https://data.rc.ufl.edu/pub/cancercenter-dept/BCB-SR"

# === Organism setup ===
ORGANISM="${organism}" 
if [[ "$ORGANISM" == "mmu" ]]; then
  GENOME_ALIAS="mm10"
else
  GENOME_ALIAS="hg38"
fi

# === Required environment variables ===
CONTRAST="${contrast}"
SEQID="${seqID}"

# === Setup paths ===
TRACKHUB_DIR="${PUBLIC_WEB_DIR}/trackhubs/${SEQID}/${CONTRAST}/${GENOME_ALIAS}"
TEMP_DIR="${PUBLIC_WEB_DIR}/igv/${SEQID}/${CONTRAST}_igv_tmp"
FINAL_ZIP="${PUBLIC_WEB_DIR}/igv/${SEQID}/${CONTRAST}_igv.zip"
mkdir -p "$TEMP_DIR"

# === Copy all bigWig and BigBed files from trackhub dir ===
cp "$TRACKHUB_DIR"/*.bigWig "$TEMP_DIR/" 2>/dev/null || echo "[WARN] No bigWig files found."
cp "$TRACKHUB_DIR"/*.bb "$TEMP_DIR/" 2>/dev/null || echo "[WARN] No BigBed file found."

# === Create zip ===
cd "$TEMP_DIR"
zip -r "$FINAL_ZIP" ./*

# === Set permissions and clean up ===
chmod o+r "$FINAL_ZIP"
rm -rf "$TEMP_DIR"

# === Output location ===
echo "[SUCCESS] IGV session zip created at:"
echo "${PUBLIC_WEB_URL}/${SEQID}/${CONTRAST}_igv.zip"
