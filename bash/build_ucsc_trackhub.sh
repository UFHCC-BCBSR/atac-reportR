#!/bin/bash

set -euo pipefail

# === CONFIG ===
PUBLIC_WEB_DIR="/orange/cancercenter-dept/web/public/BCB-SR/trackhubs"
PUBLIC_WEB_URL="https://data.rc.ufl.edu/pub/cancercenter-dept/BCB-SR/trackhubs/"
# ==============

# === Organism setup ===
ORGANISM="${organism}" 
if [[ "$ORGANISM" == "mmu" ]]; then
  GENOME_ALIAS="mm10"
  CHROM_SIZES="${CHROM_SIZES:-mm10.chrom.sizes}"
else
  GENOME_ALIAS="hg38"
  CHROM_SIZES="${CHROM_SIZES:-hg38.chrom.sizes}"
fi

# === Ensure chrom.sizes file exists ===
if [[ ! -f "$CHROM_SIZES" ]]; then
  echo "[INFO] Downloading $CHROM_SIZES..."
  wget "http://hgdownload.soe.ucsc.edu/goldenPath/${GENOME_ALIAS}/bigZips/${GENOME_ALIAS}.chrom.sizes" -O "$CHROM_SIZES"
fi


# === Required environment variables ===
CONTRAST="${contrast}"
SEQID="${seqID}"

# Optional mapping from group to sample IDs
GROUP_SAMPLE_MAP="${group_sample_map:-}"

# === Groups parsed from contrast ===
GROUP1=$(echo "$CONTRAST" | awk -F'_vs_' '{print $1}')
GROUP2=$(echo "$CONTRAST" | awk -F'_vs_' '{print $2}')

# === Setup trackhub directory ===
TRACKHUB_DIR="${PUBLIC_WEB_DIR}/${SEQID}/${CONTRAST}"
mkdir -p "${TRACKHUB_DIR}/${GENOME_ALIAS}"
TRACKDB="${TRACKHUB_DIR}/${GENOME_ALIAS}/trackDb.txt"
echo "" > "$TRACKDB"

# === Create hub.txt and genomes.txt ===
cat <<EOF > "${TRACKHUB_DIR}/hub.txt"
hub ${CONTRAST}_hub
shortLabel ATACseq ${CONTRAST}
longLabel ATACseq Differential Accessibility: ${CONTRAST}
genomesFile genomes.txt
email hkates@ufl.edu
EOF

cat <<EOF > "${TRACKHUB_DIR}/genomes.txt"
genome ${GENOME_ALIAS}
trackDb ${GENOME_ALIAS}/trackDb.txt
EOF

# === Move BigBed file from bigbed_files/ and add to trackDb ===
BB_SOURCE="${TRACKHUB_DIR}/${CONTRAST}_DA_peaks.bb"
BB_DEST="${TRACKHUB_DIR}/${GENOME_ALIAS}/${CONTRAST}_DA_peaks.bb"

if [[ -f "$BB_SOURCE" ]]; then
  echo "[INFO] Found BigBed at: $BB_SOURCE"
  mv "$BB_SOURCE" "$BB_DEST"
  cat <<EOF >> "$TRACKDB"
track da_${CONTRAST}
type bigBed 6
shortLabel DA Peaks
longLabel Differential Accessibility Peaks for ${CONTRAST}
visibility dense
color 255,0,0
priority 1
bigDataUrl $(basename "$BB_DEST")
EOF
else
  echo "[ERROR] BigBed file not found at: $BB_SOURCE"
  exit 1
fi

# === Add bigWig tracks ===
declare -A GROUP_COLORS
GROUP_COLORS["$GROUP1"]="0,0,255"
GROUP_COLORS["$GROUP2"]="255,0,0"

for group in "$GROUP1" "$GROUP2"; do
  color="${GROUP_COLORS[$group]}"
  bigwig_files=""

  if [[ -n "$GROUP_SAMPLE_MAP" ]]; then
    # Extract sample names for this group
    sample_list=$(echo "$GROUP_SAMPLE_MAP" | tr ' ' '\n' | grep "^${group}:" | cut -d: -f2 | tr ',' ' ')

    for sample in $sample_list; do
      if [[ -n "${pipeline_output:-}" ]]; then
      bw_path=$(ls "$pipeline_output"/*${sample}*.bigWig 2>/dev/null || true)
      elif [[ -n "${path_to_bigwigs:-}" ]]; then
      bw_path=$(ls "${path_to_bigwigs}"/*${sample}*.bigWig 2>/dev/null || true)
      fi
      bigwig_files="$bigwig_files $bw_path"
    done
  else
    echo "[WARN] group_sample_map not provided; falling back to wildcard search for group: $group"
    if [[ -n "${pipeline_output:-}" ]]; then
      bigwig_files=$(ls "$pipeline_output"/*${group}*.bigWig 2>/dev/null || true)
    elif [[ -n "${path_to_bigwigs:-}" ]]; then
      bigwig_files=$(ls "${path_to_bigwigs}"/*${group}*.bigWig 2>/dev/null || true)
    else
      echo "[ERROR] Must provide 'pipeline_output' or 'path_to_bigwigs'"
      exit 1
    fi
  fi

  for bw in $bigwig_files; do
    [[ -f "$bw" ]] || continue
    bw_name=$(basename "$bw")
    track_name="${bw_name%.bigWig}"
    echo "[INFO] Adding bigWig: $bw_name"
    cp "$bw" "${TRACKHUB_DIR}/${GENOME_ALIAS}/"

    cat <<EOF >> "$TRACKDB"

track $track_name
type bigWig
shortLabel $track_name
longLabel bigWig track for $track_name ($group)
visibility full
autoScale on
maxHeightPixels 100:60:8
color $color
bigDataUrl $bw_name
EOF
  done
done

# === change permissions of public web directory ===
chmod -R o+rX "$TRACKHUB_DIR"

# === Output final UCSC link ===
echo
echo "[SUCCESS] Trackhub deployed!"
echo "UCSC Genome Browser Link:"

echo "https://genome.ucsc.edu/cgi-bin/hgTracks?hubUrl=${PUBLIC_WEB_URL}$SEQID/$CONTRAST/hub.txt"
