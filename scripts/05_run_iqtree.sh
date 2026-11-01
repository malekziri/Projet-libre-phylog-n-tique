#!/bin/bash
set -e

cd "$(dirname "$0")/.."

IQTREE="/Users/malekziri/Documents/M2_Bioinformatique/BQ4CE075_Bioinformatique_des_omiques_et_de_la_métagénomique/Phylogénétique/tools/iqtree-2.4.0-macOS/bin/iqtree2"

PATIENTS=("A10" "A17" "A21" "A24" "A31" "A34")
TYPES=("SNV" "INDEL" "REARRANGEMENT" "SNV_INDEL_REARRANGEMENT")

for type in "${TYPES[@]}"; do
    mkdir -p "Results/IQTREE/${type}"

    for patient in "${PATIENTS[@]}"; do
        prefix="${patient}_${type}"
        input="Data/processed/${type}/${prefix}.nex"
        output="Results/IQTREE/${type}/${prefix}"

        echo "Analyse ${patient} — ${type}"

        "$IQTREE" \
            -s "$input" \
            -st BIN \
            -keep-ident \
            -m MFP \
            -o Normal \
            -B 1000 \
            -T 1 \
            -seed 2026 \
            --prefix "$output"
    done
done