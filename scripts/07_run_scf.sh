#!/bin/bash
set -e

cd "$(dirname "$0")/.."

IQTREE="/Users/malekziri/Documents/M2_Bioinformatique/BQ4CE075_Bioinformatique_des_omiques_et_de_la_métagénomique/Phylogénétique/tools/iqtree-2.4.0-macOS/bin/iqtree2"

PATIENTS=("A10" "A17" "A21" "A24" "A31" "A34")
TYPES=("SNV" "INDEL" "REARRANGEMENT" "SNV_INDEL_REARRANGEMENT")

for type in "${TYPES[@]}"; do
    mkdir -p "Results/SCF/${type}"

    for patient in "${PATIENTS[@]}"; do
        prefix="${patient}_${type}"
        input="Data/processed/${type}/${prefix}.nex"
        tree="Results/IQTREE/${type}/${prefix}.treefile"
        report="Results/IQTREE/${type}/${prefix}.iqtree"
        output="Results/SCF/${type}/${prefix}"

        model=$(awk '/^Best-fit model according to BIC:/ {print $NF; exit}' "$report")

        echo "Calcul du sCF : ${patient} — ${type}"

        "$IQTREE" \
            -te "$tree" \
            -s "$input" \
            -st BIN \
            -m "$model" \
            --scfl 100 \
            -T 1 \
            -seed 2026 \
            --prefix "$output"
    done
done