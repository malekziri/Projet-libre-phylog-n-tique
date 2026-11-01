#!/bin/bash
set -e

cd "$(dirname "$0")/.."

MB="mb"
TYPE="SNV_INDEL_REARRANGEMENT"
INPUT_DIR="Data/processed/${TYPE}"
OUTPUT_DIR="Results/MRBAYES_10M/${TYPE}"

PATIENTS=("A10" "A17" "A21" "A24" "A31" "A34")

mkdir -p "$OUTPUT_DIR"

for patient in "${PATIENTS[@]}"; do
    prefix="${patient}_${TYPE}"
    input="${INPUT_DIR}/${prefix}.nex"
    mb_input="${OUTPUT_DIR}/${prefix}_MrBayes.nex"

    python3 - "$input" "$mb_input" "$prefix" <<'PY'
from pathlib import Path
import re
import sys

input_file = Path(sys.argv[1])
output_file = Path(sys.argv[2])
prefix = sys.argv[3]

source = input_file.read_text(encoding="utf-8")

matrix = re.search(
    r"\bMATRIX\b(.*?);",
    source,
    flags=re.IGNORECASE | re.DOTALL,
).group(1)

sequences = [
    line.split()[1]
    for line in matrix.splitlines()
    if line.strip()
]

has_constant_sites = any(
    len(set(column)) == 1
    for column in zip(*sequences)
)

coding = "ALL" if has_constant_sites else "VARIABLE"

source = re.sub(
    r"DATATYPE\s*=\s*STANDARD[^;]*;",
    "DATATYPE=RESTRICTION MISSING=?;",
    source,
    count=1,
    flags=re.IGNORECASE,
)

source += f"""
BEGIN MRBAYES;
    SET AUTOCLOSE=YES SEED=2026 SWAPSEED=2026;
    LSET CODING={coding} RATES=EQUAL;
    PRSET STATEFREQPR=DIRICHLET(1.0,1.0);
    OUTGROUP Normal;
    MCMCP NGEN=10000000 NRUNS=2 NCHAINS=4 SAMPLEFREQ=1000
          PRINTFREQ=1000 DIAGNFREQ=10000 SAVEBRLENS=YES
          FILENAME={prefix};
    MCMC;
    SUMP RELBURNIN=YES BURNINFRAC=0.25;
    SUMT RELBURNIN=YES BURNINFRAC=0.25;
END;
"""

output_file.write_text(source, encoding="utf-8")
PY

    echo "Analyse MrBayes : ${patient}"

    (
        cd "$OUTPUT_DIR"
        "$MB" "${prefix}_MrBayes.nex"
    )
done