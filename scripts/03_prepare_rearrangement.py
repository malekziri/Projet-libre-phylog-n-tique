from pathlib import Path

import pandas as pd


BASE_DIR = Path(__file__).resolve().parent.parent
INPUT_FILE = BASE_DIR / "Data/raw/41586_2020_2581_MOESM1_ESM.xlsx"
OUTPUT_DIR = BASE_DIR / "Data/processed/REARRANGEMENT"

SAMPLES = {
    "A10": ["A", "C", "D", "E"],
    "A17": ["A", "C", "D", "E", "F"],
    "A21": ["A", "C", "D", "E", "F", "G", "H", "I", "J"],
    "A24": ["A", "C", "D", "E"],
    "A31": ["A", "C", "D", "E", "F"],
    "A34": ["A", "C", "D", "E"],
}

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

all_data = pd.read_excel(INPUT_FILE, sheet_name="rearrangements")

summaries = []
presence_tables = []

for patient, letters in SAMPLES.items():
    samples = [f"{patient}-{letter}" for letter in letters]
    data = all_data[all_data["patient"] == patient].copy()

    events = {}
    annotations = {}

    for _, row in data.iterrows():
        mutation = (
            f"{row['chr1']}:{int(row['pos1'])}:{row['str1']}"
            f"|{row['chr2']}:{int(row['pos2'])}:{row['str2']}"
        )

        if mutation not in events:
            events[mutation] = {sample: 0 for sample in samples}

        annotations[mutation] = row["paired-flag"]

        sample_letters = str(row["samples"]).upper().split("/")
        counts = str(row["tumour-counts"]).split("/")

        for letter, count in zip(sample_letters, counts):
            sample = f"{patient}-{letter.strip()}"
            support = pd.to_numeric(count, errors="coerce")

            if sample in events[mutation] and pd.notna(support) and support > 0:
                events[mutation][sample] = 1

    matrix = pd.DataFrame(events, index=samples).astype(int)
    matrix = matrix.reindex(sorted(matrix.columns), axis=1)

    present = matrix.sum(axis=0)
    absent = len(samples) - present

    summaries.append({
        "patient": patient,
        "tumor_samples": len(samples),
        "rearrangement_count": matrix.shape[1],
        "variable_sites": int(((present > 0) & (absent > 0)).sum()),
        "informative_sites": int(((present >= 2) & (absent >= 2)).sum()),
    })

    matrix.loc["Normal"] = 0

    site_numbers = {
        mutation: number
        for number, mutation in enumerate(matrix.columns, start=1)
    }

    presence = (
        matrix.rename_axis(index="sample", columns="mutation")
        .reset_index()
        .melt(
            id_vars="sample",
            var_name="mutation",
            value_name="presence",
        )
    )

    presence.insert(0, "patient", patient)
    presence["site_number"] = presence["mutation"].map(site_numbers)
    presence["annotation"] = presence["mutation"].map(annotations)

    presence_tables.append(
        presence[
            [
                "patient",
                "sample",
                "site_number",
                "mutation",
                "annotation",
                "presence",
            ]
        ]
    )

    nexus_file = OUTPUT_DIR / f"{patient}_REARRANGEMENT.nex"

    with nexus_file.open("w", encoding="utf-8") as handle:
        handle.write("#NEXUS\n\nBEGIN DATA;\n")
        handle.write(
            f"    DIMENSIONS NTAX={matrix.shape[0]} "
            f"NCHAR={matrix.shape[1]};\n"
        )
        handle.write('    FORMAT DATATYPE=STANDARD SYMBOLS="01";\n')
        handle.write("    MATRIX\n")

        for sample, row in matrix.iterrows():
            sequence = "".join(row.astype(str))
            handle.write(f"    {sample}    {sequence}\n")

        handle.write("    ;\nEND;\n")

pd.DataFrame(summaries).to_csv(
    OUTPUT_DIR / "REARRANGEMENT_summary.csv",
    index=False,
    encoding="utf-8",
)

pd.concat(presence_tables, ignore_index=True).to_csv(
    OUTPUT_DIR / "REARRANGEMENT_presence.csv",
    index=False,
    encoding="utf-8",
)

print(f"Préparation terminée pour {len(SAMPLES)} patients.")
print(f"Résultats : {OUTPUT_DIR}")