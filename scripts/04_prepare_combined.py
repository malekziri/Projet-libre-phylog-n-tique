from pathlib import Path

import pandas as pd


BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "Data/processed"
OUTPUT_DIR = DATA_DIR / "SNV_INDEL_REARRANGEMENT"

TYPES = ["SNV", "INDEL", "REARRANGEMENT"]

SAMPLES = {
    "A10": ["A", "C", "D", "E"],
    "A17": ["A", "C", "D", "E", "F"],
    "A21": ["A", "C", "D", "E", "F", "G", "H", "I", "J"],
    "A24": ["A", "C", "D", "E"],
    "A31": ["A", "C", "D", "E", "F"],
    "A34": ["A", "C", "D", "E"],
}

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

tables = {
    event_type: pd.read_csv(
        DATA_DIR / event_type / f"{event_type}_presence.csv"
    )
    for event_type in TYPES
}

summaries = []
presence_tables = []

for patient, letters in SAMPLES.items():
    samples = [f"{patient}-{letter}" for letter in letters] + ["Normal"]

    matrices = []
    event_counts = {}
    offset = 0

    for event_type in TYPES:
        data = tables[event_type]
        data = data[data["patient"] == patient].copy()

        sites = (
            data[["site_number", "mutation"]]
            .drop_duplicates()
            .sort_values("site_number")
        )

        matrix = data.pivot(
            index="sample",
            columns="mutation",
            values="presence",
        ).reindex(index=samples, columns=sites["mutation"])

        matrix = matrix.astype(int)
        matrix.columns = [
            f"{event_type}_{mutation}" for mutation in matrix.columns
        ]
        matrices.append(matrix)

        count = matrix.shape[1]
        event_counts[event_type] = count

        data["type"] = event_type
        data["mutation"] = event_type + "_" + data["mutation"]
        data["site_number"] = data["site_number"] + offset

        presence_tables.append(
            data[
                [
                    "patient",
                    "sample",
                    "site_number",
                    "type",
                    "mutation",
                    "annotation",
                    "presence",
                ]
            ]
        )

        offset += count

    combined = pd.concat(matrices, axis=1)
    tumors = combined.drop(index="Normal")

    present = tumors.sum(axis=0)
    absent = len(letters) - present

    summaries.append({
        "patient": patient,
        "tumor_samples": len(letters),
        "snv_count": event_counts["SNV"],
        "indel_count": event_counts["INDEL"],
        "rearrangement_count": event_counts["REARRANGEMENT"],
        "total_sites": combined.shape[1],
        "variable_sites": int(((present > 0) & (absent > 0)).sum()),
        "informative_sites": int(((present >= 2) & (absent >= 2)).sum()),
    })

    nexus_file = OUTPUT_DIR / f"{patient}_SNV_INDEL_REARRANGEMENT.nex"

    with nexus_file.open("w", encoding="utf-8") as handle:
        handle.write("#NEXUS\n\nBEGIN DATA;\n")
        handle.write(
            f"    DIMENSIONS NTAX={combined.shape[0]} "
            f"NCHAR={combined.shape[1]};\n"
        )
        handle.write('    FORMAT DATATYPE=STANDARD SYMBOLS="01";\n')
        handle.write("    MATRIX\n")

        for sample, row in combined.iterrows():
            sequence = "".join(row.astype(str))
            handle.write(f"    {sample}    {sequence}\n")

        handle.write("    ;\nEND;\n")

pd.DataFrame(summaries).to_csv(
    OUTPUT_DIR / "SNV_INDEL_REARRANGEMENT_summary.csv",
    index=False,
    encoding="utf-8",
)

pd.concat(presence_tables, ignore_index=True).to_csv(
    OUTPUT_DIR / "SNV_INDEL_REARRANGEMENT_presence.csv",
    index=False,
    encoding="utf-8",
)

print(f"Combinaison terminée pour {len(SAMPLES)} patients.")
print(f"Résultats : {OUTPUT_DIR}")