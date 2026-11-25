from pathlib import Path
import gzip

import numpy as np
import pandas as pd


project = Path(".")
reference_dir = project / "01_input/reference"
symbol_file = Path(__file__).with_name("reference_gene_symbols.csv")
result_dir = project / "03_results/04_source_tracing"
output_dir = result_dir / "female_reference"
p_cutoff = 0.05
log2fc_cutoff = 0.5

raw = pd.read_excel(
    reference_dir / "GSE219045_Rat_FPKM_table_Mar22.xlsx",
    sheet_name="Rat_FPKM_table_March22",
).set_index("GeneID")
symbols = pd.read_csv(symbol_file, dtype=str).set_index("Geneid")["Gene symbol"]
de = pd.read_csv(project / "03_results/02_DE/tables/DE_all.csv")
assert raw.index.is_unique and symbols.index.is_unique and de.feature_id.is_unique
assert not raw.isna().any().any() and (raw >= 0).all().all()

with gzip.open(reference_dir / "GSE219045_family.soft.gz", "rt") as handle:
    records = handle.read().split("^SAMPLE = ")[1:]
metadata = []
for record in records:
    lines = record.splitlines()
    if "!Sample_organism_ch1 = Rattus norvegicus" not in lines:
        continue
    if "!Sample_characteristics_ch1 = Sex: Female" not in lines:
        continue
    title = next(s.split(" = ", 1)[1] for s in lines if s.startswith("!Sample_title = "))
    column = title.removeprefix("Rat_")
    if lines[0] == "GSM6765843":
        column = "BrainMedulla_Female"
    metadata.append({
        "GEO accession": lines[0], "Reference column": column,
        "Tissue": column.removesuffix("_Female"), "Sex": "Female",
    })
metadata = pd.DataFrame(metadata).sort_values("Tissue")
female_columns = metadata["Reference column"].tolist()
assert len(female_columns) == 44 and len(set(female_columns)) == 44
assert set(female_columns) == {
    c for c in raw if c.endswith("_Female") or c in ["Ovary", "Uterus", "Vagina"]
}

female = raw[female_columns].rename(columns=dict(zip(female_columns, metadata.Tissue)))

selected = de.loc[
    (de.contrast == "OP_vs_N") & (de.RNA_type == "mRNA")
    & (de.pvalue < p_cutoff) & (de.log2FoldChange.abs() > log2fc_cutoff)
].copy()
selected["Expression change"] = np.where(selected.log2FoldChange > 0, "Upregulated", "Downregulated")
matched = selected.loc[selected.feature_id.isin(female.index)].set_index("feature_id")
assert matched.index.isin(symbols.index).all()
expression = female.loc[matched.index]
zero = expression.max(axis=1).eq(0)
ties = expression.eq(expression.max(axis=1), axis=0).sum(axis=1).gt(1) & ~zero
assert set(expression.index[zero]) == {"ENSRNOG00000057953", "ENSRNOG00000003227"}
assert set(expression.index[ties]) == {"ENSRNOG00000056162"}
assigned_expression = expression.loc[~zero & ~ties]
order = np.argsort(-assigned_expression.to_numpy(), axis=1, kind="stable")
rows = np.arange(len(order))
values = assigned_expression.to_numpy()
tissues = assigned_expression.columns.to_numpy()
assignment = matched.loc[assigned_expression.index, ["Expression change", "log2FoldChange", "pvalue", "padj"]].copy()
assignment.insert(0, "Gene symbol", symbols.reindex(assignment.index))
assignment["Highest-expression tissue"] = tissues[order[:, 0]]
assignment["Highest FPKM"] = values[rows, order[:, 0]]
assignment["Second-highest tissue"] = tissues[order[:, 1]]
assignment["Second-highest FPKM"] = values[rows, order[:, 1]]
assignment.index.name = "Geneid"

counts = pd.crosstab(assignment["Highest-expression tissue"], assignment["Expression change"])
counts = counts.reindex(index=female.columns, columns=["Upregulated", "Downregulated"], fill_value=0)
counts.index.name = "Tissue"
counts = counts.sort_values(["Upregulated", "Downregulated"], ascending=False)
unassigned = matched.loc[zero | ties, ["Expression change"]].join(expression)
unassigned.insert(0, "Gene symbol", symbols.reindex(unassigned.index))
unassigned.index.name = "Geneid"
bone_ids = assignment.index[assignment["Highest-expression tissue"].isin(["BoneMarrow", "Cartilage"])]
bone = assignment.loc[bone_ids, ["Gene symbol", "Expression change", "Highest-expression tissue"]].join(female)

summary = []
for change in ["Upregulated", "Downregulated"]:
    current = assignment.loc[assignment["Expression change"] == change]
    summary.append({
        "Expression change": change,
        "Selected mRNAs": int((selected["Expression change"] == change).sum()),
        "Matched mRNAs": int((matched["Expression change"] == change).sum()),
        "Assigned mRNAs": len(current),
        "Female bone marrow": int(current["Highest-expression tissue"].eq("BoneMarrow").sum()),
        "Female cartilage": int(current["Highest-expression tissue"].eq("Cartilage").sum()),
    })
summary = pd.DataFrame(summary)
assert len(assignment) + len(unassigned) == len(matched)
assert counts.to_numpy().sum() == len(assignment)

output_dir.mkdir(parents=True)
female.to_csv(output_dir / "female_reference_FPKM.tsv", sep="\t")
metadata.to_csv(output_dir / "reference_metadata.csv", index=False)
assignment.to_csv(output_dir / "tissue_assignment.csv")
counts.to_csv(output_dir / "tissue_counts.csv")
unassigned.to_csv(output_dir / "unassigned_expression.csv")
bone.to_csv(output_dir / "bone_gene_expression.csv")
summary.to_csv(output_dir / "summary.csv", index=False)
print(summary.to_string(index=False))
print(counts.head(8).to_string())
print(f"Results: {output_dir}")
