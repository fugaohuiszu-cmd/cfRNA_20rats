from pathlib import Path
import argparse
import json
import os
import time

THREADS = "64"
for key in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[key] = THREADS

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics.pairwise import cosine_similarity

import CellFreeGMF

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-counts", required=True)
    parser.add_argument("--de-results", required=True)
    parser.add_argument("--celltype-map", required=True)
    parser.add_argument("--celltype-colors", required=True)
    parser.add_argument("--reference-dir", required=True)
    parser.add_argument("--output-dir", default="./03_results/06_GMF_CellChat/output")
    return parser.parse_args()

def require(path, label):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Missing {label}: {path}")
    return path

def aggregate_scores(sample_cell, sample_names, cell_labels, cell_order):
    rows = []
    labels = np.asarray(cell_labels, dtype=str)
    for sample_index, sample_id in enumerate(sample_names):
        for cell_type in cell_order:
            cell_index = np.where(labels == cell_type)[0]
            if len(cell_index) == 0:
                raise ValueError(f"No cells for cell type: {cell_type}")
            rows.append(
                {
                    "sample_id": sample_id,
                    "cell_type": cell_type,
                    "n_reference_cells": int(len(cell_index)),
                    "GMF_score": float(sample_cell[sample_index, cell_index].sum()),
                }
            )
    return pd.DataFrame(rows)

args = parse_args()
raw_counts_path = require(args.raw_counts, "raw count matrix")
de_results_path = require(args.de_results, "DE results")
celltype_map_path = require(args.celltype_map, "11-cell-type mapping")
celltype_colors_path = require(args.celltype_colors, "cell-type colors")
reference_dir = Path(args.reference_dir)
output_dir = Path(args.output_dir)
input_dir = output_dir / "input"
gmf_dir = output_dir / "gmf"
input_dir.mkdir(parents=True, exist_ok=True)
gmf_dir.mkdir(parents=True, exist_ok=True)

full_gene_map_path = require(
    reference_dir
    / "Set7_GSE169396_human_femoral_head_full_cfRNA_mRNA_curated_res0p3_used_gene_mapping.tsv",
    "full rat-human gene mapping",
)
sc_logcpm_path = require(
    reference_dir / "GSE169396_scRNA_logCPM_Set7_cell_by_human_gene.npz",
    "GSE169396 scRNA logCPM matrix",
)
cell_meta_path = require(
    reference_dir / "GSE169396_res0p3_curated_cell_metadata_for_GMF.csv",
    "GSE169396 cell metadata",
)

params = {
    "alpha": 0.001,
    "beta": 0.001,
    "iter_num": 300,
    "random_seed": 2025,
    "threads_requested": int(THREADS),
    "sample_similarity_diagonal": 0,
}
sample_order = [f"N{i}" for i in range(1, 11)] + [f"OP{i}" for i in range(1, 11)]
group_map = {sample_id: ("N" if sample_id.startswith("N") else "OP") for sample_id in sample_order}

de = pd.read_csv(de_results_path)
required_de_columns = {
    "contrast",
    "feature_id",
    "RNA_type",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
}
missing_de_columns = required_de_columns - set(de.columns)
if missing_de_columns:
    raise ValueError(f"DE results missing columns: {sorted(missing_de_columns)}")

selected_de = de.loc[
    (de["contrast"] == "OP_vs_N")
    & (de["RNA_type"] == "mRNA")
    & de["pvalue"].notna()
    & (de["pvalue"] < 0.05)
].copy()
selected_de = selected_de.sort_values(["pvalue", "feature_id"]).reset_index(drop=True)
if selected_de["feature_id"].duplicated().any():
    raise ValueError("Selected DE mRNA feature IDs are duplicated")
if selected_de.shape[0] != 1973:
    raise ValueError(f"Expected 1973 P<0.05 mRNAs, found {selected_de.shape[0]}")

raw = pd.read_csv(raw_counts_path, index_col=0)
required_raw_columns = {"RNA_type", "geneid", *sample_order}
missing_raw_columns = required_raw_columns - set(raw.columns)
if missing_raw_columns:
    raise ValueError(f"Raw matrix missing columns: {sorted(missing_raw_columns)}")
raw.index = raw.index.astype(str)
if raw.index.duplicated().any():
    raise ValueError("Raw matrix feature IDs are duplicated")
if not (raw.index == raw["geneid"].astype(str)).all():
    raise ValueError("Raw matrix row names and geneid column are inconsistent")

raw_mrna = raw.loc[raw["RNA_type"] == "mRNA", sample_order].apply(pd.to_numeric)
library_size = raw_mrna.sum(axis=0)
if (library_size <= 0).any():
    raise ValueError("At least one sample has a non-positive mRNA library size")

full_gene_map = pd.read_csv(full_gene_map_path, sep="\t", dtype=str)
if set(full_gene_map.columns) != {"gene", "rat_gene"}:
    raise ValueError("Full gene mapping must contain gene and rat_gene columns")
if full_gene_map["gene"].duplicated().any() or full_gene_map["rat_gene"].duplicated().any():
    raise ValueError("Full gene mapping is not one-to-one")

selected_rat_genes = set(selected_de["feature_id"].astype(str))
used_gene_map = full_gene_map.loc[full_gene_map["rat_gene"].isin(selected_rat_genes)].copy()
used_gene_map = used_gene_map.loc[used_gene_map["rat_gene"].isin(raw_mrna.index)].reset_index(drop=True)
if used_gene_map.shape[0] != 1593:
    raise ValueError(f"Expected 1593 mapped genes, found {used_gene_map.shape[0]}")

human_gene_order = used_gene_map["gene"].tolist()
rat_gene_order = used_gene_map["rat_gene"].tolist()
cf_raw = raw_mrna.loc[rat_gene_order, sample_order].T
cf_raw.columns = human_gene_order
cf_logcpm = np.log2(cf_raw.div(library_size.loc[sample_order], axis=0) * 1_000_000 + 1)
if not np.isfinite(cf_logcpm.to_numpy()).all() or (cf_logcpm.to_numpy() < 0).any():
    raise ValueError("cfRNA logCPM matrix contains invalid values")

sc_full = sparse.load_npz(sc_logcpm_path).tocsr()
if sc_full.shape[1] != full_gene_map.shape[0]:
    raise ValueError("Full scRNA matrix columns do not match the full gene mapping")
full_gene_index = {gene: i for i, gene in enumerate(full_gene_map["gene"])}
subset_index = [full_gene_index[gene] for gene in human_gene_order]
sc_subset = sc_full[:, subset_index].tocsr()
if sc_subset.shape[1] != cf_logcpm.shape[1]:
    raise ValueError("cfRNA and scRNA gene dimensions do not match")

cell_meta = pd.read_csv(cell_meta_path)
required_cell_columns = {"cell_id", "raw_cluster", "fine_annotation"}
missing_cell_columns = required_cell_columns - set(cell_meta.columns)
if missing_cell_columns:
    raise ValueError(f"Cell metadata missing columns: {sorted(missing_cell_columns)}")
if cell_meta.shape[0] != sc_subset.shape[0]:
    raise ValueError("Cell metadata rows do not match scRNA matrix rows")

celltype_map = pd.read_csv(celltype_map_path)
if set(celltype_map.columns) != {"fine_annotation", "major_celltype_11_proB_split"}:
    raise ValueError("Cell-type mapping has unexpected columns")
celltype_colors = pd.read_csv(celltype_colors_path).sort_values("order")
required_color_columns = {"order", "cell_type_full", "cell_type_short", "umap_color"}
if not required_color_columns.issubset(celltype_colors.columns):
    raise ValueError("Cell-type color table has unexpected columns")

fine_to_full = dict(
    zip(celltype_map["fine_annotation"], celltype_map["major_celltype_11_proB_split"])
)
full_to_short = dict(zip(celltype_colors["cell_type_full"], celltype_colors["cell_type_short"]))
cell_meta["cell_type_full"] = cell_meta["fine_annotation"].map(fine_to_full)
cell_meta["cell_type"] = cell_meta["cell_type_full"].map(full_to_short)
if cell_meta["cell_type"].isna().any():
    missing_labels = sorted(cell_meta.loc[cell_meta["cell_type"].isna(), "fine_annotation"].unique())
    raise ValueError(f"Unmapped cell annotations: {missing_labels}")
cell_order = celltype_colors["cell_type_short"].tolist()
if cell_meta["cell_type"].nunique() != 11:
    raise ValueError("Expected 11 cell types")

selected_de.to_csv(input_dir / "DE_mRNA_pvalue_lt0.05.csv", index=False)
used_gene_map.to_csv(input_dir / "used_gene_map.tsv", sep="\t", index=False)
cf_raw.to_csv(input_dir / "cfrna_raw_counts.tsv", sep="\t")
cf_logcpm.to_csv(input_dir / "cfrna_logCPM.tsv.gz", sep="\t", compression="gzip")
sparse.save_npz(input_dir / "scrna_logCPM.npz", sc_subset, compressed=True)
cell_meta.to_csv(input_dir / "cell_metadata_11.csv", index=False)
celltype_colors.to_csv(input_dir / "celltype_colors.csv", index=False)

input_summary = {
    "DE_rule": "OP vs N; mRNA; P value < 0.05; no log2FC cutoff",
    "n_DE_mRNAs": int(selected_de.shape[0]),
    "n_mapped_human_genes": int(used_gene_map.shape[0]),
    "n_samples": int(cf_logcpm.shape[0]),
    "n_N_samples": 10,
    "n_OP_samples": 10,
    "n_scRNA_cells": int(sc_subset.shape[0]),
    "n_cell_types": 11,
    "cfRNA_logCPM_denominator": "full mRNA library size before gene subsetting",
    "params": params,
}
(input_dir / "input_summary.json").write_text(
    json.dumps(input_summary, indent=2, ensure_ascii=False), encoding="utf-8"
)

cell_matrix = sc_subset.toarray().astype(np.float64, copy=False)
if not np.isfinite(cell_matrix).all() or (cell_matrix < 0).any():
    raise ValueError("scRNA logCPM matrix contains invalid values")

score_tables = []
run_records = []
for group in ("N", "OP"):
    group_samples = [sample_id for sample_id in sample_order if group_map[sample_id] == group]
    group_matrix = cf_logcpm.loc[group_samples, human_gene_order].to_numpy(dtype=np.float64)
    sample_similarity = cosine_similarity(group_matrix)
    np.fill_diagonal(sample_similarity, 0.0)

    group_dir = gmf_dir / group
    group_dir.mkdir(parents=True, exist_ok=True)
    start = time.time()
    np.random.seed(params["random_seed"])
    result_gene_cell, result_sample_gene, result_sample_cell = CellFreeGMF.GNMF(
        sample_cfRNA=group_matrix.copy(),
        cell_cfRNA=cell_matrix.copy(),
        sample_sim=sample_similarity.copy(),
        alpha=params["alpha"],
        beta=params["beta"],
        iter_num=params["iter_num"],
        random_seed=params["random_seed"],
    )
    elapsed = time.time() - start

    np.savez_compressed(
        group_dir / "gnmf_results.npz",
        Result_cfRNA_cell=result_gene_cell,
        Result_sample_cfRNA=result_sample_gene,
        Result_sample_cell=result_sample_cell,
        sample_name=np.asarray(group_samples),
        cfRNA_name=np.asarray(human_gene_order),
        cell_id=cell_meta["cell_id"].astype(str).to_numpy(),
        fine_annotation=cell_meta["fine_annotation"].astype(str).to_numpy(),
        cell_type=cell_meta["cell_type"].astype(str).to_numpy(),
    )
    pd.DataFrame(sample_similarity, index=group_samples, columns=group_samples).to_csv(
        group_dir / "sample_similarity.csv"
    )
    pd.DataFrame(result_sample_gene, index=group_samples, columns=human_gene_order).to_csv(
        group_dir / "result_sample_cfrna.csv"
    )

    group_scores = aggregate_scores(
        result_sample_cell,
        group_samples,
        cell_meta["cell_type"].astype(str).to_numpy(),
        cell_order,
    )
    group_scores.insert(1, "group", group)
    group_scores.to_csv(group_dir / "GMF_scores.csv", index=False)
    score_tables.append(group_scores)

    record = {
        "group": group,
        "n_samples": len(group_samples),
        "n_genes": len(human_gene_order),
        "n_cells": int(cell_matrix.shape[0]),
        "n_cell_types": 11,
        "elapsed_seconds": elapsed,
        **params,
    }
    (group_dir / "run_summary.json").write_text(
        json.dumps(record, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    run_records.append(record)

all_scores = pd.concat(score_tables, ignore_index=True)
all_scores["cell_type"] = pd.Categorical(all_scores["cell_type"], categories=cell_order, ordered=True)
all_scores = all_scores.sort_values(["cell_type", "group", "sample_id"])
all_scores.to_csv(output_dir / "GMF_scores_all.csv", index=False)
pd.DataFrame(run_records).to_csv(output_dir / "run_summary.tsv", sep="\t", index=False)

final_summary = {
    "input": input_summary,
    "groups": run_records,
    "score_file": str(output_dir / "GMF_scores_all.csv"),
    "CellFreeGMF": str(CellFreeGMF.__file__),
}
(output_dir / "summary.json").write_text(
    json.dumps(final_summary, indent=2, ensure_ascii=False), encoding="utf-8"
)
