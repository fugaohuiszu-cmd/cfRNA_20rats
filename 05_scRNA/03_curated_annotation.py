from pathlib import Path
import os
import pandas as pd

PROJECT = Path(os.environ.get("CFRNA_GMF_ROOT", ".")).resolve()
REF = (
    PROJECT
    / "CellFreeGMF_work"
    / "09_Set7_GSE169396_human_femoral_head_GMF"
    / "GSE169396_human_femoral_head_reference"
)
CURATED = REF / "curated_annotation"
CURATED.mkdir(parents=True, exist_ok=True)

TOP50 = REF / "tables" / "GSE169396_cluster_top50_markers.csv"
CELL_META = REF / "tables" / "GSE169396_FHTC_reconstructed_cell_annotation.csv"

if not TOP50.exists():
    raise FileNotFoundError(f"Missing top50 marker table: {TOP50}")
if not CELL_META.exists():
    raise FileNotFoundError(f"Missing cell metadata: {CELL_META}")

fine = {
    "0": "Mature neutrophils / granulocytes",
    "1": "Immature neutrophils / promyelocyte-like granulocytes",
    "2": "Activated inflammatory neutrophils / granulocytes",
    "3": "CD4-like / CCR7+ T cells",
    "4": "Cytotoxic T/NK-like cells",
    "5": "CD1C+ monocyte/DC-like myeloid cells",
    "6": "RAG1/RAG2+ pre-B / immature B cells",
    "7": "Plasma cells",
    "8": "Erythroid precursors / erythroblasts",
    "9": "Osteoblastic lineage / skeletal stromal-chondro-osteogenic cells",
    "10": "Mature erythrocytes / RBC",
    "11": "IFN-stimulated neutrophils / granulocytes",
    "12": "C1QC+ tissue-resident-like macrophages",
    "13": "Endothelial cells",
    "14": "Mast cells",
    "15": "Plasmacytoid dendritic cells / pDCs",
    "16": "DNTT+/VPREB1+ pro-B / early B-lineage cells",
}

major = {
    "0": "Granulocytes / neutrophils",
    "1": "Granulocytes / neutrophils",
    "2": "Granulocytes / neutrophils",
    "11": "Granulocytes / neutrophils",
    "3": "T cells / cytotoxic lymphocytes",
    "4": "T cells / cytotoxic lymphocytes",
    "5": "Monocyte/macrophage/DC-like myeloid",
    "12": "Monocyte/macrophage/DC-like myeloid",
    "6": "B-lineage cells",
    "16": "B-lineage cells",
    "7": "Plasma cells",
    "8": "Erythroid cells",
    "10": "Erythroid cells",
    "9": "Osteoblastic lineage cells",
    "13": "Endothelial cells",
    "14": "Mast cells",
    "15": "pDCs",
}

confidence = {cluster: "high" for cluster in fine}

markers = pd.read_csv(TOP50)
markers["cluster"] = markers["cluster"].astype(str)

rows = []
for cluster in sorted(fine, key=lambda x: int(x)):
    sub = markers[markers["cluster"].eq(cluster)].copy()
    if sub.empty:
        raise ValueError(f"No top50 markers found for cluster {cluster}")
    sub = sub.sort_values(["p_val_adj", "avg_log2FC"], ascending=[True, False])
    evidence = ";".join(sub["gene"].astype(str).head(10).tolist())
    rows.append(
        {
            "raw_cluster": cluster,
            "fine_annotation": fine[cluster],
            "GMF_major_celltype": major[cluster],
            "confidence": confidence[cluster],
            "top_marker_evidence": evidence,
        }
    )

cluster_table = pd.DataFrame(rows)
out_cluster = CURATED / "GSE169396_res0p3_curated_cluster_annotation_for_GMF.csv"
cluster_table.to_csv(out_cluster, index=False)

cell_meta = pd.read_csv(CELL_META)
if "cluster_reconstructed" not in cell_meta.columns:
    raise ValueError("cell metadata missing cluster_reconstructed column")
cell_meta["raw_cluster"] = cell_meta["cluster_reconstructed"].astype(str)
cell_meta = cell_meta.merge(cluster_table, on="raw_cluster", how="left")
if cell_meta["fine_annotation"].isna().any():
    missing = sorted(cell_meta.loc[cell_meta["fine_annotation"].isna(), "raw_cluster"].unique())
    raise ValueError(f"Missing curated labels for clusters: {missing}")

out_cell = CURATED / "GSE169396_res0p3_curated_cell_metadata_for_GMF.csv"
cell_meta.to_csv(out_cell, index=False)
